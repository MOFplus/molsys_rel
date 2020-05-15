"""

           findR

           find reactions in a ReaxFF trajectory

           this code is inspired by chemtrayzer but focuses on heterogeneous gas phase reactions
           2019 RUB Rochus Schmid



"""

import numpy as np
import pickle
import copy
import os

import molsys
from molsys.util import pdlpio2

from molsys.util.timing import timer, Timer
from molsys.util.print_progress import print_progress
from molsys import mpiobject
from molsys.util import rotations
from molsys.util import RDB

from graph_tool import Graph, GraphView
import graph_tool.topology as gtt
import graph_tool.util as gtu
import graph_tool.draw as gtd

"""
##################### NOTES /TO BE IMPLEMENTED

add xyz coord reference to frame object for mol object generation
-> not completely implemented

add pmol also to species
move molg generation to the frame and species objects

"""

#####################  FRAME CLASS ########################################################################################

class frame:
    """container to store info on a frame

    the main attributes are:
    - fid: frame id 
    - molg:    the molecular graph (edges below threshold are filtered but the bond order is in the edge property bord)
    - xyz:     numpy array (N,3) with the xyz coords of all atoms in this frame
    - pmol:    parent mol object ( as a reference)
    - specs:   dictionary with species (indexed by their name mid -> valid only in that frame!)

    """

    def __init__(self, fid, molg, xyz, pmol, specs=None):
        self.fid = fid
        self.molg = molg
        self.xyz = xyz
        self.pmol = pmol
        if specs:
            self.specs = specs
        else:
            self.specs = {}
        return

    def add_species(self, mid, tracked=True):
        # make all frame species (observed) with a local graph representation
        self.specs[mid] = species(self.fid, mid, self.molg, make_graph=True, tracked=tracked)
        return

    def make_species(self, mid):
        # return a species object (without a local graph)
        return species(self.fid, mid, self.molg, make_graph=False, tracked=False)

    @property
    def nspecies(self):
        return len(self.specs.keys())

    def plot(self, selection=False):
        if selection:
            print ("plotting species %s of frame %d" % (str(selection), self.fid))
            for s in selection:
                gfname = "frame_%d_species_%d.png" % (self.fid, s)
                g = self.specs[s].graph
                pos = gtd.arf_layout(g, max_iter=0)
                gtd.graph_draw(g, pos=pos, vertex_text=g.vp.aid, vertex_font_size=12, vertex_size=8, \
                                           edge_text=g.ep.bord, edge_font_size=10,\
                output_size=(800, 800), output=gfname, bg_color=[1,1,1,1])
        else:            
            print ("plotting all species of frame %d" % self.fid)
            for s in self.specs.keys():
                gfname = "frame_%d_species_%d.png" % (self.fid, s)
                g = self.specs[s].graph
                pos = gtd.arf_layout(g, max_iter=0)
                gtd.graph_draw(g, pos=pos, vertex_text=g.vp.aid, vertex_font_size=12, vertex_size=8, \
                output_size=(200, 200), output=gfname, bg_color=[1,1,1,1])
        return

    def get_border(self, b):
        """get the bond order for a bond b
        
        currently no test is made if the edge exists ... you need to check

        Args:
            b (tuple of indices): atom indices of bond
        """
        e = self.molg.edge(b[0], b[1])
        border = self.molg.ep.bord[e]
        return border

    def make_mol(self, species, xyz, pmol):
        """generate a single mol object from a list of species for a TS

        we get the current coordinates as xyz and the parent mol object pmol from the pdlp file

        NOTE: the additional species added to the frame are considered to be non tracked

        TBI: use frame obejcts pmol and xyz to do it -> refactor calling code
        
        Args:
            xyz (numpy): coordinates of the atoms in the frame
            pmol (molsys object): parent mol object from pdlp file 
        """
        aids = []
        for mid in species:
            if mid not in self.specs:
                self.add_species(mid, tracked=False)
            s = self.specs[mid]
            aids += list(s.aids)
        # aids.sort() # not sure if that is really needed
        self.mol = molsys.mol.from_array(xyz[aids])
        self.mol.set_cell(pmol.get_cell())
        self.mol.set_elems([pmol.elems[i] for i in aids])
        self.mol.set_real_mass()
        self.mol.center_com(check_periodic=False)
        self.mol.apply_pbc()
        self.mol.make_nonperiodic()
        # rotate into principal axes
        xyz = self.mol.get_xyz()
        self.mol.set_xyz(rotations.align_pax(xyz, masses = self.mol.get_mass()))
        # add connectivity
        conn = []
        for i in aids:
            v = self.molg.vertex(i)
            conn_i = []
            for j in v.all_neighbors():
                assert int(j) in aids
                conn_i.append(aids.index(int(j)))
            conn.append(conn_i)
        self.mol.set_conn(conn)
        return self.mol, aids

    def make_species_mol(self):
        """get a list of mol objects for all species in the frame

        NOTE: we use here the global xyz and pmol -> should be delegated to the species object which needs refactoring
        - we should attacg 
        """
        mols = {}
        for s in self.specs:
            sp = self.specs[s]
            sxyz = self.xyz[list(sp.aids)]
            mols[s] = self.specs[s].make_mol(sxyz, self.pmol)
        return mols

    ### some DEBUG methods ###

    def write_species(self):
        foldername = "frame_%d_species" % self.fid
        os.mkdir(foldername)
        os.chdir(foldername)
        mols = self.make_species_mol()
        for s in mols:
            m = mols[s]
            m.write("spec_%d.mfpx" % s)
        os.chdir("..")

    def get_main_species_formula(self):
        """get the sumformula of the tracked species
        """
        sp = list(self.specs.keys())
        sp.sort()
        sumforms = []
        for s in sp:
            aids = list(self.specs[s].aids)
            elems = [self.pmol.elems[i] for i in aids]
            cont  = list(set(elems))
            cont.sort()
            sumform = ""
            for e in cont:
                sumform += "%s%d " % (e, elems.count(e))
            sumforms.append(sumform)
        return sumforms

    def DEBUG_write_as_xyz(self, fname):
        f = open(fname, "w")
        f.write("%d\n\n" % self.pmol.natoms)
        for i in range(self.pmol.natoms):
            x, y, z = self.xyz[i]
            f.write("%3s %12.6f %12.6f %12.6f\n" % (self.pmol.elems[i], x, y, z))
        f.close()
        return

#####################  SPECIES CLASS ########################################################################################

class species:
    """container class to keep species info (per frame!)
    """

    def __init__(self, fid, mid, molg, make_graph=False, tracked = True):
        """init species
        
        Args:
            fid (int): frame number
            mid (int): molecule id ("name" of species from label components)
            molg (graph): molg of the frame
        """
        self.fid = fid
        self.mid = mid
        self.molg = molg # parent molgraph
        # find all vertices in molg that belong to this species mid
        vs = gtu.find_vertex(molg, molg.vp.mid, mid)
        self.aids = set([int(v) for v in vs]) # atomids -> TBI do we need to sort? seems like they are sroted as they come
        self.graph = None
        if make_graph:
            # now make a view of molg for this species         
            self.make_graph()
        self.tracked = tracked
        return

    def make_graph(self):
        vfilt = self.molg.new_vertex_property("bool")
        vfilt.a[:] = False
        for v in self.aids:
            vfilt[v] = True
        self.graph = GraphView(self.molg, vfilt=vfilt)
        return

    @property
    def natoms(self):
        return len(self.aids)

    def make_mol(self, xyz, pmol):
        """generate a mol object form the species

        we get the current coordinates as xyz and the parent mol object pmol from the pdlp file
        
        Args:
            xyz (numpy): coordinates of the atoms
            pmol (molsys object): parent mol object from pdlp file 
        """
        aids = list(self.aids) # in order to map backwards
        self.mol = molsys.mol.from_array(xyz)
        self.mol.set_cell(pmol.get_cell())
        self.mol.set_elems([pmol.elems[i] for i in aids])
        self.mol.set_real_mass()
        self.mol.center_com(check_periodic=False)
        self.mol.apply_pbc()
        self.mol.make_nonperiodic()
        # rotate into principal axes
        xyz = self.mol.get_xyz()
        self.mol.set_xyz(rotations.align_pax(xyz, masses = self.mol.get_mass()))
        # add connectivity
        if self.graph is None:
            self.make_graph()
        ctab = []
        for e in self.graph.edges():
            i = aids.index(int(e.source()))
            j = aids.index(int(e.target()))
            ctab.append([i, j])
        self.mol.set_ctab(ctab, conn_flag=True)
        return self.mol

#####################  FRAME COMPARE CLASS ########################################################################################

class fcompare:
    """This class compares two frames at different levels of resolution whether species are different
       All info collected during the compairson is kept in this class in order to be exploited later
    """

    def __init__(self, f1, f2):
        """generate comparer
        
        Args:
            f1 (frame object): first frame object
            f2 (frame object): second frame object (to be compared with f1)
        """
        self.f1 = f1
        self.f2 = f2
        # defaults
        self.compare_level = 0 # 0: no comparison, 1: atom id level, 2: connectivity level
        self.umatch_f1 = [s for s in self.f1.specs if self.f1.specs[s].tracked]
        self.umatch_f2 = [s for s in self.f2.specs if self.f2.specs[s].tracked]
        self.aids_match = []
        self.bond_match = []
        self.aids_analyzed = False
        self.bonds_analyzed = False
        self.reacs = []
        self.broken_bonds = []
        self.formed_bonds = []
        self.nreacs = 0 # nomber of independent reactions for this pair of frames .. should always be 1 (??)
        return

    def report(self, all = False):
        """this method is just to implement all levels of comparison

           it does not return anything and just reports ... meant for debugging
        """
        print("##################################################")
        print("FRAMES        %5d         %5d" % (self.f1.fid, self.f2.fid))
        mf = self.check_aids()
        if all:
            for m in self.aids_match:
                print("          species %5d == species %5d" % m)
                print(" natoms:          %5d            %5d" % (self.f1.specs[m[0]].natoms, self.f2.specs[m[1]].natoms))
        if mf > 0:
            self.analyse_aids()
            for r in self.reacs:
                print ("educts  : %s" % str(list(r[0].keys())))
                print ("products: %s" % str(list(r[1].keys())))
        return

    def check_aids(self, verbose = False):
        """check on level 1 for atom id matches

        this method implements a rather complex logic:
        which species form the educts and products to define a complete reaction event
        """
        if self.compare_level < 1:
            # we need to compare
            for sk1 in self.f1.specs:
                s1 = self.f1.specs[sk1]
                # find a corresponding species in frame 2
                for sk2 in self.umatch_f2:
                    s2 = self.f2.specs[sk2]
                    if s1.aids == s2.aids:      # NOTE: this works because the vertices are always properly sorted
                        # s1 and s2 match -> add to match and to remove
                        self.aids_match.append((sk1, sk2))
                        self.umatch_f1.remove(sk1)
                        self.umatch_f2.remove(sk2)
                        break
            # now matches are found -> set compare level
            self.compare_level = 1
        match_flag = 0 # all is equal on level 1
        if len(self.umatch_f1) > 0 or len(self.umatch_f2) > 0:
            match_flag = 1 # unmatched species!!
        if verbose:
            print ("species in f1   : %s" % str(self.f1.specs.keys()))
            print ("species in f2   : %s" % str(self.f2.specs.keys()))
            print ("unmatched in f1 : %s" % str(self.umatch_f1))
            print ("unmatched in f2 : %s" % str(self.umatch_f2))

        return match_flag

    def analyse_aids(self):
        if self.aids_analyzed:
            return
        if self.compare_level == 0:
            self.check_aids()
        if len(self.umatch_f1)==0 and len(self.umatch_f2)==0:
            # there is nothing to do
            return
        # we have soem unmatched species -> there is one (or more) reaction(s) between these frames
        # find groups of species that define a reaction
        #    all atom ids in the union of the educts sets must be also in the products set
        for sk1 in self.umatch_f1:
            # first search for atoms in the unmatched f2 species
            s1 = self.f1.specs[sk1]
            educts = {sk1: s1} # dictionary of species key/species for this reaction
            educt_aids = s1.aids.copy() # set of atomids
            products = {}
            product_aids = set()
            for sk2 in self.umatch_f2:
                s2 = self.f2.specs[sk2]
                common_aids = s1.aids & s2.aids
                if len(common_aids)>0:
                    products[sk2] = s2
                    product_aids |= s2.aids
            # do the following until educt_aids and product_aids match
            while not (educt_aids == product_aids):
                # which atoms are in products that are not in the educt species? --> add them
                for a in product_aids - educt_aids:
                    # to which species in frame 1 does this atom belong to?
                    esk = self.f1.molg.vp.mid[a]
                    # is this already in educts?
                    if esk not in educts:
                        # we need to make a new species object and add it (as non-tracked)
                        educts[esk] = self.f1.make_species(esk)
                        educt_aids |= educts[esk].aids
                # which atoms are in the educts that are not in the product species? add them
                for a in educt_aids - product_aids:
                    # to which species in frame 2 does this atom belong to?
                    psk = self.f2.molg.vp.mid[a]
                    # is this already in educts?
                    if psk not in products:
                        # we need to make a new species object and add it (as non-tracked)
                        products[psk] = self.f2.make_species(psk)
                        product_aids |= products[psk].aids
            # now add the final results to the reacs list
            self.reacs.append((educts, products))
        # the above will not work if there is no tracked species in umatch_f1 (but in umatch_f2)
        # ... in other words a species has "appeared" or formed by merging two or more untracked species
        if len(self.umatch_f1)==0:
            for sk2 in self.umatch_f2:
                s2 = self.f2.specs[sk2]
                products = {sk2: s2}
                product_aids = s2.aids.copy()
                # now find all the species in frame 1 to match the atoms 
                educts = {}
                educt_aids = set()
                for a in product_aids:
                    esk = self.f1.molg.vp.mid[a]
                    if esk not in educts:
                        educts[esk] = self.f1.make_species(esk)
                        educt_aids |= educts[esk].aids
                # do the following until educt_aids and product_aids match
                while not (educt_aids == product_aids):
                    # which atoms are in products that are not in the educt species? --> add them
                    for a in product_aids - educt_aids:
                        # to which species in frame 1 does this atom belong to?
                        esk = self.f1.molg.vp.mid[a]
                        # is this already in educts?
                        if esk not in educts:
                            # we need to make a new species object and add it (as non-tracked)
                            educts[esk] = self.f1.make_species(esk)
                            educt_aids |= educts[esk].aids
                    # which atoms are in the educts that are not in the product species? add them
                    for a in educt_aids - product_aids:
                        # to which species in frame 2 does this atom belong to?
                        psk = self.f2.molg.vp.mid[a]
                        # is this already in educts?
                        if psk not in products:
                            # we need to make a new species object and add it (as non-tracked)
                            products[psk] = self.f2.make_species(psk)
                            product_aids |= products[psk].aids
                # now add the final results to the reacs list
                self.reacs.append((educts, products))   
        self.nreacs = len(self.reacs)
        self.aids_analyzed = True
        return

    def find_react_bond(self):
        """find a reactive bond in a bimolecular reaction 
        
        Args:
            r (int): index in self.reac to analyse
        """
        assert self.aids_analyzed
        g1 = self.f1.molg
        g2 = self.f2.molg
        for r in range(self.nreacs):
            broken_bonds = []
            formed_bonds = []
            educts, products = self.reacs[r]
            # get all involved atoms
            aids = set()
            for s in educts:
                aids |= educts[s].aids
            for a in aids:
                bonds1 = []
                v1 = g1.vertex(a)
                for e in v1.out_edges(): 
                    bonds1.append(int(e.target()))
                bonds2 = []
                v2 = g2.vertex(a)
                for e in v2.out_edges():
                    bonds2.append(int(e.target()))
                bs = set(bonds1) - set(bonds2)
                fs = set(bonds2) - set(bonds1)
                for b in bs:
                    if b > a:
                        broken_bonds.append((a,b))
                for b in fs:
                    if b > a:
                        formed_bonds.append((a,b))
            self.broken_bonds.append(broken_bonds)
            self.formed_bonds.append(formed_bonds)
        return

    def check_bonds(self, verbose = False):
        """check on level 2 for identical bonds

        we use the existing pairs of species in self.aids_match to identify species that have identical aids
        => now we test if they have the same bonding pattern
        """
        assert self.compare_level > 0
        self.missmatch = []
        for p in self.aids_match:
            # get the species of the pair
            s1 = self.f1.specs[p[0]]
            s2 = self.f2.specs[p[1]]
            # TBI: this might not be enough 
            #      do we need elements as vertex properties to be considered?
            #      what about a tautomerism when the initial and fianl state are symmetric?
            f = gtt.isomorphism(s1.graph, s2.graph)
            # print ("%d %d %d %d %s" % (self.f1.fid, self.f2.fid, p[0], p[1], f))
            if f == True:
                self.bond_match.append(p)
            else:
                self.missmatch.append(p)
        self.compare_level = 2
        if len(self.missmatch) > 0:
            # unmatched species on level 2
            return 2
        else:
            return 0 # all equal

    def analyse_bonds(self):
        if self.bonds_analyzed:
            return
        if self.compare_level < 2:
            self.check_bonds()
        if len(self.missmatch) == 0:
            return
        # now analyse all the pairs in self.missmatch: they have identical aids but a diffrent bond graph
        self.nreacs = len(self.missmatch)
        for p in self.missmatch:
            # get the species of the pair
            s1 = self.f1.specs[p[0]]
            s2 = self.f2.specs[p[1]]
            self.reacs.append(({p[0]: s1},{p[1]: s2}))
            # get all edges as vertex id tuples (int tuples)
            bs1 = set([(int(e.source()),int(e.target())) for e in s1.graph.edges()])
            bs2 = set([(int(e.source()),int(e.target())) for e in s2.graph.edges()])
            self.broken_bonds.append(list(bs1-bs2)) # broken: bond in frame1 but not in frame2
            self.formed_bonds.append(list(bs2-bs1)) # formed: bond in frame2 but not in frame1
        self.bonds_analyzed = True
        return


    def check(self, verbose = False):
        """check identity of species on all levels (1 and 2) 
        """
        mf = self.check_aids(verbose=verbose)
        if mf > 0:
            return mf
        # aids are equal -> chek bonds
        return self.check_bonds(verbose=verbose)

#####################  REACTIVE EVENT CLASS ########################################################################################
# this class is instantiated with a comparer (between two frames)
# and stores a reactive event 
# it knows the TS_fid and if the event is bi- or unimolecular

class revent:

    def __init__(self, comparer, fR, unimol= False):
        """generate a reaction event object

        TBI: what to do if there are more then one reaction event (in two diffrent tracked species)
             at the same time ... this is properly tracked in the comparer (nreacs>1)
             but this means there are more revents to be generated. should we do this recursively?
             how to store?
        
        Args:
            comparer (fcompare object): comparer that gave a reactive event
            fR (parent findR object): to access process_frame
            unimol (bool, optional): is unimolecular. Defaults to False.
        """
        self.unimol = unimol
        self.fR = fR
        self.comparer = comparer
        # currently we allow only for single reaction events per frame 
        # this would change if there is more than one tracked species ....
        assert comparer.nreacs == 1 , "Currently only single reaction events are processed"
        r = 0 # pick first reaction (the only one)
        educts, products = comparer.reacs[r]
        self.broken_bonds     = comparer.broken_bonds[r]
        self.formed_bonds     = comparer.formed_bonds[r]
        f1               = comparer.f1
        f2               = comparer.f2
        # choose which frame we use as TS
        # everything is referenced to f1 which is 0 (f2 is +1)
        # find avereage bond order of reactive bonds at f1/f2
        f1_averborder = 0.0
        if len(self.broken_bonds) >0:
            for b in self.broken_bonds:
                f1_averborder += f1.molg.ep.bord[f1.molg.edge(b[0], b[1])]
            f1_averborder/len(self.broken_bonds)
        f2_averborder = 0.0
        if len(self.formed_bonds) >0:
            for b in self.formed_bonds:
                f2_averborder += f2.molg.ep.bord[f2.molg.edge(b[0], b[1])]
            f2_averborder/len(self.formed_bonds)
        if f1_averborder == 0.0:
            TS_rfid = 1
        elif f2_averborder == 0.0:
            TS_rfid = 0
        else:
            if abs(f1_averborder-0.5) < abs(f2_averborder-0.5):
                # f1 closer to TS
                TS_rfid = 0
            else:
                TS_rfid = 1
        # now store data depending on relative frame id (rfid) of the TS
        if TS_rfid == 0:
            self.TS = f1
            self.PR = f2
            self.ED = self.fR.process_frame(self.TS.fid-1)
            # get corresponding species numbers
            self.TS_spec = educts
            self.PR_spec = products
            loccomp = fcompare(self.ED, self.TS)
            if loccomp.check_aids() == 0:
                # no change in atom ids .. we can use TS species for ED as well
                self.ED_spec = educts
            else:
                print ("Houston we have a problem!!! species changed between ED and TS")
        else:
            self.ED = f1
            self.TS = f2
            self.PR = self.fR.process_frame(self.TS.fid+1)
            self.ED_spec = educts
            self.TS_spec = products
            loccomp = fcompare(self.TS, self.PR)
            if loccomp.check_aids() == 0:
                # no change in atom ids .. we can use TS species for PR as well
                self.PR_spec = products
            else:
                print ("Houston we have a problem!!! species changed between TS and PR")
        # get TS_fid for ease
        self.TS_fid = self.TS.fid
        return


#####################  CENTRAL FINDR CLASS ########################################################################################

class findR(mpiobject):

    def __init__(self, pdlpfilename,  stage, rdb_path, mpi_comm = None, out = None):
        super(findR,self).__init__(mpi_comm, out)
        # To be tested: open read only on all nodes in parallel .. does it work?
        self.pdlp = pdlpio2.pdlpio2(pdlpfilename, filemode="r")
        self.mol = self.pdlp.get_mol_from_system()
        assert stage in self.pdlp.get_stages(), "Stage %s not in stages in file"
        self.traj = self.pdlp.h5file["/%s/traj" % stage]
        self.rest = self.pdlp.h5file["/%s/restart" % stage]
        data = list(self.traj.keys())
        # make sure that xyz and bondord/bondterm is available
        assert "xyz" in data
        assert "bondord" in data
        assert "bondtab" in data
        self.f_xyz = self.traj["xyz"]
        self.f_bondord = self.traj["bondord"]
        self.f_bondtab = self.traj["bondtab"]
        self.cell = np.array(self.rest["cell"])
        # get some basic info from the mol object
        self.natoms = self.mol.get_natoms()
        self.elems  = self.mol.get_elems()
        self.nframes = self.f_xyz.shape[0]
        self.maxbond = self.f_bondord.shape[1]
        # defaults and data structures
        self.frames = self.nframes*[None]
        self.timer = Timer(name="findR")
        self.bondord_cutoff = 0.5  # below this bond roder a bond is registered but swithced off
        self.min_atom = 6          # below this number of atoms a species is considered as gasphase
        self.min_elems = {"c" : 2} # below this number of atoms of the given element in the species it is not tracked
        # for debugging
        self.nunconnected = 0
        self.unconnected = []
        # for the search
        self.sstep ={
            "forward" : 200,
            "back"    : -10,
            "fine"    : 1,
        }
        self.skip_recross = 2 # number of frames to check if a reaction event was reverted
        # now open the RDB    TBI: how to deal with a parallel run?
        self.rdb = RDB.RDB(rdb_path, do_commit=False) # we open in no-commit mode -> have to call commit ourselves
        # TBI  get temp and timestep from pdlp file
        self.rdb.set_md_run(pdlpfilename, stage, nframes=self.nframes, temp=2000.0, timestep=10.0)
        return

    @timer("process_frame")
    def process_frame(self, fid):
        """process a single frame

        Args:
            fid (int): number of frame to process
            
        Returns:
            frame: frame object
        """
        if self.frames[fid] is not None:
            return self.frames[fid]
        bondord = np.array(self.f_bondord[fid])
        bondtab = np.array(self.f_bondtab[fid])
        molg = Graph(directed=False)
        molg.add_vertex(n=self.natoms)
        molg.vp.aid  = molg.new_vertex_property("int")
        molg.vp.aid.a[:] = np.arange(self.natoms, dtype="int32")
        molg.vp.filt = molg.new_vertex_property("bool")
        molg.vp.mid  = molg.new_vertex_property("int")
        molg.ep.bord = molg.new_edge_property("float")
        molg.ep.filt = molg.new_edge_property("bool")
        for j in range(self.maxbond):
            o = bondord[j]
            # invalid entries are marked with a bondord = -1.0
            if o < 0.0:
                break
            e = bondtab[j]-1
            # TO BE REMOVED ... this is a legacy check for the old incorrect pdlp files 
            if (e[0] < 0 or e[0] >= self.natoms or e[1] < 0 or e[1] >= self.natoms):
                break
            # END TO BE REMOVED
            newb = molg.add_edge(e[0], e[1])
            isbond = o > self.bondord_cutoff
            molg.ep.bord[newb] = o
            molg.ep.filt[newb] = isbond 
        # apply edge filter
        molg.set_edge_filter(molg.ep.filt)
        mid, hist = gtt.label_components(molg, vprop=molg.vp.mid)
        nspecies_all = len(hist)
        nspecies = []
        # collect number of critical elements in the species
        elem_spec = {}
        for e in self.min_elems.keys():
            spec_ecount = np.zeros(nspecies_all, dtype="int32")
            for i in range(self.natoms):
                if self.elems[i] == e:
                    spec_ecount[mid[i]] += 1
            elem_spec[e] = spec_ecount
        for i in range(nspecies_all):
            track = True
            if hist[i] < self.min_atom:
                track = False
            for e in self.min_elems.keys():
                if elem_spec[e][i] < self.min_elems[e]:
                    track = False
            if track:
                nspecies.append(i)
            # print (nspecies)
        # TBI shall we reset the vertex filter of molg again?
        # add xyz and parent mol
        xyz = np.array(self.f_xyz[fid])      
        f = frame(fid, molg, xyz, self.mol)
        for mid in nspecies:
            f.add_species(mid)
        # store the frame
        self.frames[fid] = f
        return f

    def get_comparer(self, f1, f2):
        if not self.frames[f1]:
            self.process_frame(f1)
        if not self.frames[f2]:
            self.process_frame(f2)
        return fcompare(self.frames[f1], self.frames[f2])
    
    @timer("do_frames")
    def do_frames(self, start=0, stop=None, stride=1, progress=True, plot=False):
        """generate all frames (molg, species) representations
        
        Args:
            start (integer, optional): starting frame. Defaults to None.
            stop (integer, optional): final frame (not inlcuded). Defaults to None.
            stride (integer, optional): stride. Defaults to None.
        """ 
        if not stop:
            stop = self.nframes
        if progress:
            self.pprint ("Processing frames from %d to %d (stride=%d)" % (start, stop, stride))
        for i in range(start, stop, stride):
            if progress:
                ci = (i-start)/stride
                print_progress(ci, (stop-start)/stride, prefix="Processing frames:", suffix=" done!")
            f = self.process_frame(i)
            if plot:
                f.plot()
        return

    @timer("search_frames")
    def search_frames(self, verbose=False):
        """search the frames for reactions 

        This is the core routine of findR
        """
        mode = "forward"     # search mode
        currentf = 0          # current frame
        last_event = currentf # last reaction event => product frame of a reaction (initial frame is also an event)
        open_events = []      # list of reaction events located in a forward segment (located by backtracking)
        self.process_frame(currentf)
        self.store_initial_event()
        nextf = currentf+self.sstep[mode]
        # start mainloop (just search on until we are at the end of the file)
        # init some variables
        delta_segment_end = 0
        segment_events = []
        first_event = None # this variable stores the maximum fid of a reaction event at the beginning of a segment being a recrossing
        last_event  = None # this stores the last revent object of a segment (None in the first segment) 
        stop = False
        while nextf < self.nframes and nextf >= 0 and not stop:
            if not verbose:
                print_progress(nextf/self.sstep["forward"], self.nframes/self.sstep["forward"], suffix="Scanning Frames")
            self.process_frame(nextf)
            # do comparison and report
            comparer = self.get_comparer(currentf, nextf)
            flag = comparer.check()
            if flag>0:
                # a difference was found between current and next
                if mode == "forward":
                    segment_end = nextf
                    segment_start = currentf
                    mode = "back"
                elif mode == "back":
                    subseg_end = currentf  # NOTE for the subsegment we go backwards and the role of nextf
                    subseg_start = nextf   #      and currentf changes
                    mode = "fine"
                elif mode == "fine":
                    if flag == 1:
                        # we found a bimolecular reaction
                        comparer.analyse_aids()
                        # search critical bond
                        comparer.find_react_bond()   
                        # now we are ready to make a reaction event and store it
                        revt = revent(comparer, self)
                        TS_fid = revt.TS_fid
                        if verbose:
                            print ("###########  Event at %d  #############" % TS_fid)
                    elif flag == 2:
                        # unimolecular reaction
                        comparer.analyse_bonds() # this method already finds the reactive bonds
                        # make reaction event and store
                        revt = revent(comparer, self, unimol=True)
                        TS_fid = revt.TS_fid
                        if verbose:
                            print ("###########  Unimol Event at %d  #############" % TS_fid)
                    # now add the event to the list
                    segment_events.append(revt)
                    # first we need to make sure that there is no other reaction in the subsegment
                    #   test if the product (TS_fid+1) is equal to the end of the subsegment (subseg_end)
                    #   if this is not the case then we have to go forward in fine mode further on
                    comparer_subseg = self.get_comparer(TS_fid+1, subseg_end)
                    flag = comparer_subseg.check()
                    if flag > 0:
                        if verbose:
                            print ("The product frame %d and the end of the subsegment at %d are not equal" % (TS_fid+1, subseg_end))
                            print ("should continue forward in fine mode ... ")
                        # need not to do anything (mode stay fine)
                    else:
                        # next we need to make sure that the start of the subseg species
                        #   are really idential to where we started the segment (segment_start)
                        comparer_segment = self.get_comparer(segment_start, subseg_start)
                        flag = comparer_segment.check()
                        if flag > 0:
                            if verbose:
                                print ("The subseg start %d and the segment start %d are not equal" % (subseg_start, segment_start))
                                print ("continue in back mode")
                            nextf = subseg_start
                            mode = "back"
                        else:
                            # At this point we can be sure that all events in the segment are now in the list
                            # now we should get rid of recrossings, add to the DB and connect
                            # first sort our events (because of backtracking they might not be sorted)
                            event_TSfids = [r.TS_fid for r in segment_events]
                            event_index = np.argsort(event_TSfids)
                            segment_events[:] = [segment_events[i] for i in event_index]
                            event_TSfids = [r.TS_fid for r in segment_events]
                            if verbose:
                                print ("events in segment: %s" % str(event_TSfids))
                            ### handle recrossings ################################################
                            if self.skip_recross > 0:
                                # find recrossings if skip_recross is larger than 1 -> add event indices in list recross and remove
                                recross = []
                                start_recross = 0  # if the first event is part of a recrossing over the segment boundary then set this to 1
                                # is there a possible early event?
                                if first_event is not None:
                                    # we should have a reaction event before the frame given in first event
                                    if event_TSfids[0] <= first_event:
                                        if verbose:
                                            print ("  found an early event for recrossing at %d" % event_TSfids[0])
                                        first_event = None
                                        recross = [0]
                                        start_recross = 1
                                # check the events for any recross
                                if (len(event_TSfids)+start_recross) > 1: 
                                    for i in range(start_recross, len(event_TSfids)-1):
                                        e1 = event_TSfids[i]
                                        e2 = event_TSfids[i+1]
                                        if (e2-e1) <= self.skip_recross:
                                            # two events are in recrossing distance 
                                            if verbose:
                                                print ("possible recross between %d and %d" % (e1, e2))
                                            compare_recross = self.get_comparer(e1-1, e2+1) # compare ED of event1 and PR of event2
                                            flag = compare_recross.check()
                                            if flag == 0:
                                                # this is a recrossing
                                                recross += [i, i+1]
                                                if verbose:
                                                    print ("recrossing found")
                                            else:
                                                if verbose:
                                                    print ("no recrossing")
                                # check if the last event is close enough to the segment bound 
                                e = event_TSfids[-1]
                                if segment_end-e <= self.skip_recross:
                                    compare_recross = self.get_comparer(e-1, e+self.skip_recross+1)
                                    flag = compare_recross.check(verbose=True)
                                    if flag == 0:
                                        # this is a recrossing over the segment bound
                                        if verbose:
                                            print ("recrossing over segment boundary from %d to %d" % (e-1, e+self.skip_recross+1))
                                            # DEBUG DEBUG ... analysie what happens over the segment boundary
                                            for ie in range(e-3, e+self.skip_recross+4):
                                                if self.frames[ie] is None:
                                                    self.process_frame(ie)
                                                print ("Sumform at frame %5d : %s" % (ie, self.frames[ie].get_main_species_formula()))
                                            # DEBUG DEBUG END
                                        recross += [len(event_TSfids)-1]
                                        first_event = e+self.skip_recross
                                        # set segment_end to the end of the recrossing .. if we move forward we want to compare to the end of the recrossing event
                                        print ("segment end was: %d" % segment_end)
                                        segment_end_new = e+self.skip_recross+1
                                        delta_segment_end = segment_end_new - segment_end
                                        print ("delta is %d" % delta_segment_end)
                                        segment_end = segment_end_new
                                if verbose:
                                    if len(recross) > 0:
                                        print ("The following event indices are marked to be removed because of recrossing")
                                        print (str(recross))
                                # now recrossing events are complete and we can remove them
                                if len(recross) > 0:
                                    segment_events = [segment_events[e] for e in range(len(segment_events)) if e not in recross]
                                    event_TSfids = [r.TS_fid for r in segment_events]
                                    if verbose:
                                        print ("remaining segments to be stored: %s" % str(event_TSfids))
                            ###### end of handle recrossings ####################################
                            # store events in DB now
                            for revt in segment_events:
                                if revt.unimol:
                                    self.store_uni_event(revt)
                                else:
                                    self.store_bim_event(revt)
                            # connect events (start with connecting last_event with the first in the segment)
                            if len(segment_events) > 0:
                                for revt in segment_events:
                                    self.connect_events(last_event, revt, verbose=verbose)
                                    last_event = revt
                            # now move forward again stating from segment_end
                            mode = "forward"
                            nextf = segment_end
                            # clear segment_events
                            segment_events = []
                            if verbose:
                                print ("#########################################################")
                                print ("all events in segment stored and connected .. moving forward again"   )
            currentf = nextf
            nextf = currentf+self.sstep[mode]-delta_segment_end
            delta_segment_end = 0
            # capture problems -- this should not happen
            if mode == "back":
                if nextf < segment_start:
                    back_delt = segment_start - nextf
                    if back_delt > self.skip_recross+1:
                        print ("segement_start is %d" % segment_start)
                        print ("nextf is %d" % nextf)
                        print ("backtracking went wrong -> abort")
                        raise
                    else:
                        if verbose:
                            print ("backtracking to a shifted segment bound (due to recross)")
                            print ("We shift nextf from %d to %d" % (nextf, segment_start))
                        nextf = segment_start
            if mode == "fine":
                if nextf > subseg_end:
                    print ("dinetracking went wrong -> abort")
                    raise 
            if verbose:
                sumform = self.frames[currentf].get_main_species_formula()
                print ("&& current %10d  next %10d  mode %8s   (current frame sumformula %s" % (currentf, nextf, mode, sumform))
            # check if the currentf has zero tracked species -> then stop searching
            if len(self.frames[currentf].specs) == 0:
                print ("Current frame %d has zero species to track ... stop searching")
                stop = True
        # end of mainloop
        self.rdb.commit()
        # DEBUG
        if self.nunconnected > 0:
            print ("NUMBER OF UNCONNECTED EVENTS: %d" % self.nunconnected)
            for r in self.unconnected:
                print ("%5d --> %5d" % r)
        return

    def store_initial_event(self):
        """for consistency we need an event for frame 0 just registering the initial tracked species as products

        Note: we do not store any ED or TS species. 
        """
        f = self.frames[0]
        spec = list(f.specs.keys()) # get all tracked species
        # store reaction event
        revID = self.rdb.register_revent(
            -1,          # frame ID of TS is -1 (because 0 is the PR!)
            [],          # no ED
            [],          # no TS
            spec,        # spec IDs of "products"
            0,           # no educts
            len(spec),   # tracked species (all)
            []           # no rbonds
        )
        # make mol objects and store
        xyz = np.array(self.f_xyz[0])
        for s_id in spec:
            s = f.specs[s_id]
            m = s.make_mol(xyz[list(s.aids)], self.mol)
            self.rdb.add_md_species(
                revID,
                m,
                s_id,
                1,            # PR
                tracked=True
            )
        return

    @timer("store_bim_event")
    def store_bim_event(self, revt):
        """store species and further info in database

        Args:
            revt (revent object): reaction event object
        """
        # now we should have all data
        # generate mol objectes with coordinates
        xyz_ed = np.array(self.f_xyz[revt.ED.fid])
        xyz_ts = np.array(self.f_xyz[revt.TS.fid])
        xyz_pr = np.array(self.f_xyz[revt.PR.fid])
        ED_mol = []
        ED_spec_tracked = []
        # ED/PR_spec are dicitionaries of species .. we need a sorted list of integers
        ED_spec_id = list(revt.ED_spec.keys())
        ED_spec_id.sort()
        for s in ED_spec_id:
            aids = list(revt.ED_spec[s].aids)
            m = revt.ED_spec[s].make_mol(xyz_ed[aids], self.mol)
            ED_mol.append(m)
            if revt.ED_spec[s].tracked:
                ED_spec_tracked.append(s)
        PR_mol = []
        PR_spec_tracked = []
        PR_spec_id = list(revt.PR_spec.keys())
        PR_spec_id.sort()
        for s in PR_spec_id:
            aids = list(revt.PR_spec[s].aids)
            m = revt.PR_spec[s].make_mol(xyz_pr[aids], self.mol)
            PR_mol.append(m)
            if revt.PR_spec[s].tracked:
                PR_spec_tracked.append(s)
        # only one TS (all species) .. use make_mol of frame
        TS_spec_id = list(revt.TS_spec.keys())
        TS_spec_id.sort()
        TS_mol, TS_aids = revt.TS.make_mol(revt.TS_spec, xyz_ts, self.mol)
        # map the broken/formed bonds to atom ids of the TS subsystem
        rbonds_global = revt.formed_bonds+revt.broken_bonds
        rbonds = []
        for b in rbonds_global:
            rbonds.append(TS_aids.index(b[0]))
            rbonds.append(TS_aids.index(b[1]))
        # now let us register the reaction event in the database
        revID = self.rdb.register_revent(
            revt.TS.fid,
            ED_spec_id,
            TS_spec_id,
            PR_spec_id,
            len(ED_spec_tracked),
            len(PR_spec_tracked),
            rbonds
        )
        # add the md species (mol objects etc) to this entry
        for i,s in enumerate(ED_spec_id):
            self.rdb.add_md_species(
                revID,
                ED_mol[i],
                s,
                -1,           # ED
                tracked= s in ED_spec_tracked
            )
        for i,s in enumerate(PR_spec_id):
            self.rdb.add_md_species(
                revID,
                PR_mol[i],
                s,
                1,            # PR
                tracked= s in PR_spec_tracked
            )
        self.rdb.add_md_species(
            revID,
            TS_mol,
            TS_spec_id[0],    # species of TS are stored in revent, here only first
            0
        )
        return

    @timer("store_uni_event")
    def store_uni_event(self, revt):
        """store species and further info in database for an unimol reaction event
        (one ED/TS and PR species only)

        Args:
            revt (revent object): reaction event object
        """
        assert revt.unimol==True, "store_uni_event called for a bimolecular event!"
        # now we should have all data
        # generate mol objectes with coordinates
        xyz_ed = np.array(self.f_xyz[revt.ED.fid])
        xyz_ts = np.array(self.f_xyz[revt.TS.fid])
        xyz_pr = np.array(self.f_xyz[revt.PR.fid])
        # there must be only one species in ED, PR and TS
        # ED
        assert len(revt.ED_spec) == 1
        ED_spec_id = next(iter(revt.ED_spec))
        ED_spec = revt.ED_spec[ED_spec_id]
        ED_mol = ED_spec.make_mol(xyz_ed[list(ED_spec.aids)], self.mol)
        # TS
        assert len(revt.TS_spec) == 1
        TS_spec_id = next(iter(revt.TS_spec))
        TS_spec = revt.TS_spec[TS_spec_id]
        TS_mol = TS_spec.make_mol(xyz_ts[list(TS_spec.aids)], self.mol)
        #PR
        assert len(revt.PR_spec) == 1
        PR_spec_id = next(iter(revt.PR_spec))
        PR_spec = revt.PR_spec[PR_spec_id]
        PR_mol = PR_spec.make_mol(xyz_pr[list(PR_spec.aids)], self.mol)
        # map the broken/formed bonds to atom ids of the TS subsystem
        rbonds_global = revt.formed_bonds+revt.broken_bonds
        rbonds = []
        TS_aids = list(TS_spec.aids)
        for b in rbonds_global:
            if (b[0] not in TS_aids) or (b[1] not in TS_aids):
                import pdb; pdb.set_trace() 
            rbonds.append(TS_aids.index(b[0]))
            rbonds.append(TS_aids.index(b[1]))
        # now let us register the unimol reaction event in the database
        revID = self.rdb.register_revent(
            revt.TS.fid,
            [ED_spec_id],
            [TS_spec_id],
            [PR_spec_id],
            1, # number of tracked ed species .. for uni always 1
            1, # number of tracked pr species .. for uni always 1
            rbonds,
            uni = True
        )
        # add the md species (mol objects etc) to this entry
        self.rdb.add_md_species(
            revID,
            ED_mol,
            ED_spec_id,
            -1,   # ED
        )
        self.rdb.add_md_species(
            revID,
            PR_mol,
            PR_spec_id,
            1,            # PR
        )
        self.rdb.add_md_species(
            revID,
            TS_mol,
            TS_spec_id,  
            0
        )
        return

    @timer("connect_events")
    def connect_events(self, revt1, revt2, verbose=False):
        """connecting two reaction events 
        
        Args:
            revt1 (revent object): first reaction event (can be None to indicate initial frame)
            revt2 (revent obejct): second reaction event
        """
        if revt1 is None:
            fid1 = 0
        else:
            fid1 = revt1.TS_fid+1
        fid2 = revt2.TS_fid-1
        comparer = self.get_comparer(fid1, fid2)
        flag = comparer.check()
        # if the two events are side by side (revt2.TS_fid = revt1._TSfid+1) then fid1 > fid2
        #  => this will always lead to a fail ... this should never happen.
        # in this case we can not connect from PR1 to ED2. instead ED2 is TS1!
        if (fid1 >= fid2) or (flag != 0):
            if verbose:
                print ("##########################################################")
                print ("connecting between %d and %d" % (fid1, fid2))
                print ("products of revent 1: %s" % str(revt1.PR_spec))
                print ("educts of revent   2: %s" % str(revt2.ED_spec))
                print ("sumforms frame 1 %s" % str(self.frames[fid1].get_main_species_formula()))
                print ("sumforms frame 2 %s" % str(self.frames[fid2].get_main_species_formula()))
                print ("##########################################################")
        # assert flag == 0, "ED of revt2 not equal PR revt1 .. this should never happen"
        if flag != 0:
            print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print ("This should never happen!! Could not Connect!!")
            print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            import pdb; pdb.set_trace()

            self.nunconnected += 1
            self.unconnected.append((fid1, fid2))
            return # DEBUG DEBUG -- just ignore connection in this case
        if len(comparer.bond_match) == 0:
            print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print ("This is unexpected!!")
            print ("the comparer did not find any matching species")
            import pdb; pdb.set_trace()
            self.nunconnected += 1
            self.unconnected.append((fid1, fid2))

        for match in comparer.bond_match:
            self.rdb.set_react(fid1, fid2, match[0], match[1])
        return

    ##########################  DEBUG stuff #################################################################

    def deb_check_every(self, every):
        """run over the frames and check for changes between every n frame
        
        Args:
            every (int): stride to check frames
        """
        self.process_frame(0)
        for i in range(every, self.nframes, every):
            self.process_frame(i)
            comp = self.get_comparer(i-every, i)
            flag = comp.check()
            taetae = ""
            if flag > 0:
                taetae = "!!!!!!!!"
            print("flag between %4d and %4d : %1d  %s" % (i-every, i, flag, taetae))
        return


    ########################## in case of restart ###########################################################

    @timer("store_frames")
    def store_frames(self, fname):
        f = open(fname, "wb")
        pickle.dump(self.frames, f)
        f.close()
        return

    @timer("load_frames")
    def load_frames(self, fname):
        f = open(fname, "rb")
        self.frames = pickle.load(f)
        f.close()
        return

    def report(self):
        self.timer.write()










