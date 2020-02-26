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

#####################  FRAME CLASS ########################################################################################

class frame:
    """container to store info on a frame

    the main attributes are:
    - fid: frame id 
    - molg:    the molecular graph (edges below threshold are filtered but the bond order is in the edge property bord)
    - specs:   dictionary with species (indexed by their name mid -> valid only in that frame!)
    """

    def __init__(self, fid, molg, specs=None):
        self.fid = fid
        self.molg = molg
        if specs:
            self.specs = specs
        else:
            self.specs = {}
        return

    def add_species(self, mid):
        # make all frame species (observed) witha local graph representation
        self.specs[mid] = species(self.fid, mid, self.molg, make_graph=True)
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
        """generate a mol object from a list of species

        we get the current coordinates as xyz and the parent mol object pmol from the pdlp file
        
        Args:
            xyz (numpy): coordinates of the atoms in the frame
            pmol (molsys object): parent mol object from pdlp file 
        """
        aids = []
        for mid in species:
            if mid not in self.specs:
                self.add_species(mid)
            s = self.specs[mid]
            aids += list(s.aids)
        aids.sort() # not sure if that is really needed
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
        aids.sort() # not sure if that is really needed
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
    """This class compares two frames at different levels of reolution whether species are different
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
        self.umatch_f1 = list(self.f1.specs.keys())
        self.umatch_f2 = list(self.f2.specs.keys())
        self.aids_match = []
        self.bond_match = []
        self.aids_analyzed = False
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

    def check_aids(self):
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
        return match_flag

    def analyse_aids(self):
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

    def check_bonds(self):
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
            if True:
                self.bond_match.append(p)
            else:
                self.missmatch.append(p)
        self.compare_level = 2
        if len(self.missmatch) > 0:
            # unmatched species on level 2
            return 2
        else:
            return 0 # all equal

    def check(self):
        """check identity of species on all levels (1 and 2) 
        """
        mf = self.check_aids()
        if mf > 0:
            return mf
        # aids are equal -> chek bonds
        return self.check_bonds()

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
        self.min_atom = 10         # below this number of atoms a species is considered as gasphase
        # for the search
        self.sstep ={
            "forward" : 200,
            "back"    : -10,
            "fine"    : 1,
        }
        # now open the RDB    TBI: how to deal with a paralle run?
        self.rdb = RDB.RDB(rdb_path)
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
            e = bondtab[j]-1
            o = bondord[j]
            # WARNING: this is not safe ... we assume that veriex numbers are =>0 and <natoms tp be valid ... this is likely but not entirely safe
            #                     TBI: adda table with number of bonds in the pdlp write out and use that for the loop
            if (e[0] >= 0 and e[0]<self.natoms and e[1] >= 0 and e[1] < self.natoms):
                newb = molg.add_edge(e[0], e[1])
                isbond = o > self.bondord_cutoff
                molg.ep.bord[newb] = o
                molg.ep.filt[newb] = isbond 
            else:
                # invalid vertices in e -> stop (this is NOT safe!! see comment above)
                break
        # apply edge filter
        molg.set_edge_filter(molg.ep.filt)
        mid, hist = gtt.label_components(molg, vprop=molg.vp.mid)
        nspecies_all = len(hist)
        nspecies = []
        # TBI: currently unmonitored species are filtered only on size .. there could be a per element filtering etc.
        for i in range(nspecies_all):
            if hist[i] >= self.min_atom:
                nspecies.append(i)
        # TBI shall we reset the vertex filter of molg again?      
        f = frame(fid, molg)
        for mid in nspecies:
            f.add_species(mid)
        # store the frame
        self.frames[fid] = f
        return f

    def get_comparer(self, f1, f2):
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
        while nextf < self.nframes and nextf > 0:
            if not verbose:
                print_progress(nextf/self.sstep["forward"], self.nframes/self.sstep["forward"], suffix="Scanning Frames")
            self.process_frame(nextf)
            # do comparison and report
            comparer = self.get_comparer(currentf, nextf)
            flag = comparer.check()
            if flag>0:
                # a difference was found between current and next
                if mode == "forward":
                    farthestf = nextf
                    mode = "back"
                elif mode == "back":
                    back_curr = nextf
                    mode = "fine"
                elif mode == "fine":
                    if flag == 1:
                        # bimolecular reaction
                        comparer.analyse_aids()
                        # search critical bond
                        comparer.find_react_bond()   
                        # now store the bimolecular event
                        TS_fid = self.store_bim_event(comparer)
                        if verbose:
                            print ("###########  Event at %d stored in DB #############" % TS_fid)
                    elif flag == 2:
                        # TBI unimolecular reaction
                        print ("UNIMOL REACTION EVENT!!")
                        raise
                    # we need to make sure that the educts of the reaction event (TS_fid-1)
                    #               are really idential to where we started the forward interval (last reaction event)
                    comparer2 = self.get_comparer(last_event, TS_fid-1)
                    flag = comparer2.check()
                    if flag>0:
                        # there is still an event missing so we need to backtrack again and keep the current event in the open_events
                        open_events.append(TS_fid)   # open events are stored in inverse chronological order ... the latest is first
                        nextf = back_curr
                        mode = "back"
                    else:
                        # These frames are equal which means the last product frame and the current educt frame are identical
                        # => we can continue and connect the tracked species by and edge in the reaction graph
                        mode = "forward"
                        self.connect_events(comparer2)
                        last_event = TS_fid+1
                        while len(open_events) > 0:
                            TS_fid = open_events.pop()                             # take the next evnt to be connected from the registered events
                            comparer2 = self.get_comparer(last_event, TS_fid-1)    # and compare to previous in the chain
                            flag = comparer2.check()
                            assert flag == 0, "THis should never happen ... connecting chain of events and found an unexpected deviation between frames %d and %d" % (last_event, TS_fid-1) 
                            self.connect_events(comparer2)
                            last_event = TS_fid+1
                        # now we are done with the connecting of events -> last_event is set and we can continue
                        nextf = farthestf
                        if verbose:
                            print ("#########################################################")
                            print ("all events stored and connected .. moving forward again"   )
            currentf = nextf
            nextf = currentf+self.sstep[mode]
            if verbose:
                print ("&& current %5d  next %5d  mode %s" % (currentf, nextf, mode))
        # end of mainloop
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
    def store_bim_event(self, comparer):
        """store species and further info in database

        DEBUG: currently only printing results and writing mol objects
        
        Args:
            comparer (fcompare object): current comparer between to frames (f1 is educts)
        """
        for r in range(comparer.nreacs):
            educts, products = comparer.reacs[r]
            broken_bonds     = comparer.broken_bonds[r]
            formed_bonds     = comparer.formed_bonds[r]
            f1               = comparer.f1
            f2               = comparer.f2
            # choose which frame we use as TS
            # everything is referenced to f1 which is 0 (f2 is +1)
            # find avereage bond order of reactive bonds at f1/f2
            f1_averborder = 0.0
            if len(broken_bonds) >0:
                for b in broken_bonds:
                    f1_averborder += f1.molg.ep.bord[f1.molg.edge(b[0], b[1])]
                f1_averborder/len(broken_bonds)
            f2_averborder = 0.0
            if len(formed_bonds) >0:
                for b in formed_bonds:
                    f2_averborder += f2.molg.ep.bord[f2.molg.edge(b[0], b[1])]
                f2_averborder/len(formed_bonds)
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
                TS = f1
                PR = f2
                ED = self.process_frame(TS.fid-1)
                # get corresponding species numbers
                TS_spec = educts
                PR_spec = products
                loccomp = fcompare(ED, TS)
                if loccomp.check_aids() == 0:
                    # no change in atom ids .. we can use TS species for ED as well
                    ED_spec = educts
                else:
                    print ("Houston we have a problem!!! species changed between ED and TS")
            else:
                ED = f1
                TS = f2
                PR = self.process_frame(TS.fid+1)
                ED_spec = educts
                TS_spec = products
                loccomp = fcompare(TS, PR)
                if loccomp.check_aids() == 0:
                    # no change in atom ids .. we can use TS species for PR as well
                    PR_spec = products
                else:
                    print ("Houston we have a problem!!! species changed between TS and PR")
            # now we should have all data
            # generate mol objectes with coordinates
            xyz_ed = np.array(self.f_xyz[ED.fid])
            xyz_ts = np.array(self.f_xyz[TS.fid])
            xyz_pr = np.array(self.f_xyz[PR.fid])
            ED_mol = []
            ED_spec_tracked = []
            # ED/PR_spec are dicitionaries of species .. we need a sorted list of integers
            ED_spec_id = list(ED_spec.keys())
            ED_spec_id.sort()
            for s in ED_spec_id:
                aids = list(ED_spec[s].aids)
                m = ED_spec[s].make_mol(xyz_ed[aids], self.mol)
                ED_mol.append(m)
                if ED_spec[s].tracked:
                    ED_spec_tracked.append(s)
            PR_mol = []
            PR_spec_tracked = []
            PR_spec_id = list(PR_spec.keys())
            PR_spec_id.sort()
            for s in PR_spec_id:
                aids = list(PR_spec[s].aids)
                m = PR_spec[s].make_mol(xyz_pr[aids], self.mol)
                PR_mol.append(m)
                if PR_spec[s].tracked:
                    PR_spec_tracked.append(s)
            # only one TS (all species) .. use make_mol of frame
            TS_spec_id = list(TS_spec.keys())
            TS_spec_id.sort()
            TS_mol, TS_aids = TS.make_mol(TS_spec, xyz_ts, self.mol)
            # map the broken/formed bonds to atom ids of the TS subsystem
            rbonds_global = formed_bonds+broken_bonds
            rbonds = []
            for b in rbonds_global:
                rbonds.append(TS_aids.index(b[0]))
                rbonds.append(TS_aids.index(b[1]))
            # now let us register the reaction event in the database
            revID = self.rdb.register_revent(
                TS.fid,
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
            return TS.fid  # return the TS frame ID determined in this process

    def connect_events(self, comparer):
        """connecting a reaction event 
        
        Args:
            comparer (fcompare object): comparer between last_event (products) and TS_fid-1 (educts of next event)
        """
        for match in comparer.bond_match:
            self.rdb.set_react(comparer.f1.fid, comparer.f2.fid, match[0], match[1])
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










