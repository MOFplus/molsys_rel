"""

           findR

           find reactions in a ReaxFF trajectory

           this code is inspired by chemtrayzer but focuses on heterogeneous gas phase reactions
           2019 RUB Rochus Schmid



"""

import numpy as np
import pickle
import copy

import molsys
from molsys.util import pdlpio2

from molsys.util.timing import timer, Timer
from molsys.util.print_progress import print_progress
from molsys import mpiobject

from graph_tool import Graph, GraphView
import graph_tool.topology as gtt
import graph_tool.util as gtu
import graph_tool.draw as gtd

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
        return species(self.fid, mid, self.molg, make_graph=True) # for DEBUG .. remove make_graph

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


class species:
    """container class to keep species info (per frame!)
    """

    def __init__(self, fid, mid, molg, make_graph=False):
        """init species
        
        Args:
            fid (int): frame number
            mid (int): molecule id ("name" of species from label components)
            molg (graph): molg of the frame
        """
        self.fid = fid
        self.mid = mid
        # find all vertices in molg that belong to this species mid
        vs = gtu.find_vertex(molg, molg.vp.mid, mid)
        self.aids = set([int(v) for v in vs]) # atomids -> TBI do we need to sort? seems like they are sroted as they come
        self.graph = None
        if make_graph:
            # now make a view of molg for this species         
            molg.vp.filt.a[:] = False
            for v in vs:
                molg.vp.filt[v] = True
            self.graph = GraphView(molg, vfilt=molg.vp.filt)
        return

    @property
    def natoms(self):
        return len(self.aids)

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
        self.reacs = []
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
                        # we need to make a new species object and add it
                        educts[esk] = self.f1.make_species(esk)
                        educt_aids |= educts[esk].aids
                # which atoms are in the educts that are not in the product species? add them
                for a in educt_aids - product_aids:
                    # to which species in frame 2 does this atom belong to?
                    psk = self.f2.molg.vp.mid[a]
                    # is this already in educts?
                    if psk not in products:
                        # we need to make a new species object and add it
                        products[psk] = self.f2.make_species(psk)
                        product_aids |= products[psk].aids
            # now add the final results to the reacs list
            self.reacs.append((educts, products))
        return

    def find_react_bond(self, r):
        """find a reactive bond in a bimolecular reaction 
        
        Args:
            r (int): index in self.reac to analyse
        """
        educts, products = self.reacs[r]
        # get all involved atoms
        aids = set()
        for s in educts:
            aids |= educts[s].aids
        g1 = self.f1.molg
        g2 = self.f2.molg
        reactive_bond = []
        for a in aids:
            bonds1 = []
            v1 = g1.vertex(a)
            for e in v1.out_edges(): 
                bonds1.append(int(e.target()))
            bonds2 = []
            v2 = g2.vertex(a)
            for e in v2.out_edges():
                bonds2.append(int(e.target()))
            s = set(bonds1) ^ set(bonds2)
            for b in s:
                if b > a:
                    reactive_bond.append((a,b))
        return reactive_bond

               

class findR(mpiobject):

    def __init__(self, pdlpfilename,  stage, mpi_comm = None, out = None):
        super(findR,self).__init__(mpi_comm, out)
        # To be tested: open read only on all nodes in parallel .. does it work?
        self.pdlp = pdlpio2.pdlpio2(pdlpfilename, filemode="r")
        self.mol = self.pdlp.get_mol_from_system()
        assert stage in self.pdlp.get_stages(), "Stage %s not in stages in file"
        self.traj = self.pdlp.h5file["/%s/traj" % stage]
        data = list(self.traj.keys())
        # make sure that xyz and bondord/bondterm is available
        assert "xyz" in data
        assert "bondord" in data
        assert "bondtab" in data
        self.f_xyz = self.traj["xyz"]
        self.f_bondord = self.traj["bondord"]
        self.f_bondtab = self.traj["bondtab"]
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
        self.smode = "forward" # must be either forward, back or fine
        self.scurrent = 0
        self.snext    = 0
        self.sstep ={
            "forward" : 200,
            "back"    : -10,
            "fine"    : 1,
        }
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
    def search_frames(self):
        """search the frames for reactions 

        This is the core routine of findR
        """
        self.smode = "forward"
        self.scurrent = 0
        self.process_frame(self.scurrent)
        self.snext = self.scurrent+self.sstep[self.smode]
        # start mainloop (just search on until we are at the end of the file)
        while self.snext < self.nframes and self.snext > 0:
            self.process_frame(self.snext)
            # do comparison and report
            comparer = self.get_comparer(self.scurrent, self.snext)
            flag = comparer.check_aids()
            print ("&&FLAG %d" % flag)
            #DEBUG DEBUG
            # comparer.report(all=True)
            #DEBUG DEBUG
            if flag>0:
                # a difference was found between current and next
                if self.smode == "forward":
                    forw_next = self.snext
                    forw_curr = self.scurrent
                    self.smode = "back"
                elif self.smode == "back":
                    back_next = self.snext
                    back_curr = self.scurrent
                    self.smode = "fine"
                elif self.smode == "fine":
                    comparer.analyse_aids()
                    print ("##### REACTION EVENT FOUND!!!!")
                    comparer.report()
                    #
                    # TBI -> search critical bond
                    assert len(comparer.reacs)==1
                    print (comparer.find_react_bond(0))
                    # TBI -> go back and forth to locate safe educt and product
                    # TBI -> store in database
                    #
                    # we need to make sure that the educts of the reaction event (current frame)
                    #               are really idential to where we started the forward interval
                    comparer2 = self.get_comparer(forw_curr, self.scurrent)
                    flag = comparer2.check_aids()
                    if flag>0:
                        # there is still an event missing so we need to backtrack again
                        self.snext = back_next
                        self.smode = "back"
                    else:
                        print ("##### Continue searching")
                        self.smode = "forward"
                        self.snext = forw_next
            self.scurrent = self.snext
            self.snext = self.scurrent+self.sstep[self.smode]
            print ("&& current %5d  next %5d  step %5d  mode %s" % (self.scurrent, self.snext, self.sstep[self.smode], self.smode))
        # end of mainloop
        return



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










