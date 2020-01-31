"""

           findR

           find reactions in a ReaxFF trajectory

           this code is inspired by chemtrayzer but focuses on heterogeneous gas phase reactions
           2019 RUB Rochus Schmid



"""

import numpy as np
import pickle

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

    def __init__(self, fid, molg, specs):
        self.fid = fid
        self.molg = molg
        self.specs = specs
        return

    @property
    def nspecies(self):
        return len(self.specs.keys())

    def plot(self):
        print ("plotting species of frame %d" % self.fid)
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

    def __init__(self, fid, mid, molg):
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
        self.aids = [int(v) for v in vs] # atomids -> TBI do we need to sort? seems like they are sroted as they come
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
        return


    def check_species(self):
        """this method is just to implement all levels of comparison

           it does not return anything and just reports ... meant for debugging
        """
        print("##################################################")
        print("FRAMES        %5d         %5d" % (self.f1.fid, self.f2.fid))
        print("level 0:      %5d         %5d species" % (self.f1.nspecies, self.f2.nspecies))
        # map species ... anything that remains has changed
        sk1l = list(self.f2.specs.keys())
        sk2l = list(self.f2.specs.keys())
        aid_match = [] 
        for sk1 in self.f1.specs:
            s1 = self.f1.specs[sk1]
            # find a corresponding species in sk2
            for sk2 in sk2l:
                s2 = self.f2.specs[sk2]
                if s1.aids == s2.aids:
                    # s1 and s2 match -> add to match and to remove
                    aid_match.append((sk1, sk2))
                    sk1l.remove(sk1)
                    sk2l.remove(sk2)
                    break
        print("level 1: (atom ID match)")
        for m in aid_match:
            print("          species %5d == species %5d" % m)
        for s in sk1l:
            print("          species %5d in frame %d unmatched" % (s, self.f1.fid))
        for s in sk2l:
            print("          species %5d in frame %d unmatched" % (s, self.f2.fid))
        return


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
            "back"    : 10,
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
        specs = {}
        for mid in nspecies:
            specs[mid] = species(fid, mid, molg)
        # TBI shall we reset the vertex filter of molg again?      
        f = frame(fid, molg, specs)
        # store the frame
        self.frames[fid] = f
        return f
    
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
        while self.snext < self.nframes:
            self.process_frame(self.snext)
            # do comparison and report
            comparer = fcompare(self.frames[self.scurrent], self.frames[self.snext])
            comparer.check_species()
            self.scurrent = self.snext
            self.snext = self.scurrent+self.sstep[self.smode]
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










