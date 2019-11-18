"""

           findR

           find reactions in a ReaxFF trajectory

           this code is inspired by chemtrayzer but focuses on heterogeneous gas phase reactions



"""

import numpy as np

import molsys
from molsys.util import pdlpio2

from molsys.util.timing import timer, Timer
from molsys import mpiobject

from graph_tool import Graph, GraphView
import graph_tool.topology as gtt
import graph_tool.util as gtu

class frame:
    """container to store info on a frame

    the main attributes are:
    - fid: frame id 
    - molg:    the molecular graph (edges below threshold are filtered but the bond order is in the edge property bord)
    - species: dictionary with lists of atom indices, the key of the dict is the molecule id (mid vertex property of molg)
    - speciesg: view of molg for a specific species 
    """

    def __init__(self, fid, molg, species, speciesg):
        self.fid = fid
        self.molg = molg
        self.species = species
        self.speciesg = speciesg
        return

    @property
    def nspecies(self):
        return len(self.species.keys())



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
        self.min_atom = 4          # below this number of atoms a species is considered as gasphase
        return

    @timer("process_frame")
    def process_frame(self, i):
        """process a single frame

        Args:
            i (int): number of frame to process
            
        Returns:
            frame: frame object
        """
        bondord = np.array(self.f_bondord[i])
        bondtab = np.array(self.f_bondtab[i])
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
            if e[0] >= 0:
                newb = molg.add_edge(e[0], e[1])
                isbond = o > self.bondord_cutoff
                molg.ep.bord[newb] = o
                molg.ep.filt[newb] = isbond 
            else:
                break
        # apply edge filter
        molg.set_edge_filter(molg.ep.filt)
        mid, hist = gtt.label_components(molg, vprop=molg.vp.mid)
        nspecies_all = len(hist)
        nspecies = []
        for i in range(nspecies_all):
            if hist[i] >= self.min_atom:
                nspecies.append(i)
        species  = {}
        speciesg = {}
        for s in nspecies:
            vs = gtu.find_vertex(molg, molg.vp.mid, s)
            species[s] = [int(v) for v in vs]        
            molg.vp.filt.a[:] = False
            for v in vs:
                molg.vp.filt[v] = True
            speciesg[s] = GraphView(molg, vfilt=molg.vp.filt)     
        f = frame(i, molg, species, speciesg)
        return f
    
    def do_frames(self, start=0, stop=None, stride=1):
        """generate all frames (molg, species) representations
        
        Args:
            start (integer, optional): starting frame. Defaults to None.
            stop (integer, optional): final frame (not inlcuded). Defaults to None.
            stride (integer, optional): stride. Defaults to None.
        """ 
        if not stop:
            stop = self.nframes
        for i in range(start, stop, stride):
            f = self.process_frame(i)
            self.frames.append(f)
        return

    def report(self):
        self.timer.write()










