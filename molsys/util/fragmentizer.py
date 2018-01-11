# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 18:29:19 2016

@author: rochus

          Fragmentizer class

          depends on graph addon (this means graph_tool must be installed)
          
          MPI-safe version .. Warning: use prnt function which is overloaded

"""
from __future__ import print_function
import os
import numpy
import logging
import glob
import molsys
import csv
from . import atomtyper

try:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    mpi_rank = MPI.COMM_WORLD.Get_rank()
    mpi_size = MPI.COMM_WORLD.Get_size()
except ImportError as e:
    mpi_comm = None
    mpi_size = 1
    mpi_rank = 0
    mpi_err = e

# overload print function in parallel case
try:
    import __builtin__
except ImportError:
    import builtins as __builtin__
def print(*args, **kwargs):
    if mpi_rank == 0:
        return __builtin__.print(*args, **kwargs)
    else:
        return


import logging

logger = logging.getLogger("molsys.fragmentizer")

if mpi_comm is None:
    logger.error("MPI NOT IMPORTED DUE TO ImportError")
    logger.error(mpi_err)

class fragmentizer:

    def __init__(self, source="mofp"):
        """
        fragmentizer gets a catalog of fragments
        if source is "file" it reads from local disk either
              from the current directory or from $MOLSYS_FRAGS
              a catalog is expected to be in fragments.csv
        if source is "mofp" it will use the API to download from MOF+

        :Paramters:

            - source: either "file" or "mofp" [default "mofp"]
        """
        # default
        self.fragments = {}
        self.frag_vtypes = {}
        self.frag_prio = {}
        self.source = source
        if source == "file":
            if os.environ.has_key("MOLSYS_FRAGS"):
                self.frag_path = os.environ["MOLSYS_FRAGS"]
            else:
                self.frag_path = "."
            self.read_catalog()
        elif source == "mofp":
            # API calls are done on the master only
            if mpi_rank == 0:
                from mofplus import FF_api
                self.api = FF_api()
            else:
                self.api = None
            self.catalog_from_API()
        else:
            raise ValueError("Unknown source specified")
        return

    def read_catalog(self):
        """
        file-mode: read available fragments from csv file
        """
        f = open(self.frag_path + "/fragments.csv", "rb")
        csvf = csv.reader(f, delimiter=",")
        for row in csvf:
            fname  = row[0]
            vtypes = row[1].split()
            prio   = int(row[2])
            self.fragments[fname] = None
            self.frag_vtypes[fname] = vtypes
            self.frag_prio[fname] = prio
        f.close()
        return

    def catalog_from_API(self):
        """
        API call on master only ... broadcasted to other nodes
        """
        if mpi_rank == 0:
            frags = self.api.list_FFfrags()
        else:
            frags = None
        if mpi_size > 1:
            frags = mpi_comm.bcast(frags, root=0)
        for f in frags:
            self.fragments[f[0]]= None
            self.frag_vtypes[f[0]] = f[2]
            self.frag_prio[f[0]] = f[1]
        return

    def read_frag(self, fname):
        """
        file-mode: read a fragment and convert to a graph
        """
        m = molsys.mol()
        m.read(self.frag_path + "/" + fname + ".mfpx", ftype="mfpx")
        m.addon("graph")
        m.graph.make_graph()
        self.fragments[fname] = m
        return

    def read_frag_from_API(self,fname):
        """
        API call on master only ... broadcsted to other nodes 
        """
        if mpi_rank == 0:
            m = self.api.get_FFfrag(fname, mol = True)
        else:
            m = []
        if mpi_size > 1:
            m = mpi_comm.bcast(m, root=0)
        m.addon("graph")
        m.graph.make_graph()
        self.fragments[fname] = m
        return

    def __call__(self, mol, man = False, plot=False):
        """
        tries to assign all fragmnets in the catalog to the mol object

        :Parameters:

            - mol : mol object to be fragmentized

        """
        # set all fragment info to none
        mol.set_nofrags()
        #
        mol.addon("graph")
        mol.graph.make_graph()
        if plot:
            mol.graph.plot_graph(plot, ptype="png", vsize=20, fsize=20)
        mol.set_nofrags()
        # get list of atypes
        atypes = mol.get_atypelist()
        vtype = map(lambda e: e.split("_")[0], atypes)
        vtype = filter(lambda e: (e[0] != "x") and (e[0] != "h"), vtype)
        vtype = list(set(vtype))
        # print(vtype)
        # scan for relevant fragments
        scan_frag = []
        scan_prio = []
        for fname in self.fragments.keys():
            # check if all vtypes in frag appear in the systems vtype
            if all(v in vtype for v in self.frag_vtypes[fname]):
                scan_frag.append(fname)
                scan_prio.append(self.frag_prio[fname])
                if self.fragments[fname] == None:
                    # not read in yet
                    if self.source == "file":
                        self.read_frag(fname)
                    elif self.source == "mofp":
                        self.read_frag_from_API(fname)
                    else:
                        raise ValueError("unknwon source for fragments")
        # now sort according to prio
        sorted_scan_frag = [scan_frag[i] for i in numpy.argsort(scan_prio)]
        sorted_scan_frag.reverse()
        # now run over the system and test the fragments
        atypes = mol.get_atypes()
        fi = 0
        for fname in sorted_scan_frag:
            fidx = mol.graph.find_fragment(self.fragments[fname],add_hydrogen=True)
            for alist in fidx:
                # if any of the atoms in alist is already in a fragment we can skip
                assigned_already = any(mol.fragnumbers[i] >= 0 for i in alist)
                if not assigned_already:
                    if plot:
                        self.fragments[fname].graph.plot_graph(fname, ptype="png", vsize=20, fsize=20, size=400)
                    for i in alist:
                        mol.fragtypes[i]   = fname
                        mol.fragnumbers[i] = fi
                        # print "atom %s set to fragment %s" % (atypes[i], fname)
                    fi += 1
        ### patch for working with pyr/dabco ligands concerning 
        if self.pure_check(mol):
            logger.info("Fragmentation was successful")
            if man == True:
                atyper = atomtyper.atomtyper(mol)
                atyper()
        else:
            if man == True:
                logger.error("Fragmentation failed!!!")
            else:
                logger.error("Fragmentation failed --> manipulate atypes")
                ### manipulate atypes
                manipulated = False
                for i,at in enumerate(atypes):
                    if at == "n2_c2" and mol.fragtypes[i] == "-1": 
                        manipulated = True
                        mol.atypes[i] = "n3_m"
                    elif at == "n3_c3" and mol.fragtypes[i] == "-1": 
                        manipulated = True
                        mol.atypes[i] = "n4_m"
                if manipulated: self.__call__(mol, man = manipulated)
        return

    @staticmethod
    def pure_check(mol):
        """
        check if all atoms are in a fragment and all is consistent
        """
        fragnames = []        # this a list of all existing fragment names
        frag_atoms = []       # a list of the fragments with their atoms
        nfrags =  max(mol.fragnumbers)+1
        fraglist  = [None]*(nfrags) # additional list to hold the fragments with their name
        for i in range(nfrags):
            frag_atoms.append([])
        for i in range(mol.natoms):
            ft = mol.fragtypes[i]
            fn = mol.fragnumbers[i]
            if ft == "0":
                return False
            else:
                if fraglist[fn] == None:
                    # set name of fragment
                    fraglist[fn] = ft
                else:
                    # check if this is the same name
                    if fraglist[fn] != ft: return False
                if ft not in fragnames:
                    fragnames.append(ft)
                frag_atoms[mol.fragnumbers[i]].append(i)
        # in the end make sure that all fragments have been named
        if None in fraglist: return False
        return True


