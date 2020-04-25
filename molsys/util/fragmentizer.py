# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 18:29:19 2016

@author: rochus

          Fragmentizer class

          depends on graph addon (this means graph_tool must be installed)
          
          MPI-safe version .. Warning: use prnt function which is overloaded


        RS (May 2020, corona times)
        This is a complete revision of the fragmentizer. we allow only fragments with more than 1 atom.
        All single atom fragments like halides Me groups, ether Os or thio SHs are assigned "by hand" in a hardcoded way
        Note that some names of fragments are thus defined here and NOT in the MOFplus database.

        - Option for local fragments is removed!!
        - fragments are cached if you do multiple fragmentizations

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

############## single atom fragment rules #################################
# sarules_1 first attempt rules with full atomtype
# sarules_2 second attempt for atoms with only a root atomtype
sarules_1 = {
    "f1_c1"   : "f",         # we define the halogens explicit ... if they are not at C but at what else? N-Cl ?? 
    "cl1_c1"  : "cl",        #                                     --> should be a bigger fragment then
    "br1_c1"  : "br",        #
    "i1_c1"   : "i",         #
    "o2_c2"   : "eth",        # this is an ether and NOT an OH
    "o2_c1h1" : "oh",         # alcohol .. NOT a carboxylic acid .. found as co2-one anyway
    "s2_c1h1" : "sh",         # thiol
    "s2_c2"   : "thio",       # thioether
    "n3_c3"   : "tamin",      # tertiary amine
    "n3_c1h2" : "nh2",        # primary amine 
}

sarules_2 = {
    "c4"      : "me",         # this is a methyl group .. sp3 carbon 
}

class fragmentizer:

    def __init__(self):
        """
        fragmentizer gets a catalog of fragments
        it will use the API to download from MOF+
        """
        # default
        self.fragments = {}          # this is a cache and can be reused if multiple calls are done
        self.frag_vtypes = {}
        self.frag_prio = {}
        # API calls are done on the master only
        if mpi_rank == 0:
            from mofplus import FF_api
            self.api = FF_api()
        else:
            self.api = None
        self.catalog_from_API()
        for f in self.fragments.keys():
            print (f, self.frag_vtypes[f])
        return


    def catalog_from_API(self):
        """
        API call on master only ... broadcasted to other nodes

        Ignore prio == 1 fragments (just to be sure)
        """
        if mpi_rank == 0:
            frags = self.api.list_FFfrags()
        else:
            frags = None
        if mpi_size > 1:
            frags = mpi_comm.bcast(frags, root=0)
        for f in frags:
            # ignore if prio is one
            if f[1] > 1:
                self.fragments[f[0]]= None
                self.frag_vtypes[f[0]] = [vt for vt in f[2] if vt != "h1"]
                self.frag_prio[f[0]] = f[1]
        return

    def read_frags_from_API(self,fnames):
        """
        API call on master only ... broadcsted to other nodes 

        New (2020):
        gets only those that are not in the cache already
        """
        fnames_get = [f for f in fnames if self.fragments[f] is None]
        if len(fnames_get) > 0:
            print ("Getting fragments from MOF+:")
            print (fnames_get)
            if mpi_rank == 0:
                mstr = self.api.get_FFfrags(fnames_get, out = "str")
            else:
                mstr = {}
            if mpi_size > 1:
                mstr = mpi_comm.bcast(mstr, root=0)
            for name,lines in mstr.items():
                m = molsys.mol.from_string(lines)
                m.addon("graph")
                m.graph.make_graph(hashes=False)
                self.fragments[name] = m
        return

    def __call__(self, mol, plot=False, get_missing=False, verbose=True):
        """
        tries to assign all fragmnets in the catalog to the mol object

        :Parameters:

            - mol  (mol object): mol object to be fragmentized
            - plot (bool, opt): write fragment graphs as png file (defualts to False)
            - get_missing (bool, opt): if True analyze in case of failure and propose missing fragments 
            - verbose (bool, opt): talk to the user (defaults to True) 

        """
        if verbose:
            print ("This is Fragmentizer")
            print ("====================")
        # set all fragment info to none
        mol.set_nofrags()
        #
        mol.addon("graph")
        mol.graph.make_graph(hashes=False)
        if plot:
            mol.graph.plot_graph(plot, ptype="png", vsize=20, fsize=20)
        # get list of atypes
        atypes = mol.get_atypelist()
        vtype = map(lambda e: e.split("_")[0], atypes)
        vtype = filter(lambda e: (e[0] != "x") and (e[0] != "h"), vtype)
        vtype = list(set(vtype))
        if verbose:
            print ("The system contains the following root atomtypes for which we test fragments:")
            print (vtype)
        # scan for relevant fragments
        scan_frag = []
        scan_prio = []
        tmpnames = []
        for fname in self.fragments.keys():
            # check if all vtypes in frag appear in the systems vtype
            if all(v in vtype for v in self.frag_vtypes[fname]):
                scan_frag.append(fname)
                scan_prio.append(self.frag_prio[fname])
                if self.fragments[fname] is None:
                    tmpnames.append(fname)
        self.read_frags_from_API(tmpnames)
        # now sort according to prio
        sorted_scan_frag = [scan_frag[i] for i in numpy.argsort(scan_prio)]
        sorted_scan_frag.reverse()
        # now run over the system and test the fragments
        atypes = mol.get_atypes()
        fragnames = []        # this a list of all existing fragment names
        frag_atoms = []       # a list of the fragments with their atoms
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
                    fragnames.append(fname)
                    frag_atoms.append(alist)
                    fi += 1
        # now all retrieved frags are tested .. now try to assign the remaining atoms with fragments if possible
        for i in range(mol.natoms):
            if mol.fragtypes[i] == '-1':
                at  = atypes[i] # atomtype
                rat = at.split("_")[0] # root atomtype
                # this is unassigned .. is it hydrogen?
                if rat != "h1":
                    fname = None
                    # first test explicit sarules_1
                    if at in sarules_1:
                        fname = sarules_1[at]
                    else:
                        if rat in sarules_2:
                            fname = sarules_2[rat]
                    # did we find something that applies?
                    if fname is not None:
                        # all hydogens connected to this atom are also part of the fragment (check if already assigned .. should not be)
                        alist = [i]
                        for j in mol.conn[i]:
                            if atypes[j].split("_")[0] == "h1":
                                alist.append(j)
                                assert mol.fragtypes[j] == "-1"
                        for i in alist:
                            mol.fragtypes[i] = fname
                            mol.fragnumbers[i] = fi
                        fragnames.append(fname)
                        frag_atoms.append(alist)
                        fi += 1
        nfrags =  max(mol.fragnumbers)+1
        # now analyse if fragmentation was successful and propse missing fragments if get_missing is True
        if verbose:
            print ("The following fragments have been found in the system")
            print (set(fragnames))
            nuassigned = mol.fragtypes.count('-1')
            if nuassigned > 0:
                print ("!!!! Unassigned atoms: %d    mol %s" % (nuassigned, mol.name))
        if get_missing:
            # try to identify groups of atoms
            pass

        return

    # this method should go evetually and is left here only for legacy reasons (who needs it as a static method?)
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
                if fraglist[fn] is None:
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

    