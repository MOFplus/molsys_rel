# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 11:25:43 2017

@author: rochus


        addon module FF to implement force field infrastructure to the molsys

"""


import numpy as np
from molsys.util.timing import timer, Timer
from mofplus import aftype, aftype_sort

import itertools
import copy
import string
import json

import logging
logger = logging.getLogger("molsys.ff")


class ric:
    """
    class to detect and keep all the (redundant) internal coordinates of the system
    """


    def __init__(self, mol):
        """
        find all bonds, angles, oops and dihedrals in the molsys mol object
        """

        self.timer = Timer()

        self.conn   = mol.get_conn()
        self.natoms = mol.get_natoms()

        self.bnd    = self.find_bonds()
        self.ang    = self.find_angles()
        self.oop    = self.find_oops()
        self.dih    = self.find_dihedrals()
        self.report()
        self.timer.write_logger(logger.info)
        return


    # the follwoing code is adapted from pydlpoly/py/assign_FF.py
    # instead of bond objects we use simpler lists of lists with indices

    @timer("find bonds")
    def find_bonds(self):
        bonds=[]
        for a1 in xrange(self.natoms):
            for a2 in self.conn[a1]:
                if a2 > a1: bonds.append([a1, a2])
        return bonds

    @timer("find angles")
    def find_angles(self):
        angles=[]
        for ca in xrange(self.natoms):
            apex_atoms = self.conn[ca]
            naa = len(apex_atoms)
            for ia in xrange(naa):
                aa1 = apex_atoms[ia]
                other_apex_atoms = apex_atoms[ia+1:]
                for aa2 in other_apex_atoms:
                    angles.append([aa1, ca, aa2])
        return angles

    @timer("find oops")
    def find_oops(self):
        oops=[]
        # there are a lot of ways to find oops ...
        # we assume that only atoms with 3 partners can be an oop center
        for ta in xrange(self.natoms):
            if (len(self.conn[ta]) == 3):
                # ah! we have an oop
                a1, a2, a3 = tuple(self.conn[ta])
                oops.append([ta, a1, a2, a3])
                #oops.append([ta, a2, a1, a3])
                #oops.append([ta, a3, a1, a2])
        return oops

    @timer("find dihedrals")
    def find_dihedrals(self):
        dihedrals=[]
        for a2 in xrange(self.natoms):
            for a3 in self.conn[a2]:
                # avoid counting central bonds twice
                if a3 > a2:
                    endatom1 = list(self.conn[a2])
                    endatom4 = list(self.conn[a3])
                    endatom1.remove(a3)
                    endatom4.remove(a2)
                    for a1 in endatom1:
                        con1 = list(self.conn[a1])
                        for a4 in endatom4:
                            # code to detect 4 and 5 rings
                            #smallring = 0
                            #if con1.count(a4):
                            #    smallring = 4
                            #else:
                            #    con4 = list(self.conn[a4])
                            #    for c1 in con1:
                            #        if con4.count(c1):
                            #            smallring = 5
                            #            break
                            if a1 != a4: dihedrals.append([a1,a2,a3,a4])
        return dihedrals

    def report(self):
        logger.info("Reporting RICs")
        logger.info("%7d bonds"     % len(self.bnd))
        logger.info("%7d angles"    % len(self.ang))
        logger.info("%7d oops"      % len(self.oop))
        logger.info("%7d dihedrals" % len(self.dih))
        return



class ff:

    def __init__(self, mol):
        """
        instantiate a ff object which will be attached to the parent mol

        :Parameter:

             - mol : a mol type object (can be a derived type like bb or topo as well)
        """

        self._mol = mol
        self.ric = ric(mol)
        self.timer = Timer()
        logger.debug("generated the ff addon")
        return

    @timer("assign parameter")
    def assign_params(self, FF, source="mofp"):
        """
        method to orchestrate the parameter assignment for this system using a force field defined with
        FF

        :Parameter:

            - FF :    [string] name of the force field to be used in the parameter search
            - source: [string, default=mofp] where to get the data from
        """
        self.timer.start("connect to DB")
        self.FF = FF
        self.source = source
        if self.source == "mofp":
            from mofplus import FF_api
            self.api = FF_api()
        elif self.source == "file":
            raise ValueError, "to be implemented"
        else:
            raise ValueError, "unknown source for assigning parameters"
        self.timer.stop()
        # as a first step we need to generate the fragment graph
        self.timer.start("fragment graph")
        self._mol.addon("fragments")
        self.fragments = self._mol.fragments
        self.fragments.make_frag_graph()
        self.timer.stop()
        # now make a private list of atom types including the fragment name
        self.timer.start("make atypes")
        self.aftypes = []
        for i, a in enumerate(self._mol.get_atypes()):
            self.aftypes.append(aftype(a, self._mol.fragtypes[i]))
        self.timer.stop()
        # detect refsystems
        self.find_refsystems()
        # make data structures
        with self.timer("make data structures"):
            ric_type = {"cha":[[i] for i in range(self._mol.natoms)],
                    "vdw":[[i] for i in range(self._mol.natoms)], 
                    "bnd":self.ric.bnd, 
                    "ang":self.ric.ang, 
                    "dih":self.ric.dih, 
                    "oop":self.ric.oop}
            self.parind = {
                "bnd": [None]*len(self.ric.bnd),
                "ang": [None]*len(self.ric.ang),
                "dih": [None]*len(self.ric.dih),
                "oop": [None]*len(self.ric.oop),
                "cha": [None]*self._mol.natoms,
                "vdw": [None]*self._mol.natoms,
                }
            self.par = {
                "bnd": {},
                "ang": {},
                "dih": {},
                "oop": {},
                "vdw": {},
                "cha": {},
                }
        with self.timer("parameter assignement loop"):
            for ref in self.scan_ref:
                logger.info("assigning params for ref system %s" % ref)
                curr_fraglist = self.ref_fraglists[ref]
                curr_par = {\
                    "bnd" : self.ref_params[ref]["twobody"]["bnd"],\
                    "ang" : self.ref_params[ref]["threebody"]["ang"],\
                    "dih" : self.ref_params[ref]["fourbody"]["dih"],\
                    "oop" : self.ref_params[ref]["fourbody"]["oop"],
                    "cha" : self.ref_params[ref]["onebody"]["charge"],
                    "vdw" : self.ref_params[ref]["onebody"]["vdw"]
                    }
                for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
                    for i, r in enumerate(ric_type[ic]):
                        if self.parind[ic][i] == None:
                            if self.atoms_in_subsys(r, curr_fraglist):
                                # no params yet and in current refsystem => check for params
                                params = []
                                ptypes = []
                                full_parname_list = []
                                # check unsorted list first
                                parname = self.get_parname(r)
                                if parname in curr_par[ic]:
                                    for par in curr_par[ic][parname]:
                                        ptypes.append(par[0])
                                        full_parname = par[0] + "->" + str(parname) + "|" + ref
                                        full_parname_list.append(full_parname)
                                        if not full_parname in self.par[ic]:
                                            self.par[ic][full_parname] = par
                                # now check sorted list (if already in ptype skip), only for manybody ics
                                if ic not in ["cha", "vdw"]:
                                    parname = self.get_parname_sort(r, ic)
                                    if parname in curr_par[ic]:
                                        for par in curr_par[ic][parname]:
                                            if not par[0] in ptypes:
                                                ptypes.append(par[0])
                                                full_parname = par[0] + "->" + str(parname) + "|" + ref
                                                full_parname_list.append(full_parname)
                                                if not full_parname in self.par[ic]:
                                                    self.par[ic][full_parname] = par
                                if full_parname_list != []:
                                    self.parind[ic][i] = full_parname_list
                                #else:
                                #    print "DEBUG DEBUG DEBUG %s" % ic
                                #    print self.get_parname(r)
                                #    print self.get_parname_sort(r, ic)
                #EQUIVALENCE
                # now all params for this ref have been assigned ... any equivalnce will be renamed now in aftypes
                curr_equi_par = self.ref_params[ref]["onebody"]["equil"]
                for i,a in enumerate(copy.copy(self.aftypes)):
                    if self.atoms_in_subsys([i], curr_fraglist):
                        # NOTE: new aftype keys are now always tuples of aftypes
                        if (a,) in curr_equi_par:
                            # revised aftype : curr_equi_par[a] is a list of a tuple with ("equiv", ("new aftype"))
                            at, ft = curr_equi_par[(a,)][0][1][0].split("@")
                            self.aftypes[i] = aftype(at, ft)
        # DEBUG DEBUG
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            print "Located Paramters for %3s" % ic
            for k in self.par[ic].keys(): print k
            unknown_par = []
            for i, p in enumerate(ric_type[ic]):
                if self.parind[ic][i] == None:
                    parname = self.get_parname(p)
                    if not parname in unknown_par:
                        unknown_par.append(parname)
            for p in unknown_par: print "No params for %3s %s" % (ic, p)
        self.setup_pair_potentials()
        self.timer.write_logger(logger.info)
        return
    
    def setup_pair_potentials(self, radfact = 1.0, radrule="arithmetic", epsrule="geometric"):
        """
        Method to setup the pair potentials based on the per atom type assigned parameters
        :Parameters:
            - radfact (int): factor to be multiplied during radius generation, default to 1.0
            - radrule (str): radiusrule, default to arithmetic
            - epsrule (str): epsilonrule, default to geometric
        """
        self.vdwdata = {}
        self.types2numbers = {} #equivalent to self.dlp_types
        types = self.par["vdw"].keys()
        for i, t in enumerate(types):
            if t not in self.types2numbers.keys():
                self.types2numbers[t]=str(i)
        ntypes = len(types)
        for i in xrange(ntypes):
            for j in xrange(i, ntypes):
                pair = types[i]+":"+types[j]
                #TODO check availability of an explicit paramerter
                par_i = self.par["vdw"][types[i]][1]
                par_j = self.par["vdw"][types[j]][1]
                if radrule == "arithmetic":
                    rad = radfact+(par_i[0]+par_j[0])
                elif radrule == "geometric":
                    rad = 2.0*radfact*np.sqrt(par_i[0]*par_j[0])
                else:
                    raise IOError("Unknown radius rule %s specified" % radrule)
                if epsrule == "arithmetic":
                    eps = 0.5 * (par_i[1]+par_j[1])
                elif epsrule == "geometric":
                    eps = np.sqrt(par_i[1]*par_j[1])
                else:
                    raise IOError("Unknown radius rule %s specified" % radrule)
                self.vdwdata[pair] = [rad,eps]
        return


    @timer("find reference systems")
    def find_refsystems(self):
        """
        function to detect the reference systems:
            - self.scan_ref      : list of ref names in the order to be searched
            - self.ref_systems   : dictionary of mol objects
            - self.ref_fraglist  : list of fragment indices belonging to this refsystem
            - self.ref_params    : paramtere dictionaries per refsystem (n-body/type)
        """
        if self.source == "mofp":
            self.timer.start("get reference systems")
            scan_ref  = []
            scan_prio = []
            ref_list = self.api.list_FFrefs(self.FF)
            for ref in ref_list:
                refname, prio, reffrags = ref
                if len(reffrags) > 0 and all(f in self.fragments.get_fragnames() for f in reffrags):
                    scan_ref.append(refname)
                    scan_prio.append(prio)
            # sort to be scanned referecnce systems by their prio
            self.scan_ref = [scan_ref[i] for i in np.argsort(scan_prio)]
            self.scan_ref.reverse()
            self.timer.stop()
            # now get the refsystems and make their fraggraphs
            self.timer.start("make ref frag graphs")
            self.ref_systems = {}
            for ref in self.scan_ref:
                ref_mol = self.api.get_FFref_graph(ref, mol=True)
                ref_mol.addon("fragments")
                ref_mol.fragments.make_frag_graph()
                self.ref_systems[ref] = ref_mol
            self.timer.stop()
            # now search in the fraggrpah for the reference systems
            self.timer.start("scan for ref systems")
            logger.info("Searching for reference systems:")
            self.ref_fraglists = {}
            for ref in copy.copy(self.scan_ref):
                # TODO: if a ref system has only one fragment we do not need to do a substructure search but
                #       could pick it from self.fragemnts.fraglist
                subs = self._mol.graph.find_subgraph(self.fragments.frag_graph, self.ref_systems[ref].fragments.frag_graph)
                logger.info("   -> found %5d occurences of reference system %s" % (len(subs), ref))
                if len(subs) == 0:
                    # this ref system does not appear => discard
                    self.scan_ref.remove(ref)
                    del(self.ref_systems[ref])
                # join all fragments
                subs_flat = itertools.chain.from_iterable(subs)
                self.ref_fraglists[ref] = list(set(subs_flat))
            self.timer.stop()
            # get the parameters
            self.timer.start("get ref parmeter sets")
            self.ref_params = {}
            for ref in self.scan_ref:
                logger.info("Getting params for %s" % ref)
                self.ref_params[ref] = self.api.get_params_from_ref(self.FF, ref)
                print "DEBUG DEBUG Ref system %s" % ref
                print self.ref_params[ref]
            self.timer.stop()
        elif source == "file":
            raise ValueError, "assigning reference systems from file needs to be implemeted"
        else:
            raise ValueError, "unknown source %s" % source
        return


    def atoms_in_subsys(self, alist, fsubsys):
        """
        this helper function checks if all fragments of atoms (indices) in alist
        appear in the list of fragemnts (indices) in fsubsys
        """
        return all(f in fsubsys for f in map(lambda a: self._mol.fragnumbers[a], alist))

    def get_parname(self, alist):
        """
        helper function to produce the name string using the self.aftypes
        """
        l = map(lambda a: self.aftypes[a], alist)
        return tuple(l)

    def get_parname_sort(self, alist, ic):
        """
        helper function to produce the name string using the self.aftypes
        """
        l = map(lambda a: self.aftypes[a], alist)
        return tuple(aftype_sort(l, ic))
