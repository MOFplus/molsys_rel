# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 11:25:43 2017

@author: rochus


        addon module FF to implement force field infrastructure to the molsys

"""


import numpy as np
from molsys.util.timing import timer, Timer

import itertools

import logging
logger = logging.getLogger("molsys.ff")


class ric:
    """
    class to detact and keep all the (redundant) internal coordinates of the system
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

    def assign_params(self, FF, source="mofp"):
        """
        method to orchestrate the parameter assigenemntfor this sytem using a force ifled defined with
        FF

        :Parameter:

            - FF :    [string] name of the force field to be used in the parameter search
            - source: [string, default=mofp] where to get the data from
        """
        self.source = source
        if self.source == "mofp":
            from mofplus import FF_api
            self.api = FF_api()
        elif self.source == "file":
            raise ValueError, "to be implemented"
        else:
            raise ValueError, "unknown source for assigning parameters"
        # as a first step we need to generate the fragment graph
        self.timer.start("fragment graph")
        self._mol.addon("fragments")
        self.fragments = self._mol.fragments
        self.fragments.make_frag_graph()
        fragnames  = self.fragments.get_fragnames()
        self.timer.stop()
        # now we load all the reference systems of relevance
        if source == "mofp":
            self.timer.start("get reference systems")
            scan_ref  = []
            scan_prio = []
            ref_list = self.api.list_FFrefs(FF)
            for ref in ref_list:
                refname, prio, reffrags = ref
                if len(reffrags) > 0 and all(f in fragnames for f in reffrags):
                    scan_ref.append(refname)
                    scan_prio.append(prio)
            # sort to be scanned referecnce systems by their prio
            scan_ref = [scan_ref[i] for i in np.argsort(scan_prio)]
            scan_ref.reverse()
            self.timer.stop()
            # now get the refsystems and make their fraggraphs
            self.timer.start("get ref frag graphs")
            self.ref_systems = {}
            for ref in scan_ref:
                ref_mol = self.api.get_FFref_graph(ref, mol=True)
                ref_mol.addon("fragments")
                ref_mol.fragments.make_frag_graph()
                self.ref_systems[ref] = ref_mol
            self.timer.stop()
            # now search in the fraggrpah for the reference systems
            self.timer.start("scan for ref systems")
            self.ref_fraglists = {}
            for ref in self.ref_systems.keys():
                # TODO: if a ref system has only one fragment we do not need to do a substructure search but
                #       could pick it from self.fragemnts.fraglist
                subs = self._mol.graph.find_subgraph(self.fragments.frag_graph, self.ref_systems[ref].fragments.frag_graph)
                if len(subs) == 0:
                    # this ref system does not appear => discard
                    scan_ref.remove(ref)
                    del(self.ref_systems[ref])
                # join all fragments
                subs_flat = itertools.chain.from_iterable(subs)
                self.ref_fraglists[ref] = list(set(subs_flat))
            self.timer.stop()

            for r in scan_ref:
                print r
                print self.ref_fraglists[r]

            self.timer.write_logger(logger.info)
            return

