# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 21:30:26 2016

@author: rochus
"""

import logging
logger = logging.getLogger("molsys.fragments")


class fragments:
    """
    fragments is an addon class to support advanced fragment
    handling
    """


    def __init__(self, mol):
        self._mol = mol
        self.setup = False
        self.check()
        self.make_frag_conn()
        return

    def check(self):
        """
        check if all atoms are in a fragment and all is consistent
        """
        self.fragnames = [] # this a list of all existing fragment names
        self.frag_atoms = [] # a list of the fragments with their atoms
        self.nfrags =  max(self._mol.fragnumbers)+1               
        self.fraglist  = [None]*(self.nfrags) # additional list to hold the fragments with their name
        for i in xrange(self.nfrags):
            self.frag_atoms.append([])
        self.setup = True
        for i in xrange(self._mol.natoms):
            ft = self._mol.fragtypes[i]
            fn = self._mol.fragnumbers[i]
            if ft == "0":
                logger.error("atom %d (%s) is not in fragment" % (i, self._mol.atypes[i]))
                self.setup=False
            else:
                if self.fraglist[fn] == None:
                    # set name of fragment
                    self.fraglist[fn] = ft
                else:
                    # check if this is the same name
                    assert self.fraglist[fn] == ft, \
                         "The fragmentname %s of atom %d does not match with the prior definition %s" % (ft, fn, self.fraglist[fn]) 
                if ft not in self.fragnames:
                    self.fragnames.append(ft)
                self.frag_atoms[self._mol.fragnumbers[i]].append(i)
        # in the end make sure that all fragments have been named
        if None in self.fraglist:
            raise ValueError, "A fragment name is missing"
        return

    def make_frag_conn(self):
        """
        generate a fragment connectivity
        """
        assert self.setup
        # prepare the atom list for the fragments
        self.frag_conn = []
        self.frag_conn_atoms = []
        for i in xrange(self.nfrags):
            self.frag_conn.append([])
            self.frag_conn_atoms.append([])
        for i,f in enumerate(self.fraglist):
            # determine all external bonds of this fragment
            for ia in self.frag_atoms[i]:
                for ja in self._mol.conn[ia]:
                    j = self._mol.fragnumbers[ja]
                    if i != j:
                        # this is an external bond
                        self.frag_conn[i].append(j)
                        self.frag_conn_atoms[i].append((ia,ja))
                        logger.debug("fragment %d (%s) bonds to fragment %d (%s) %s - %s" %\
                                      (i, f, j, self.fraglist[j], self._mol.atypes[ia], self._mol.atypes[ja]))
        return
            
            


