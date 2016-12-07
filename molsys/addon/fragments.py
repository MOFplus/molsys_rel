# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 21:30:26 2016

@author: rochus
"""

import logging
logger = logging.getLogger("molsys.graph")


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
        self.nfrags =  max(self._mol.fragnumbers)+1               
        self.fraglist  = [None]*(self.nfrags) # additional list to hold the fragments with their name
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
        self.frag_atoms = []
        self.frag_conn = []
        for i in xrange(self.nfrags):
            self.frag_atoms.append([])
        # now sort in the atoms
        for i in xrange(self._mol.natoms):
            self.frag_atoms[self._mol.fragnumbers[i]].append(i)
        for i,f in enumerate(self.fraglist):
            # determine all external bonds of this fragment
            ext_bond = []
            fraga = self.frag_atoms[i]
            for ia in fraga:
                for j in self._mol.conn[ia]:
                    if j not in fraga:
                        # this is an external bond
                        ext_bond.append(j)
            logger.debug("fragment %s has %d external bonds" % (f, len(ext_bond)))
            # now check to which fragments these external bonds belong to
            # we connect both ways and test consistency at the end
            # if a fragment has multiple bonds to another this is ignored 
            fconn = []
            for ea in ext_bond:
                for j,fj in enumerate(self.fraglist):
                    if ea in self.frag_atoms[j]:
                        if j not in fconn:
                            logger.debug(" -> bonded to fragment %s" % fj)
                            fconn.append(j)
                        break
            self.frag_conn.append(fconn)
        # now test if the fragemnt connectivity is two way consistent
        for i in xrange(self.nfrags):
            for j in self.frag_conn[i]:
                if j>i:
                    if not i in self.frag_conn[j]:
                        print "Fragment topology is inconsistent! This should never happen"
                        print "%d is bonded to %d but not reverse" % (i, j)
                        raise ValeError, "DEBUG this code!"
        return
            



