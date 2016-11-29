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
        return

    def check(self):
        """
        check if all atoms are in a fragment
        """
        self.fragnames = []
        self.setup = True
        for i in xrange(self._mol.natoms):
            ft = self._mol.fragtypes
            if ft == "0":
                logger.error("atom %d (%s) is not in fragment" % (i, self._mol.atypes[i]))
                self.setup=False
            else:
                if ft not in self.fragnames:
                    self.fragnames.append(ft)
        return

    def make_frag_conn(self):
        """
        generate a fragment connectivity
        """
        assert self.setup



        return

