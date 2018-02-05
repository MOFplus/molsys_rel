#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from molsys.addon import base

class hydrogen(base):
    """
    Class offering methods to add hydrogens to experimental structures. Up to now
    only a method for sp2 carbons is implemented.
    """

    def add_sp2(self, idx):
        """
        Method to add a hydrogen to a sp2 carbon.

        :Parameters:
          - idx (int): index of the carbon atom where the  hydrogen should be added.
        """
        # get central atom
        ca = self._mol.xyz[idx]
        # check number of bonded partners
        assert len(self._mol.conn[idx]) == 2
        # get bonded atom pos, map them to the image of the first one
        mapped = self._mol.map2image(self._mol.xyz[[idx, self._mol.conn[idx][0], self._mol.conn[idx][1]]])
        vecs = mapped[1:]-ca
        vec = np.sum(vecs, axis=0)
        nvec = -vec/np.linalg.norm(vec)
        # add hydrogen
        self._mol.add_atom('h','h', ca+1.098*nvec)
        self._mol.add_bonds(idx, self._mol.natoms - 1)
        return
