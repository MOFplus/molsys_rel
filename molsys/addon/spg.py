# -*- coding: utf-8 -*-
"""

    spg

    implements an addon to access the features of the spglib within molsys
    https://atztogo.github.io/spglib/

    you need a recent spglib (>= 1.9.0) because the python import has changed at some point
    it can be installed via pip (pip install spglib)


Created on Wed Dec  7 15:44:36 2016

@author: rochus
"""

import spglib
import numpy
from molsys.util import elems

class spg:

    def __init__(self, mol):
        """
        generate a spg object

        :Parameters:

            - mol: mol object to be kept as a parent ref
        """
        self._mol = mol
        self.spgcell = None # tuple of (lattice, position, numbers) as used in spglib
        #
        self.spg_version = spglib.get_version()
        self.symprec = 1.0e-4
        return

    def set_symprec(self, thresh):
        self.symprec = thresh
        return

    def generate(self, omit=[]):
        # convert element symbols into atomic numbers
        if len(omit)>0:
            new_omit = []
            for e in omit:
                if type(e) != type(1):
                    new_omit.append(elems.number[e.lower()])
                else:
                    new_omit.append(e)
            omit = new_omit
        lattice = numpy.array(self._mol.get_cell(), order="C", dtype="double")
        pos     = self._mol.get_frac_xyz()
        pos     = pos%1.0
        # pos     = pos%1.0
        num     = self._mol.get_elems_number()
        pos_rem = []
        num_rem = []
        for i in xrange(self._mol.natoms):
            if num[i] not in omit:
                pos_rem.append(pos[i])
                num_rem.append(num[i])
        pos_rem = numpy.array(pos_rem, order="C", dtype="double")
        num_rem = numpy.array(num_rem, dtype="intc")
        self.spgcell = (lattice, pos_rem, num_rem)
        return

    def get_spacegroup(self):
        result = spglib.get_spacegroup(self.spgcell, symprec=self.symprec)
        return result


