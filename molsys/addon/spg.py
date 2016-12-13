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
import molsys

import logging
logger = logging.getLogger("molsys.spg")

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
        self.symprec = 1.0e-2
        logger.info("Addon spg loaded (version %d.%d.%d)" % self.spg_version)
        return

    def set_symprec(self, thresh):
        """
        set the symmetry threshold
        """
        self.symprec = thresh
        return

    def get_symprec(self, thresh):
        """
        get the symmetry threshold
        """
        return self.symprec

    def generate(self, omit=[]):
        """
        Generate the spglib specific representation of the structure
        (Needs to be called before any other call to spg methods)
         :Parameters:

             - omit : a list of either integers (atomic numbers) or element strings to be omited [optional]
        """
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
        """
        determine the space group of the current system
        returns a tuple with the symbol and the integer number
        """
        assert self.spgcell != None
        result = spglib.get_spacegroup(self.spgcell, symprec=self.symprec)
        result = result.split()
        symbol = result[0]
        number = int(result[1][1:-1])
        return (symbol, number)

    def make_P1(self, spgnum):
        """
        to be implemented by Julian from his topo tools

        :Parameters:

            - spgnum : integer space group number
        """
        # how to convert international spgnum to hall number
        dataset = spglib.get_symmetry_from_database(spgnum)
        # apply operations to self._mol and generate a new mol object
        # use detect_conn etc to complete it.
        return

    def get_primitive_cell(self):
        """
        get the primitve cell as a new mol object
        """
        assert self.spgcell != None
        new_spgcell = spglib.find_primitive(self.spgcell)
        if new_spgcell == None:
            logger.error("Search for primitive cell failed with symprec %f" % self.symprec)
            return
        print new_spgcell[0]
        print new_spgcell[2]
        new_mol = molsys.mol()
        new_mol.set_natoms(len(new_spgcell[2]))
        new_mol.set_cell(new_spgcell[0])
        new_mol.set_xyz(new_mol.get_real_from_frac(new_spgcell[1]))
        new_mol.set_elems_number(new_spgcell[2])
        # now add the connectivity
        new_mol.detect_conn()
        # RS: we could do atomtyping ... but this would have to be a method of mol ...
        new_mol.set_atypes(["0"]*new_mol.get_natoms())
        new_mol.set_nofrags()
        return new_mol


