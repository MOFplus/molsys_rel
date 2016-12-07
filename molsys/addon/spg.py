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
        self.symprec = 1.0e-5
        return
        
    def set_symprec(self, thresh):
        self_symprec = thresh
        return

    def generate(self, omit_H=False):
        lattice = self._mol.get_cell()
        pos     = self._mol.get_frac_xyz()
        num     = self._mol.get_elems_number()        
        self.spgcell = (lattice, pos, num)
        return

    def get_spacegroup(self):
        result = spglib.get_spacegroup(self.spgcell, symprec=self.symprec)
        return result                  
        
        
