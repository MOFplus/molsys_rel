#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 17:09:25 2017

        slicer
        
        class to take a mol object and truncate it by a number predefined slicing planes
        - can be used to make slabs 
        - or nano-particles of arbitrary shape

        planes are defined by a vector (to which these planes are orthogonal) and a distance from the origin 
        the vector is given in hkl indices and a relative distance.
         example: [1,0,0] is the 100 plane .. [2,0,0] the 200 plane
         but with [1,0,0] and a distance of 2 also the 200 plane can be defined

@author: rochus
"""

import molsys
import numpy as np

class slicer:
    
    def __init__(self, mol, orig_atoms=None):
        self.mol = mol
        self.planes = []
        self.stubs = {}
        # calibrate what is the origin of the initial system
        if orig_atoms:
            shift = np.zeros([3])
            for i in orig_atoms:
                shift += self.mol.get_xyz()[i]
            shift /= len(orig_atoms)
            self.mol.translate(-shift)
            self.mol.wrap_in_box()
        return
    
    def set_stub(self, atype, new_elem, new_atype, dist):
        self.stubs[atype] = (new_elem, new_atype, dist)
        return
    
    def set_plane(self, hkl, dist=None, symm=False):
        """ set a plane by hkl and dist, if symm is true a second plane with dist = -dist is added
        """
        hkl = np.array(hkl)
        assert hkl.shape == (3,)
        if dist == None:
            dist=1.0
        cell = self.mol.get_cell()
        vect = np.dot(cell, hkl)
        self.planes.append(vect*dist)
        if symm:
            self.planes.append(vect*-dist)
        return
    
    def __call__(self, supercell=None, copy=False, orig=[0,0,0], max_dist=2.0):
        if copy:
            mol = self.mol.copy()
        else:
            mol = self.mol
        shift = np.dot(mol.get_cell(), np.array(orig))
        if supercell:
            mol.make_supercell(supercell)
        delete = []
        stubs  = []
        xyz = mol.get_xyz() - shift
        for p in self.planes:
            # run over all planes and find those atoms to be sliced 
            # at the same time detect those to be fixed with a stub
            distp = np.sqrt((p*p).sum())
            normp = p/distp
            proj = np.dot(xyz, normp)
            flags = np.greater(proj, distp)
            for i,f in enumerate(flags):
                if f:
                    # atom i is outside the slicing panes and needs to be removed (if it is not alread to be deleted)
                    if delete.count(i) == 0:
                        delete.append(i)
                        # now test all atoms bonded to i if they are to be deleted as well
                        for j in mol.conn[i]:
                            if not flags[j]:
                                # compute normal vector from j to i 
                                vect = xyz[i]-xyz[j]
                                dvect = np.sqrt((vect*vect).sum())
                                vect = vect/dvect
                                if dvect < max_dist:
                                    stubs.append((j, vect))
        # now add stubs (only those registered)
        for s in stubs:
            if s[0] not in delete:
                at = mol.get_atypes()[s[0]]
                if at in self.stubs:
                    sel, sat, sdist = self.stubs[at]
                    # to compute the new position we need sdist times the unit vect plus the overall shift
                    spos = xyz[s[0]]+sdist*s[1]+shift
                    k = mol.add_atom(sel, sat, spos)
                    mol.add_conn(s[0], k)
        # now all is parsed and we can remove the corresponding atoms
        mol.delete_atoms(delete) 
        return mol
    
    
        
        
