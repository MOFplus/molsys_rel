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
    
    def __init__(self, mol, orig_atoms=None, cell_factor = 1.5):
        """Generate a slicer object for a given mol obejct
        
        Args:
            mol (mol object): periodic system to be sliced
            orig_atoms (list of integers, optional): Defaults to None. COM of these atoms defines the origin of the cell  
            cell_factor (float, optional): Defaults to 1.5. [description]
        """
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
        self.cell_factor = cell_factor
        return
    
    def set_stub(self, atype, new_elem, new_atype, dist):
        """set a new stub definition (used when slicing)
        
        Args:
            atype (string): atom type that remains to which a bond has been sliced
            new_elem (string): Element of the stub
            new_atype (string): atomtype of the stub
            dist (float): distance from the remaining element
        """
        self.stubs[atype] = (new_elem, new_atype, dist)
        return
    
    def set_plane(self, hkl, dist=None, symm=False):
        """set a new plane to cut
        
        Args:
            hkl (list of integers): defines the hkl plane to be sliced
            dist (float, optional): Defaults to None. distance in multiples of cell params from origin where to cut
            symm (bool, optional): Defaults to False. if True a second plane at negative dist is applied
        """
        hkl = np.array(hkl)
        assert hkl.shape == (3,)
        if dist is None:
            dist=1.0
        cell = self.mol.get_cell()
        vect = np.dot(cell, hkl)
        self.planes.append(vect*dist)
        if symm:
            self.planes.append(vect*-dist)
        return
    
    def __call__(self, supercell=None, copy=False, orig=[0,0,0], max_dist=2.0):
        """perform the slicing operation
            supercell (list of three integers, optional): Defaults to None. size of the supercell of the initial system to be generated before performing the slicing
            copy (bool, optional): Defaults to False. If True a new mol obejct is generated and returned (by default it is modified)
            orig (list of integers, optional): Defaults to [0,0,0]. origin in sizes of the initial cell
            max_dist (float, optional): Defaults to 2.0. maximum distance of a bond for setting a stub
        
        Returns:
            [type]: [description]
        """
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
        periodic = np.array([True,True,True])
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
            # no test which dimension this plane cuts
            periodic = np.logical_and(np.equal(p, 0.0),periodic)
        # now add stubs (only those registered)
        for s in stubs:
            if s[0] not in delete:
                at = mol.get_atypes()[s[0]]
                if at in self.stubs:
                    sel, sat, sdist = self.stubs[at]
                    # to compute the new position we need sdist times the unit vect plus the overall shift
                    spos = xyz[s[0]]+sdist*s[1]+shift
                    k = mol.add_atom(sel, sat, spos)
                    # print("adding bond between %d %s and %d %s" % (s[0], self.mol.elems[s[0]], k, self.mol.elems[k]))
                    mol.add_conn(s[0], k)
        # now all is parsed and we can remove the corresponding atoms
        mol.delete_atoms(delete) 
        # now add an offset to the cell for all nonperiodic dims
        cell = self.mol.get_cell()
        for i,d in enumerate(periodic):
            if not d:
                cell[i]*= self.cell_factor
        self.mol.set_cell(cell)
        # now regenerate the pconn
        self.mol.add_pconn()
        return mol
    
    
        
        
