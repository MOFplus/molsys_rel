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


class plane(object):

    def __init__(self, name, hkl, dist):
        self.name = name
        self.hkl = np.array(hkl)
        self.dist = dist
        return

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
        if not orig_atoms is None:
            shift = np.zeros([3])
            for i in orig_atoms:
                shift += self.mol.get_xyz()[i]
            shift /= len(orig_atoms)
            self.mol.translate(-shift)
            self.mol.wrap_in_box()
        self.cell_factor = cell_factor
        return
    
    def set_stub(self, atype, new_elem, new_atype, dist, plane_name="plane"):
        """set a new stub definition (used when slicing)
        
        Args:
            atype (string): atom type that remains to which a bond has been sliced
            new_elem (string): Element of the stub
            new_atype (string): atomtype of the stub
            dist (float): distance from the remaining element
            plane_name (string, optional): defaults to "plane". Name of plane to use this stub for.
        """
        self.stubs[atype+"_"+plane_name] = (new_elem, new_atype, dist)
        return
    
    def set_plane(self, hkl, dist=0.5, name="plane", symm=False):
        """set a plane by hkl and dist, if symm is true a second plane with dist = -dist is added

        Args:
            hkl (list): hkl indices as a list
            dist (float, optional): Defaults to 0.5. Distance along this hkl vector in multiples 
            name (string, optional): Defaults to "plane". name of the plane
            symm (bool, optional): Defaults to False. add another plane at -dist if True
        """        
        hkl = np.array(hkl)
        assert hkl.shape == (3,)
        p = plane(name, hkl, dist)
        self.planes.append(p)
        if symm:
            p = plane(name, -hkl, dist)
            self.planes.append(p)
        return
    
    def __call__(self, supercell=None, copy=False, origshift=None, max_dist=2.0):
        """perform the slicing operation
            supercell (list of three integers, optional): Defaults to None. size of the supercell of the initial system to be generated before performing the slicing
            copy (bool, optional): Defaults to False. If True a new mol obejct is generated and returned (by default it is modified)
            origshift (list of floats, optional): Defaults to None. originshift in sizes of the initial cell
            max_dist (float, optional): Defaults to 2.0. maximum distance of a bond for setting a stub
        
        Returns:
            [type]: [description]
        """
        if copy:
            mol = self.mol.copy()
        else:
            mol = self.mol
        # the origin is by default at [0.5, 0.5, 0.5] in fractional coordinates
        # but we shift this by origin if this is provided. Note that this shift 
        # is with respect to the original!!! unit cell 
        orig = np.array([0.5, 0.5, 0.5],"d")
        if not (origshift is None):
            orig += origshift/np.array(supercell)
        if supercell:
            mol.make_supercell(supercell)
        delete = []
        stubs  = []
        xyz  = mol.get_xyz()
        fxyz = mol.get_frac_from_xyz()
        fxyz -= orig
        fxyz *= np.array(supercell)
        periodic = np.array([True,True,True])
        for ip, p in enumerate(self.planes):
            pname = p.name
            pvec  = p.hkl
            # run over all planes and find those atoms to be sliced 
            # at the same time detect those to be fixed with a stub
            normp = pvec/np.sqrt((pvec*pvec).sum())
            proj = np.dot(fxyz, normp)
            flags = np.greater(proj, p.dist)
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
                                    stubs.append((j, vect, pname))
            # no test which dimension this plane cuts
            periodic = np.logical_and(np.equal(pvec, 0.0),periodic)
        # now add stubs (only those registered)
        print self.stubs.keys()
        for s in stubs:
            if s[0] not in delete:
                at = mol.get_atypes()[s[0]]
                stub_name = at+"_"+s[2]
                print stub_name
                if stub_name in self.stubs:
                    sel, sat, sdist = self.stubs[stub_name]
                    # to compute the new position we need sdist times the unit vect plus the overall shift
                    spos = xyz[s[0]]+sdist*s[1]
                    k = mol.add_atom(sel, sat, spos)
                    print("adding bond between %d %s and %d %s" % (s[0], self.mol.elems[s[0]], k, self.mol.elems[k]))
                    mol.add_bonds(s[0], k)
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
    
    
        
        
