# -*- coding: utf-8 -*-

"""

       module to implement an addon feature: chemcoord zmatrix manipulation 

       NOTE: this is only imported by __init__.py if chemcoords is present

       based on Marco Dygas routines in the old IOmod

"""
import chemcoord

import logging
import numpy
logger = logging.getLogger("molsys.zmat")

class zmat:

    def __init__(self, mol):
        """
        instantiate a graph object which will be attached to the parent mol

        :Parameter:

             - mol : a mol type object (can be a derived type like bb or topo as well)
        """
        self._mol = mol
        logger.debug("generated the zmat object")
        return

    def to_Cartesian(self):
        """
        transforms the mol object into a chemcoord.Cartesian object
        """
        natoms = self._mol.natoms
        elems = copy.deepcopy(self._mol.elems)
        for i, j in enumerate(elems): 
            elems[i] = j.strip().capitalize()
        xyz = pandas.DataFrame(self._mol.xyz, columns=["x","y","z"], dtype='float64')
        elems = pandas.DataFrame(elems, columns=["atom"], dtype='str')
        output = chemcoord.Cartesian(pandas.concat([elems, xyz], axis=1))
        output.index = range(1, natoms+1)
        return output
        
    def from_Cartesian(self, cartesian):
        """
        loads mol xyz data from a chemcoord.Cartesian object
        """
        #self._mol.xyz = cartesian[:, ['x', 'y', 'z']].as_matrix()
        self._mol.set_xyz(cartesian[:, ['x', 'y', 'z']].as_matrix())
        #elements = cartesian[:, 'atom'].as_matrix()
        #for i in range(len(elements)):
        #    if len(elements[i]) == 1:
        #        elements[i] += ' '
        #self.natoms = self.xyz.shape[0]
        #self.elements = elements
        #self.set_masses()
        #self.set_vdwr()
        return
    
    def rotate_dihedral(self, idx, deg):
        """ 
        Rotates a dihedral angle
        Parameters:
          - idx : List of atom indices of the atoms spanning the dihedral
          - deg : target angle in degrees
        """
        if self.xyz.shape[0] < 4:
            raise IOError('The amount of atoms in the molecule is smaller than 4!')
        if len(idx) != 4:
            raise IOError('The amount of indices is not 4!')
        xyz = self.to_Cartesian()
        idx_array = [[idx[0],      0,      0,      0], \
                     [idx[1], idx[0],      0,      0], \
                     [idx[2], idx[1], idx[0],      0], \
                     [idx[3], idx[2], idx[1], idx[0]]]
        idx_array = numpy.array(idx_array)
        buildlist = xyz._get_buildlist(fixed_buildlist = idx_array)
        zmat = xyz.to_zmat(buildlist)
        zmat[idx[3], 'dihedral'] = deg
        xyz = zmat.to_xyz()
        self.from_Cartesian(xyz)
        return

    def add_fragment(self, amol, pc, dist, dihed):
        natoms = self._mol.natoms
        avec = self._mol.xyz[pc[1],:] - self._mol.xyz[pc[0],:]
        bvec = self._mol.xyz[pc[2],:] - self._mol.xyz[pc[0],:]
        avec /= numpy.linalg.norm(avec)
        bvec /= numpy.linalg.norm(bvec)
        nvec = numpy.cross(avec,bvec)
        cvec = -1.0 * (avec+bvec)
        cvec = (cvec/numpy.linalg.norm(cvec))*dist
        trans = self._mol.xyz[pc[0],:]+cvec
        print trans
        self._mol.add_mol(amol,translate=trans)
        self._mol.conn[pc[0]].append(natoms)
        self._mol.conn[natoms].append(pc[0])
        return


