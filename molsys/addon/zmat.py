# -*- coding: utf-8 -*-

"""

       module to implement an addon feature: chemcoord zmatrix manipulation 

       NOTE: this is only imported by __init__.py if chemcoords is present

       based on Marco Dygas routines in the old IOmod

"""
import chemcoord

import logging
import copy
import pandas
import numpy
import molsys.util.rotations as rotations
import scipy.optimize as opt
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
        xyz = cartesian[:, ['x', 'y', 'z']].as_matrix()
        idx = list(cartesian.index-1)
        ordered_xyz = []
        for i in range(self._mol.natoms):
            ordered_xyz.append(xyz[idx.index(i),:])
        self._mol.xyz = numpy.array(ordered_xyz)
        return
    
    def rotate_dihedral(self, idx, deg):
        """ 
        Rotates a dihedral angle
        Parameters:
          - idx : List of atom indices of the atoms spanning the dihedral
          - deg : target angle in degrees
        """
        if self._mol.xyz.shape[0] < 4:
            raise IOError('The amount of atoms in the molecule is smaller than 4!')
        if len(idx) != 4:
            raise IOError('The amount of indices is not 4!')
        xyz = self.to_Cartesian()
        idx = (numpy.array(idx)+1).tolist()
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

    def change_angle(self, idx, deg):
        """ 
        Changes the value of an angle
        Parameters:
          - idx : List of atom indices of the atoms spanning the angle
          - deg : target angle in degrees
        """
        if self._mol.xyz.shape[0] < 3:
            raise IOError('The amount of atoms in the molecule is smaller than 3!')
        if len(idx) != 3:
            raise IOError('The amount of indices is not 3!')
        xyz = self.to_Cartesian()
        idx = (numpy.array(idx)+1).tolist()
        idx_array = [[idx[0],      0,      0,      0], \
                     [idx[1], idx[0],      0,      0], \
                     [idx[2], idx[1], idx[0],      0]]
        idx_array = numpy.array(idx_array)
        buildlist = xyz._get_buildlist(fixed_buildlist = idx_array)
        zmat = xyz.to_zmat(buildlist)
        zmat[idx[2], 'angle'] = deg
        xyz = zmat.to_xyz()
        self.from_Cartesian(xyz)
        return

    def add_fragment(self, amol, pc, dist):
        def bound(triple):
            triple=numpy.array(triple,dtype='float64')
            floor = float(numpy.floor(triple[0])) % 2.0
            triple[0:3] %= 1.0
            if floor >= 0.5: triple[0] = 1.0 - triple[0]
            return triple
        def penalty(t, args):
            t = bound(t)
            n1,n2,c1,c2 = args[0],args[1],args[2],args[3]
            c2n = rotations.rotate_by_triple(c2,t)
            n2n = rotations.rotate_by_triple(n2,t)
            return 1-abs(numpy.dot(n1,n2n))+(1+numpy.dot(c1,c2n))
        natoms = self._mol.natoms
        ### mol vecs
        a1vec = self._mol.xyz[pc[1],:] - self._mol.xyz[pc[0],:]
        b1vec = self._mol.xyz[pc[2],:] - self._mol.xyz[pc[0],:]
        a1vec /= numpy.linalg.norm(a1vec)
        b1vec /= numpy.linalg.norm(b1vec)
        n1vec = numpy.cross(a1vec,b1vec)
        n1vec /= numpy.linalg.norm(n1vec)
        c1vec = -1.0 * (a1vec+b1vec)
        c1vec = c1vec/numpy.linalg.norm(c1vec)
        ### amol vecs
        a2vec = amol.xyz[amol.conn[0][0]] - amol.xyz[0]
        b2vec = amol.xyz[amol.conn[0][1]] - amol.xyz[0]
        a2vec /= numpy.linalg.norm(a2vec)
        b2vec /= numpy.linalg.norm(b2vec)
        c2vec = -1.0 * (a2vec+b2vec)
        c2vec = c2vec/numpy.linalg.norm(c2vec)
        n2vec = numpy.cross(a2vec,b2vec)
        n2vec /= numpy.linalg.norm(n2vec)
        ### rots
        count = 0
        while count < 1000:
            a = opt.minimize(penalty, numpy.random.uniform(0,1,(3,)),
            args = [n1vec, n2vec, c1vec, c2vec], method="SLSQP")
            if a.fun < 0.1:
                count = 100000000000
            else:
                print(a)
        ### trans
        c1vec = c1vec*dist
        trans = self._mol.xyz[pc[0],:]+c1vec
        self._mol.add_mol(amol,translate=trans, rotate = bound(a.x))
        self._mol.conn[pc[0]].append(natoms)
        self._mol.conn[natoms].append(pc[0])
        return


