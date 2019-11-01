# -*- coding: utf-8 -*-
from __future__ import absolute_import

import molsys.util.elems as elements
import molsys.util.rotations as rotations
import string
import copy
from collections import Counter
import numpy as np


import logging
logger = logging.getLogger("molsys.bb")

class bb:

    def __init__(self,mol):
        """        '''
        initialize data structure for mol object so it can handle building block files

        Args:
            mol (molsys.mol): mol instance where this addon is added to
        """
        self.mol = mol
        assert self.mol.bcond == 0, "BBs must be non-periodic. mol object has bcond=%d" % self.mol.bcond
        self.connector = []            # list of connector atoms (can be dummies) used for orientation (TBI: COM of multiple atoms)
        self.connector_atoms = []       # list of lists: atoms that actually bond to the other BB
        self.connector_types = []       # list of integers: type of a connector (TBI: connectors_atype should contain the atype of the OTHER atom bonded to .. same layout as connector_atoms)
        self.connector_dummies = []       # list of integers: index of connectors that are dummies
        self.connector_atypes = None
        self.center_type = None
        self.mol.is_bb = True
        return

    def __mildcopy__(self, memo):
        """self.mol instance is kept the same, the rest is deepcopied
        __mildcopy__ is meant as an auxiliary method of mol.__deepcopy__
        to prevent recursion error. Deepcopying a bb instance works
        as usual because the bb.mol.__deepcopy__ stops the recursion
        with bb.mol.bb.__mildcopy__"""
        try: #python3
            newone = type(self)(self.mol.__class__())
        except: #python2
            newone = type(self)(bb, None)
        newdict = newone.__dict__
        newdict.update(self.__dict__)
        for key, val in newdict.items():
            if key != "mol":
                newdict[copy.deepcopy(key, memo)] = copy.deepcopy(val, memo)
        return newone

    def setup(self, connector,
                    connector_atoms=None,
                    connector_types= None,
                    connector_atypes=None,
                    center_type="coc",
                    center_xyz=None,
                    name=None,
                    rotate_on_z=False,
                    align_to_pax=False,
                    pax_use_connonly=False):
        """setup the BB data (only use this method to setup BBs!)

        TBI: possibility to have a sequence of indices for a connector -> set a flag
        
        Args:
            connector (list of ints): list of connector atoms
            connector_atoms (list of lists of ints, optional): actual atoms bonding. Defaults to None.
            connector_types (list of ints, optional): if present then special connectors exist. Defaults to None.
            connector_atypes (list of lists of strings, optional): atomtype of the atom to which the connector is bonded to
                                                 same layout as connector_atoms, if present then _types are generted. Defaults to None.
            center_type (string, optional): either "com" or "coc" or "special". if "special" center_xyz must be given. Defaults to "coc".
            center_xyz (numpy array, optional): coordinates of the center. Defaults to None.
            name (string, optional) : The name of the BB
        """
        assert not (rotate_on_z and align_to_pax) 
        self.connector = connector
        nc = len(connector)
        if connector_atoms is not None:
            assert len(connector_atoms)==nc
            self.connector_atoms = connector_atoms
        else:
            self.connector_atoms = [[i] for i in self.connector]
        if connector_types is not None:
            assert len(connector_types)==nc
            self.connector_types = connector_types
        else:
            if connector_atypes is not None:
                # if atypes are given then no types should be given .. this is determined
                assert len(connector_atypes)==nc
                for a,at in zip(self.connector_atoms,connector_atypes):
                    assert len(a)==len(at)
                self.connector_atypes = connector_atypes
                self.connector_types = []
                known_atypes = []
                for at in self.connector_atypes:
                    if at not in known_atypes:
                        known_atypes.append(at)
                    self.connector_types.append(known_atypes.index(at))
            else:
                self.connector_types = [0 for i in range(nc)]
        # detect connector dummies by element 'x'
        for i,cats in enumerate(self.connector_atoms):
            for j,cat in enumerate(cats):
                if self.mol.elems[cat-1].lower() == 'x':
                    self.connector_dummies.append(cat-1)
        assert center_type in ["com", "coc", "special"]
        self.center_type = center_type
        if self.center_type == "special":
            assert center_xyz is not None
            assert center_xyz.shape == (3,)
        elif self.center_type == "com":
            self.mol.set_real_mass()
            center_xyz = self.mol.get_com()
        else:
            mass, masstype = self.mol.get_mass(return_masstype=True)
            self.mol.set_unit_mass()            
            center_xyz = self.mol.get_com(idx=self.connector)
            self.mol.set_mass(mass, masstype)
        self.mol.translate(-center_xyz)
        if name is not None:
            self.name = name
        if rotate_on_z:
            self.rotate_on_z()
        if align_to_pax:
            self.align_pax_to_xyz(use_connxyz=pax_use_connonly)
        # finally sort for rising type in order to allow proper wrting to mfpx
        self.sort_connector_type()
        return

    def setup_with_bb_info(self,conn_identifier = 'He',center_point='coc'):
        """ Converts a mol object with a given atom as conn_identifier label to a BB.
        It can currently only work with a single conn_identifier, which means, that no special_connectors can be defined in this way.        
        Args:
            conn_identifier (str, optional): Atom type of the Atom used as connector label The connected atoms of the conn_identified atom will become the new connectors. Defaults to 'He'.
            center_point (str, optional): Definition of which center to use in the BB. Defaults to 'coc'.
        """
        # get indices of atoms conencted to conn_identifier
        cident_idx = [i for i,e in enumerate(self.mol.elems) if e.lower() == conn_identifier.lower()]
        logger.debug(cident_idx)
        connector = []
        for i,c in enumerate(self.mol.conn):
            for j,ci in enumerate(c):
                if cident_idx.count(ci) != 0: # atom i is connected to to an identifier_atom
                    connector.append(i)
        logger.debug('connectors',self.mol.connectors)
        # remove identifier atoms
        for ci,c in enumerate(connector):
            offset = [True for cidx in cident_idx if cidx < c].count(True)
            connector[ci] -= offset
        self.mol.delete_atoms(cident_idx)
        self.setup(connector)
        return
    
    @property
    def connector_xyz(self):
        """ get the xyz coords of the connector atoms
        TBI: if a connector is two atoms compute the center of mass """
        return self.mol.get_xyz(idx=self.connector)

    @property
    def connector_dist(self):
        """ get the distance from the center (origin 0,0,0) to the connector positions
        """
        return np.linalg.norm(self.connector_xyz, axis=1)

    def rotate_on_z(self):
        """ especially if this is a linker (2 connectors) we want it to lie on the z-axis
        we always use the first connector (could also be a regular SBU!) to be on the z-axis """
        c1_xyz = self.mol.xyz[self.connectors[0]]
        z_axis = np.array([0.0,0.0,1.0],"d")
        theta = rotations.angle(z_axis,c1_xyz) # angle to rotate
        if (theta > 1.0e-10) and (theta < (np.pi-1.0e-10)):
            axis  = rotations.normalize(rotations.cross_prod(z_axis,c1_xyz)) # axis around which we rotate
            self.mol.xyz = rotations.rotate(self.mol.xyz, axis, -theta)
        return

    def align_pax_to_xyz(self,use_connxyz = False):
        """ Aligns the coordinates to match the three principal axes of the intertial tensor to the cartesian coordinate system.
        Does not yet use weights
        
        Args:
            use_connxyz (bool, optional): If True, use only the coordinates of the connectors. Defaults to False.
        """
        ### TODO has been tested for m-bdc and some others to work, test for others aswell! 
        if use_connxyz == False:
            xyz = self.mol.get_xyz()
            self.mol.set_xyz(rotations.align_pax(xyz))
            return
        xyz = self.connector_xyz
        eigval,eigvec = rotations.pax(xyz)
        eigorder = np.argsort(eigval)
        rotmat = eigvec[:,eigorder] #  sort the column vectors in the order of the eigenvalues to have largest on x, second largest on y, ... 
        self.mol.set_xyz(rotations.apply_mat(rotmat,self.mol.get_xyz()))
        return

    def sort_connector_type(self):
        """This internal method sorts the connectors according to a rsing type

        this is necessary for writing mfpx files, where the different types are seperated by a slash
        """
        order = np.argsort(self.connector_types)
        self.connector      = [self.connector[i] for i in order]
        self.connector_atoms = [self.connector_atoms[i] for i in order]
        self.connector_types = [self.connector_types[i] for i in order]
        if self.connector_atypes is not None:
            self.connector_atypes = [self.connector_atypes[i] for i in order]
        return
