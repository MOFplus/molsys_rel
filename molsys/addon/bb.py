# -*- coding: utf-8 -*-
from __future__ import absolute_import

#from molsys import mol
import molsys.util.elems as elements
import molsys.util.rotations as rotations
import string
import copy
import numpy as np
numpy = np


import logging
logger = logging.getLogger("molsys.bb")

class bb:

    def __init__(self,mol):
        '''
        initialize data structure for mol object so it can handle building block files
        beware, there is a litte chaos in the writing of bb files. what is needed is:
        - is_mol=True i.o. to write
        - mol.connectors as a list of 'connectors'
        - mol.connector_atoms as a list of lists of the atoms (inner list) belonging to the connector (outer list)
        - mol.connectors_type a list if type indices, nconn*[0] for no special connectors, otherwise something like [0,0,0,0,1,1]. needs to be ordered ascendingly:w
        '''
#        if mol.is_bb == False:
#            raise IOError,("No bb info available!")
        self.mol = mol
        self.mol.dummies_hidden=False
        self.mol.connectors = []
        self.mol.connector_dummies=[]
        self.mol.connector_atoms = []
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



    def setup(self,name='default',specific_conn=None, linker=False, zflip=False, nrot=2, label = None):
        self.mol.specific_conn = specific_conn  # that should be obtained from the file itself ?!!?
        self.mol.linker = linker
        self.mol.name = name
        self.mol.zflip  = zflip
        self.mol.nrot   = nrot
        if not linker:
            if self.mol.zflip: logger.warning("zflip only supported for linkers")
            if self.mol.nrot>1: logger.debug("rotations only supported for linkers")
        if linker: self.rotate_on_z()
        self.mol.label = label
        #self.find_dummies()
        self.center()
        self.extract_connector_xyz()
#        self.hide_dummy_atoms()
        return

    def center(self):
        if self.mol.center_point == "com":
            self.mol.set_real_mass()
            center = self.mol.get_com()
        elif self.mol.center_point == "coc":
            self.mol.set_unit_mass()
            center = self.mol.get_com(idx=self.mol.connectors)
        elif self.mol.center_point == "special":
            center = self.mol.special_center_point
        else:
            raise IOError("unknown center point option")
        self.mol.translate(-center)
        return

    def get_coc(self):
        """redundant if center_point == 'coc' """
        mass, masstype = self.mol.get_mass(return_masstype=True)
        self.mol.set_unit_mass()
        coc = self.mol.get_com(idx=self.mol.connectors)
        self.mol.set_mass(mass, masstype)
        return coc

    #def get_radius(self):
    #    coc = self.get_coc()
    #    self.mol.connector_xyz[:,0:]

    def hide_dummy_atoms(self):
        ''' depreciated, has been used to remove dummies, requires them to be the last atoms
            we now remove dummies after construction using remove_dummies() in mol.py'''
        self.mol.dummies_hidden=True
        self.mol.bb = copy.deepcopy(self.mol)
        self.mol.natoms = self.mol.natoms - len(self.mol.connector_dummies)
        self.mol.xyz = self.mol.xyz[0:self.mol.natoms,:]
        self.mol.conn = self.mol.conn[0:self.mol.natoms]
        self.mol.elems = self.mol.elems[0:self.mol.natoms]
        self.mol.atypes =self.mol.atypes[0:self.mol.natoms]
        return

    def extract_connector_xyz(self):
        conn_xyz = []
        self.mol.conn_elems = []
        for c in self.mol.connectors:
            conn_xyz.append(self.mol.xyz[c].tolist())
            self.mol.conn_elems.append(np.array(self.mol.elems)[c].tolist())
        try:
            self.mol.connector_xyz = np.array(conn_xyz,"d")
        except ValueError:
            conn_xyz = [np.mean(cc,axis=0) for cc in conn_xyz]
            self.mol.connector_xyz = np.array(conn_xyz,"d")
        self.mol.conn_dist = np.sqrt(np.sum(self.mol.connector_xyz*self.mol.connector_xyz,axis=1))



    def rotate_on_z(self):
        """ especially if this is a linker (2 connectors) we want it to lie on the z-axis
        do this AFTER center but BEFORE extract_connector_xyz
        we always use the first connector (could also be a regular SBU!) to be on the z-axis """
        c1_xyz = self.mol.xyz[self.mol.connectors[0]]
        z_axis = np.array([0.0,0.0,1.0],"d")
        theta = rotations.angle(z_axis,c1_xyz) # angle to rotate
        if (theta > 1.0e-10) and (theta < (np.pi-1.0e-10)):
            axis  = rotations.normalize(rotations.cross_prod(z_axis,c1_xyz)) # axis around which we rotate
            self.mol.xyz = rotations.rotate(self.mol.xyz, axis, -theta)
        return

    def is_superpose(self, other, thresh=1.0e-1):
        """ we test if two molecular systems are equal (superimpose) by way of calculating the rmsd
        :Parameters:
            - other      : mol instance of the system in question
            - thresh=0.1 : allowed deviation of rmsd between self and other mconnecl
        """
        if self.mol.natoms != other.natoms: return False
        rmsd = 0.0
        for i in range(self.mol.natoms):
            sxyz = self.mol.xyz[i]
            r = other.xyz-sxyz
            d = np.sqrt(np.sum(r*r, axis=1))
            closest = np.argsort(d)[0]
            if d[closest] > thresh: return False, 0.0
            if self.mol.elems[i] != other.elems[closest]: return False, 0.0
            rmsd += d[closest]
        rmsd = np.sqrt(np.sum(rmsd*rmsd))/self.mol.natoms
        return True, rmsd


    def calc_centers(self,shift_to_original=True):
        centers = []
        for i,d  in enumerate(self.mol.connecting_atoms):
            ci = np.sum(self.mol.xyz[d],axis=0)/float(len(d))
            centers.append(ci)
        if shift_to_original==True:
            self.mol.centers = centers+self.mol.center_xyz
        else:
            self.mol.centers = centers
        return centers
    
    def add_bb_info(self,conn_identifier = 'He',center_point='coc'):
        '''
        use for quick conversion of a chemdraw-built MM3 tinker txyz file into a BB file type
        the connecting atom is next to a He atom, it is beign deleted and 
        '''
        self.mol.center_point = center_point
        self.mol.is_bb=True
        
        # get indices of atoms conencted to conn_identifier
        cident_idx = [i for i,e in enumerate(self.mol.elems) if e.lower() == conn_identifier.lower()]
        logger.debug(cident_idx)
        self.mol.connectors = []
        for i,c in enumerate(self.mol.conn):
            for j,ci in enumerate(c):
                if cident_idx.count(ci) != 0: # atom i is connected to to an identifier_atom
                    self.mol.connectors.append(i)
        logger.debug('connectors',self.mol.connectors)
        # remove identifier atoms
        for ci,connector in enumerate(self.mol.connectors):
            offset = [True for cidx in cident_idx if cidx < connector].count(True)
            self.mol.connectors[ci] -= offset
        self.mol.connector_atoms = [[i] for i in self.mol.connectors]
        self.mol.connectors_type = [0 for i in range(len(self.mol.connectors))]
        self.mol.delete_atoms(cident_idx)
        self.align_pax_to_xyz(use_connxyz=False)
        return
    
    def align_pax_to_xyz(self,use_connxyz = False):
        '''
        use to align the principal axes of the building block with x,y,and z.
        does not yet use weights, but assumes w=1 for all  atoms
        '''
        ### TODO has been tested for m-bdc and some others to work, test for others aswell! 
        import molsys.util.rotations as rotations
        self.center()
        if use_connxyz == True:
            xyz = numpy.array(self.mol.xyz[self.mol.connectors])
        else:
            xyz = numpy.array(self.mol.xyz)
            self.mol.xyz = rotations.align_pax(xyz)
            #self.mol.view(program='vmdx')
            return
        eigval,eigvec = rotations.pax(xyz)
        eigorder = numpy.argsort(eigval)
        rotmat = eigvec[:,eigorder] #  sort the column vectors in the order of the eigenvalues to have largest on x, second largest on y, ... 
        self.mol.xyz = rotations.apply_mat(rotmat,self.mol.xyz)
        #import pdb; pdb.set_trace()
        return


