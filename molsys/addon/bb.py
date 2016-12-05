# -*- coding: utf-8 -*-

#from molsys import *
import molsys.util.elems as elements
import molsys.util.rotations as rotations
import string
import copy
import numpy as np


import logging
logger = logging.getLogger("molsys.bb")

class bb:

    def __init__(self,mol):
#        if mol.is_bb == False:
#            raise IOError, "No bb info available!"
        self.mol = mol
        self.mol.dummies_hidden=False
        self.mol.connectors = []
        self.mol.connector_dummies=[]
        self.mol.connector_atoms = []
        return

    def setup(self,name='default',specific_conn=None, linker=False, zflip=False, nrot=2, label = None):
        self.mol.specific_conn = specific_conn  # that should be obtained from the file itself ?!!?
        self.mol.linker = linker
        self.mol.name = name
        self.mol.zflip  = zflip
        self.mol.nrot   = nrot
        if not linker:
            if self.mol.zflip: logger.warning("zflip only supported for linkers")
            if self.mol.nrot>1: logger.warning("rotations only supported for linkers")
        if linker: self.rotate_on_z()
        self.mol.label = label
        #self.find_dummies()
        self.center()
        self.extract_connector_xyz()
        self.hide_dummy_atoms()
        return

    def center(self):
        if self.mol.center_point == "com":
            # compute center of mass
            amass = []
            for e in self.mol.elems: amass.append(elements.mass[e])
            amass = np.array(amass,"d")
            molmass = np.sum(amass)
            center = np.sum((self.mol.xyz*amass[:,np.newaxis]),axis=0)
            center /= molmass
        elif self.mol.center_point == "coc":
            # compute center of connectors (all equal mass)
            convec = []
            for c in self.mol.connectors: convec.append(self.mol.xyz[c])
            center = np.sum(np.array(convec),axis=0)
            center /= float(len(self.mol.connectors))
        elif self.mol.center_point == "special":
            center = self.mol.special_center_point
        else:
            print "unknown center point option"
            raise IOError
        #self.center_xyz = center
        self.mol.translate(-center)
        return


    def hide_dummy_atoms(self):
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
            self.mol.conn_elems.append(self.mol.elems[c])
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

    def get_coc(self,conns):
        amass = []
        for e in conns: amass.append(1.0)
        amass = np.array(amass,dtype='float64')
        #print amass, amass[:,np.newaxis].shape,self.xyz.shape
        conns=np.array(conns,dtype='int')
        center = np.sum((amass[:,np.newaxis]*self.mol.xyz[conns,:]),axis=0)/np.sum(amass)
        #center = np.sum((amass[:,np.newaxis]*self.xyz),axis=0)/np.sum(amass)
        return center


    def is_superpose(self, other, thresh=1.0e-1):
        """ we test if two molecular systems are equal (superimpose) by way of calculating the rmsd
        :Parameters:
            - other      : mol instance of the system in question
            - thresh=0.1 : allowed deviation of rmsd between self and other mol
        """
        if self.mol.natoms != other.natoms: return False
        rmsd = 0.0
        for i in xrange(self.mol.natoms):
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


