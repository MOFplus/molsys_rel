# -*- coding: utf-8 -*-

#from molsys import *
import string
import copy
import mol
import molsys.util.images as images
from io import formats
import numpy as np

class bb(mol.mol):

    def __init__(self):
        mol.mol.__init__(self)
        self.dummies_hidden=False
        self.connectors = []
        self.connector_dummies=[]
        self.connecting_atoms = []
        return

    def write(self,fname,ftype='mfpx',**kwargs):
        ''' generic writer for the mol class
        :Parameters:
            - fname        : the filename to be written
            - ftype="mfpx" : the parser type that is used to writen the file
            - **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        formats.write[ftype](self,fname,bb=True,**kwargs)
        return

    def setup(self,name='default',specific_conn=None, linker=False, zflip=False, nrot=2, label = None):
        self.specific_conn = specific_conn  # that should be obtained from the file itself ?!!?
        self.linker = linker
        self.name = name
        self.zflip  = zflip
        self.nrot   = nrot
        if not linker:
            if self.zflip: print "Warning: zflip only supported for linkers"
            if self.nrot>1: print "Warning: rotations only supported for linkers"
        if linker: self.rotate_on_z()
        self.label = label
        #self.find_dummies()
        self.center()
        self.extract_connector_xyz()
        self.hide_dummy_atoms()
        return

    def center(self):
        if self.center_point == "com":
            # compute center of mass
            amass = []
            for e in self.elems: amass.append(elements.mass[e])
            amass = np.array(amass,"d")
            molmass = np.sum(amass)
            center = np.sum((self.xyz*amass[:,np.newaxis]),axis=0)
            center /= molmass
        elif self.center_point == "coc":
            # compute center of connectors (all equal mass)
            convec = []
            for c in self.connectors: convec.append(self.xyz[c])
            center = np.sum(np.array(convec),axis=0)
            center /= float(len(self.connectors))
        elif self.center_point == "special":
            center = self.special_center_point
        else:
            print "unknown center point option"
            raise IOError
        #self.center_xyz = center
        self.translate(-center)
        return
    
    #def find_dummies(self,dummy_label='x'):
        

    def hide_dummy_atoms(self):
        self.dummies_hidden=True
        self.bb = copy.deepcopy(self)
        self.natoms = self.natoms - len(self.connector_dummies)
        self.xyz = self.xyz[0:self.natoms,:]
        self.conn = self.conn[0:self.natoms]
        self.elems = self.elems[0:self.natoms]
        self.atypes =self.atypes[0:self.natoms]
        return
    
    def extract_connector_xyz(self):
        conn_xyz = []
        self.conn_elems = []
        for c in self.connectors:
            conn_xyz.append(self.xyz[c].tolist())
            self.conn_elems.append(self.elems[c])
        self.connector_xyz = np.array(conn_xyz,"d")
        self.conn_dist = np.sqrt(np.sum(self.connector_xyz*self.connector_xyz,axis=1))



    def rotate_on_z(self):
        """ especially if this is a linker (2 connectors) we want it to lie on the z-axis
        do this AFTER center but BEFORE extract_connector_xyz
        we always use the first connector (could also be a regular SBU!) to be on the z-axis """
        c1_xyz = self.xyz[self.connectors[0]]
        z_axis = np.array([0.0,0.0,1.0],"d")
        theta = rotations.angle(z_axis,c1_xyz) # angle to rotate
        if (theta > 1.0e-10) and (theta < (np.pi-1.0e-10)):
            axis  = rotations.normalize(rotations.cross_prod(z_axis,c1_xyz)) # axis around which we rotate
            self.xyz = rotations.rotate(self.xyz, axis, -theta)
        return

    def get_coc(self,conns):
        amass = []
        for e in conns: amass.append(1.0)
        amass = np.array(amass,dtype='float64')
        #print amass, amass[:,np.newaxis].shape,self.xyz.shape
        conns=np.array(conns,dtype='int')
        center = np.sum((amass[:,np.newaxis]*self.xyz[conns,:]),axis=0)/np.sum(amass)
        #center = np.sum((amass[:,np.newaxis]*self.xyz),axis=0)/np.sum(amass)
        return center


    def is_superpose(self, other, thresh=1.0e-1):
        """ we test if two molecular systems are equal (superimpose) by way of calculating the rmsd
        :Parameters:
            - other      : mol instance of the system in question
            - thresh=0.1 : allowed deviation of rmsd between self and other mol
        """
        if self.natoms != other.natoms: return False
        rmsd = 0.0
        for i in xrange(self.natoms):
            sxyz = self.xyz[i]
            r = other.xyz-sxyz
            d = np.sqrt(np.sum(r*r, axis=1))
            closest = np.argsort(d)[0]
            if d[closest] > thresh: return False, 0.0
            if self.elems[i] != other.elems[closest]: return False, 0.0
            rmsd += d[closest]
        rmsd = np.sqrt(np.sum(rmsd*rmsd))/self.natoms
        return True, rmsd

#### file conversion stuff

    def calc_centers(self,shift_to_original=True):
        centers = []
        for i,d  in enumerate(self.connecting_atoms):
            ci = np.sum(self.xyz[d],axis=0)/float(len(d))
            centers.append(ci)
        if shift_to_original==True:
            self.centers = centers+self.center_xyz
        else:
            self.centers = centers
        return centers



    #def connectors_to_dummy_neighbors(self):




