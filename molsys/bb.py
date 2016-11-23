# -*- coding: utf-8 -*-

#from molsys import *
import string
import copy
import numpy as np
import molsys

class bb(molsys.molsys):
    
    def __init__(self, name, specific_conn=None, linker=False, zflip=False, nrot=2, label = None, filetype='mfpx'):
        molsys.molsys.__init__(self)
        self.specific_conn = specific_conn  # that should be obtained from the file itself ?!!?
        self.linker = linker
        self.zflip  = zflip
        self.nrot   = nrot
        if not linker:
            if self.zflip: print "Warning: zflip only supported for linkers"
            if self.nrot>1: print "Warning: rotations only supported for linkers"
        self.connectors = []
        self.dummies=[]
        self.connecting_atoms = []
        self.center()
        if linker: self.rotate_on_z()
        self.extract_connector_xyz()
        self.hide_dummy_atoms()
        self.name = name
        self.label = label
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
        self.center_xyz = center
        self.translate(-center)
        return
        
    def extract_connector_xyz(self):
        conn_xyz = []
        self.conn_elems = []
        for c in self.connectors:
            conn_xyz.append(self.xyz[c].tolist())
            self.conn_elems.append(self.elems[c])
        self.connector_xyz = np.array(conn_xyz,"d")
        self.conn_dist = np.sqrt(np.sum(self.connector_xyz*self.connector_xyz,axis=1))
        return

    def hide_dummy_atoms(self):
        self.sbu = copy.deepcopy(self)
        self.natoms = self.natoms - len(self.dummies)
        self.xyz = self.xyz[0:self.natoms,:]
        self.conn = self.conn[0:self.natoms]
        self.elems = self.elems[0:self.natoms]


    def rotate_on_z(self):
        """ especially if this is a linker (2 connectors) we want it to lie on the z-axis
        do this AFTER center but BEFORE extract_connector_xyz
        we always use the first connector (could also be a regular SBU!) to be on the z-axis """
        c1_xyz = self.xyz[self.connectors[0]]
        z_axis = np.array([0.0,0.0,1.0],"d")
        theta = vector.angle(z_axis,c1_xyz) # angle to rotate
        if (theta > 1.0e-10) and (theta < (np.pi-1.0e-10)): 
            axis  = vector.normalize(vector.cross_prod(z_axis,c1_xyz)) # axis around which we rotate
            self.xyz = vector.rotate(self.xyz, axis, -theta)
        return
    
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
        
        
        
        
