# -*- coding: utf-8 -*-

#from molsys import *
import string
import copy
import numpy as np

class bb(molsys):
    
    def __init__(self, name, specific_conn=None, linker=False, zflip=False, nrot=2, label = None):
        molsys.__init__(self)
        self.specific_conn = specific_conn
        self.linker = linker
        self.zflip  = zflip
        self.nrot   = nrot
        if not linker:
            if self.zflip: print "Warning: zflip only supported for linkers"
            if self.nrot>1: print "Warning: rotations only supported for linkers"
        self.connectors = []
        self.dummies=[]
        self.dummy_neighbors =None
        self.version = self.get_version(name)
        if self.version == 2: 
            self.read_v2txyz(name)
        else:
            self.read_txyz(name)
        self.center()
        if linker: self.rotate_on_z()
        self.extract_connector_xyz()
        self.hide_dummy_atoms()
        self.name = name
        self.label = label
        return
    
    def get_version(self,name):
        f=open(name,'r')
        line = f.readline().split()
        if line[1] == 'v2.0':
            return 2
        return 1
        
    def read_txyz(self, name):
        f = open(name, "r")
        line = st.split(f.readline())
        self.natoms = st.atoi(line[0])
        self.periodic=None
        self.center_point = line[1]
        con_info = line[2:]
        self.connectors = []
        self.connectors_type = []
        contype_count = 0
        for c in con_info:
            if c == "/":
                contype_count += 1
            else:
                self.connectors.append(st.atoi(c)-1)
                self.connectors_type.append(contype_count)
        if contype_count == 0: self.connectors_type = None
        xyz = []
        self.conn  =[]
        self.elems =[]
        self.atypes=[]
        for i in xrange(self.natoms):
            line = st.split(f.readline())
            self.elems.append(st.lower(line[1]))
            xyz.append(map(st.atof, line[2:5]))
            self.atypes.append(line[5])
            self.conn.append((np.array(map(st.atoi, line[6:]),"i")-1).tolist())
        self.xyz = np.array(xyz, "d")
        if self.center_point == "special":
            line = st.split(f.readline())
            self.special_center_point = np.array(map(st.atof, line[0:3]),"d")
        try:
            line = st.split(f.readline())
            if line != [] and line[0][:5] == 'angle':
                self.angleterm = line
        except:
            pass
        f.close()
        # sanity test on connector types
        if self.specific_conn:
            if self.connectors_type == None:
                raise ValueError, "SBU file does not define different connector types"
            if (contype_count+1) != len(self.specific_conn):
                raise ValueError, "Wrong number of connector types in SBU file"
        return
    
    def read_v2txyz(self, name):
        f = open(name, "r")
        line = f.readline().split()
        self.natoms = int(line[0])
        self.periodic=None
        self.center_point = line[2]
        con_info = line[3:]
        self.connectors = []
        self.connectors_type = []
        contype_count = 0
        if contype_count == 0: self.connectors_type = None
        xyz = []
        self.conn  =[]
        self.elems =[]
        self.atypes=[]
        for i in xrange(self.natoms):
            line = f.readline().split()
            self.elems.append(string.lower(line[1]))
            xyz.append(map(float, line[2:5]))
            self.atypes.append(line[5])
            self.conn.append((np.array(map(int, line[6:]),"i")-1).tolist())
        self.xyz = np.array(xyz, "d")
        if self.center_point == "special":
            line = f.readline().split()
            self.special_center_point = np.array(map(float, line[0:3]),"d")
        try:
            line = f.readline().split()
            if line != [] and line[0][:5] == 'angle':
                self.angleterm = line
        except:
            pass
        f.close()
        # sanity test on connector types
        if self.specific_conn:
            if self.connectors_type == None:
                raise ValueError, "SBU file does not define different connector types"
            if (contype_count+1) != len(self.specific_conn):
                raise ValueError, "Wrong number of connector types in SBU file"
            

        self.dummies = []
        self.dummy_neighbors=[]
        #print 'con_info', con_info
        for c in con_info:
            ss = c.split('*') # ss[0] is the dummy neighbors, ss[1] is the connector atom
            if len(ss) != 2: raise IOError('This is not a proper BB file, convert with script before!')
            stt = ss[0].split(',')
            #print 'c', c,'st',stt, ss
            self.connectors.append(int(ss[1])-1)
            #self.connectors_type.append(contype_count)
            # now we check whether the connector is a dummy or not: identify by 'X'(elems) or '666'(atypes)
            if string.lower(self.elems[int(ss[1])-1]) == 'x':
                #self.connectors.append(self.natoms)
                self.dummies.append(int(ss[1])-1) # simplest case only with two atoms being the connecting atoms
                #self.natoms += 1
            self.dummy_neighbors.append((np.array(map(int,stt)) -1).tolist())
            #print self.dummy_neighbors
            #center = (self.xyz[self.dummy_neighbors[-1][0]] + self.xyz[self.dummy_neighbors[-1][1]])/2.0
            #self.xyz = np.vstack([self.xyz,center])
            #self.conn.append([])
            #self.elems.append(str(len(self.dummies)))
            #self.atypes.append(str(len(self.dummies)))
                #print 'x1', self.xyz[self.dummy_neighbors[-1][0]]
                #print 'x2', self.xyz[self.dummy_neighbors[-1][1]]
                #print center
                
            #ctype = st[-1]
        #print self.dummies
        #print self.connectors
        #print self.dummy_neighbors
        #print self.natoms
        #print self.conn
        #print self.xyz.shape
        #print self.xyz
            
        return    

    def read_xyz(self, fname):
        print "Reading from xyz not supported for SBUs!"
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
    
    #def hide_dummy_atoms(self):
        #self.sbu = copy.depcopy(self)
        #self.sbu.natoms = self.natoms - len(self.dummies)
        #self.sbu.xyz = self.xyz[0:self.sbu.natoms,:]
        #self.sbu.conn = self.conn[0:self.sbu.natoms]
        #self.sbu.elems = self.elems[0:self.sbu.natoms]

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
        for i,d  in enumerate(self.dummy_neighbors):
            ci = np.sum(self.xyz[d],axis=0)/float(len(d))
            centers.append(ci)
        if shift_to_original==True:
            self.centers = centers+self.center_xyz
        else:
            self.centers = centers
        return centers
    
    def write_bb(self,fname=None):
        if not fname: fname = self.name.split('.')[0]+'.bb'
        
        if self.sbu.dummy_neighbors == None: self.sbu.dummy_neighbors = [[i] for i in self.connectors]
        
        f = open(fname, 'w')
        f.write('# type bb\n')
        connstrings = ''
        for i,d in enumerate(self.sbu.dummy_neighbors):   #we're using the backup here, beware!
            for j in d:
                connstrings = connstrings + str(j+1) + ','
            connstrings = connstrings[0:-1] + '*' + str(self.sbu.connectors[i]+1)+' ' # remove trailing comma and and star
        f.write('# bbcenter '+self.center_point+'\n')
        lastline = None
        if self.center_point == 'special':
            lastline = '%10.7f %10.7f %10.7f\n' & (self.special_center_point[0],self.special_center_point[1],self.special_center_point[2])
        f.write('# bbconn '+connstrings+'\n')
        f.write("%5d\n" % (self.sbu.natoms))
        for i in xrange(self.sbu.natoms):
            line = ("%3d %-3s" + 3*"%12.6f" + " %20s %15s %4d") % \
               tuple([i+1]+[self.sbu.elems[i]]+ self.sbu.xyz[i].tolist() + [self.sbu.atypes[i]]+['0']+[0])
            conn = (np.array(self.sbu.conn[i])+1).tolist()
            if len(conn) != 0:
                line += (len(conn)*"%7d") % tuple(conn)
            f.write("%s \n" % line)
        if lastline:
            f.write(lastline)
        f.close()
        return
    
    #def connectors_to_dummy_neighbors(self):
        
        
        
        
