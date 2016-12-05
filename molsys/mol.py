# -*- coding: utf-8 -*-
import string as st
import numpy as np
import types
import copy
import string

from util import unit_cell
from util import elems as elements
from util import rotations
from util import images
from io import formats

import addon

# set up logging using a logger
# note that this is module level because there is one logger for molsys
# DEBUG/LOG goes to logfile, whereas WARNIGNS/ERRORS go to stdout
#
# NOTE: this needs to be done once only here for the root logger molsys
# any other module can use either this logger or a child logger
# no need to redo this config in the other modules!
import logging
logger    = logging.getLogger("molsys")
logger.setLevel(logging.DEBUG)
fhandler  = logging.FileHandler("molsys.log")
fhandler.setLevel(logging.DEBUG)
shandler  = logging.StreamHandler()
shandler.setLevel(logging.WARNING)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m-%d %H:%M')
fhandler.setFormatter(formatter)
shandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.addHandler(shandler)




np.set_printoptions(threshold=20000)

SMALL_DIST = 1.0e-3

class mol:

    def __init__(self):
        self.natoms=0
        self.cell=None
        self.cellparams=None
        self.images_cellvec=None
        self.xyz=None
        self.elems=[]
        self.atypes=[]
        self.conn=[]
        self.fragtypes=[]
        self.fragnumbers=[]
        self.periodic= None
        self.is_bb=False
        return

    #####  I/O stuff ############################

    def read(self,fname,ftype='mfpx',**kwargs):
        ''' generic reader for the mol class
        :Parameters:
            - fname        : the filename to be red
            - ftype="mfpx" : the parser type that is used to read the file
            - **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        formats.read[ftype](self,fname,**kwargs)
        return

    def write(self,fname,ftype='mfpx',**kwargs):
        ''' generic writer for the mol class
        :Parameters:
            - fname        : the filename to be written
            - ftype="mfpx" : the parser type that is used to writen the file
            - **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        formats.write[ftype](self,fname,**kwargs)
        return


    ##### addons ##################################

    def addon(self, addmod):
        """
        add an addon module to this object

        :Parameters:

            - addmol: string name of the addon module
        """
        if addmod == "graph":

            if addon.graph != None:
                # ok, it is present and imported ...
                self.graph = addon.graph(self)
            else:
                logger.error("graph_tool is not installed! This addon can not be used")
        elif addmod  == "fragments":
            self.fragments = addon.fragments(self)
        elif addmod  == "bb":
            self.bb = addon.bb(self)
        else:
            logger.error("the addon %s is unknown")
        return

    ##### connectivity ########################

    def detect_conn(self, tresh = 0.1,remove_duplicates = False):
        """
        detects the connectivity of the system, based on covalent radii.

        :Parameters:
            - tresh  (float): additive treshhold 
            - remove_duplicates  (bool): flag for the detection of duplicates
        """
        xyz = self.xyz
        elems = self.elems
        if type(self.cell) != type(None):
            cell_abc = self.cellparams[:3]
            cell_angles = self.cellparams[3:]
            if cell_angles[0] != 90.0 or cell_angles[1] != 90.0 or cell_angles[2] != 90.0:
                inv_cell = np.linalg.inv(self.cell)
        natoms = self.natoms
        conn = []
        duplicates = []
        for i in xrange(natoms):
            a = xyz - xyz[i]
            if type(self.cell) != type(None):
                if cell_angles[0] == 90.0 and cell_angles[1] == 90.0 and cell_angles[2] == 90.0:
                    a -= cell_abc * np.around(a/cell_abc)
                else:
                    frac = np.dot(a, inv_cell)
                    frac -= np.around(frac)
                    a = np.dot(frac, self.cell)
            dist = ((a**2).sum(axis=1))**0.5 # distances from i to all other atoms
            conn_local = []
            if remove_duplicates == True:
                for j in xrange(i,natoms):
                    if i != j and dist[j] < tresh:
                        logger.warning("atom %i is duplicate of atom %i" % (j,i))
                        duplicates.append(j)
            else:
                for j in xrange(natoms):
                    if i != j and dist[j] <= elements.get_covdistance([elems[i],elems[j]])+tresh:
                        conn_local.append(j)
            if remove_duplicates == False: conn.append(conn_local)
        if remove_duplicates:
            if len(duplicates)>0:
                logger.warning("Found %d duplicates" % len(duplicates))
                self.natoms -= len(duplicates)
                self.set_xyz(np.delete(xyz, duplicates,0))
                self.set_elements(np.delete(elems, duplicates))
                self.set_atypes(np.delete(self.atypes,duplicates))
                self.set_fragtypes(np.delete(self.fragtypes,duplicates))
                self.set_fragnumbers(np.delete(self.fragnumbers,duplicates))
            self.detect_conn(tresh = tresh)
        else:
            self.set_conn(conn)
        return

    def report_conn(self):
        ''' Print infomration on current connectivity, coordination number
            and the respective atomic distances '''
        print "REPORTING CONNECTIVITY"
        for i in xrange(self.natoms):
            conn = self.conn[i]
            print "atom %3d   %2s coordination number: %3d" % (i, self.elems[i], len(conn))
            for j in xrange(len(conn)):
                d = self.get_neighb_dist(i,j)
                print "   -> %3d %2s : dist %10.5f " % (conn[j], self.elems[conn[j]], d)
        return

    ###  periodic systems .. cell manipulation ############

    def make_supercell(self,supercell):
        ''' Extends the periodic system in all directions by the factors given in the
            supercell upon preserving the connectivity of the initial system
            :Parameters:
                - supercell: List of integers, e.g. [3,2,1] extends the cell three times in x and two times in y'''
        img = [np.array(i) for i in images.tolist()]
        ntot = np.prod(supercell)
        nat = copy.deepcopy(self.natoms)
        nx,ny,nz = supercell[0],supercell[1],supercell[2]
        #pconn = [copy.deepcopy(self.pconn) for i in range(ntot)]
        conn =  [copy.deepcopy(self.conn) for i in range(ntot)]
        xyz =   [copy.deepcopy(self.xyz) for i in range(ntot)]
        elems = copy.deepcopy(self.elems)
        left,right,front,back,bot,top =  [],[],[],[],[],[]
        neighs = [[] for i in range(6)]
        iii = []
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    ixyz = ix+nx*iy+nx*ny*iz
                    iii.append(ixyz)
                    if ix == 0   : left.append(ixyz)
                    if ix == nx-1: right.append(ixyz)
                    if iy == 0   : bot.append(ixyz)
                    if iy == ny-1: top.append(ixyz)
                    if iz == 0   : front.append(ixyz)
                    if iz == nz-1: back.append(ixyz)
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    ixyz = ix+nx*iy+nx*ny*iz
                    dispvect = np.sum(self.cell*np.array([ix,iy,iz])[:,np.newaxis],axis=0)
                    xyz[ixyz] += dispvect
                    i = copy.copy(ixyz)
                    for cc in range(len(conn[i])):
                        for c in range(len(conn[i][cc])):
                            pc = self.get_distvec(cc,conn[i][cc][c])[2]
                            if len(pc) != 1:
                                print self.get_distvec(cc,conn[i][cc][c])
                                print c,conn[i][cc][c]
                                raise ValueError('an Atom is connected to the same atom twice in different cells! \n requires pconn!! use topo molsys instead!')
                            pc = pc[0]
                            if pc == 13:
                                conn[i][cc][c] = int( conn[i][cc][c] + ixyz*nat )
                            else:
                                px,py,pz     = img[pc][0],img[pc][1],img[pc][2]
                                iix,iiy,iiz  = (ix+px)%nx, (iy+py)%ny, (iz+pz)%nz
                                iixyz= iix+nx*iiy+nx*ny*iiz
                                conn[i][cc][c] = int( conn[i][cc][c] + iixyz*nat )

        self.conn, self.xyz = [],[]
        for cc in conn:
            for c in cc:
                self.conn.append(c)
        self.natoms = nat*ntot
        self.xyz = np.array(xyz).reshape(nat*ntot,3)
        self.cellparams[0:3] *= np.array(supercell)
        self.cell *= np.array(supercell)[:,np.newaxis]
        self.elems *= ntot
        self.atypes*=ntot
        self.fragtypes*=ntot
        self.fragnumbers*=ntot
        self.images_cellvec = np.dot(images, self.cell)
        #print xyz
        return xyz,conn

    def wrap_in_box(self, thresh=SMALL_DIST):
        ''' In case atoms are outside the box defined by the cell,
            this routine finds and shifts them into the box'''
        if not self.periodic: return
        # puts all atoms into the box again
        frac_xyz = self.get_frac_xyz()
        # now add 1 where the frac coord is negative and subtract where it is larger then 1
        frac_xyz += np.where(np.less(frac_xyz, 0.0), 1.0, 0.0)
        frac_xyz -= np.where(np.greater_equal(frac_xyz, 1.0+thresh*0.1), 1.0, 0.0)
        # convert back
        self.set_xyz_from_frac(frac_xyz)
        return

    def get_frac_xyz(self):
        ''' Returns the fractional atomic coordinates'''
        if not self.periodic: return None
        cell_inv = np.linalg.inv(self.cell)
        return np.dot(self.xyz, cell_inv)

    def get_frac_from_real(self,real_xyz):
        ''' same as get_frac_xyz, but uses input xyz coordinates
        :Parameters:
            - real_xyz: the xyz coordinates for which the fractional coordinates are retrieved'''
        if not self.periodic: return None
        cell_inv = np.linalg.inv(self.cell)
        return np.dot(real_xyz, cell_inv)

    def get_real_from_frac(self,frac_xyz):
        ''' returns real coordinates from an array of fractional coordinates using the current cell info '''
        return np.dot(np.array(frac_xyz),self.cell)

    def set_xyz_from_frac(self, frac_xyz):
        ''' Sets atomic coordinates based on input fractional coordinates
        Parameter:
            - '''
        if not self.periodic: return
        self.xyz = np.dot(frac_xyz,self.cell)

    def get_image(self,xyz, img):
        ''' returns the xyz coordinates of a set of coordinates in a specific cell
        :Parameters:
            - xyz   : xyz coordinates for which the image coordinates are to be retrieved
            - img   : descriptor of the image, either an "images" integer (see molsys.util.images)
                      or the unit direction vector, e.g. [1,-1,0]'''
        xyz = np.array(xyz)
        try:
            l = len(img)
            dispvec = np.sum(self.cell*np.array(img)[:,np.newaxis],axis=0)
        except TypeError:
            dispvec = np.sum(self.cell*np.array(images[img])[:,np.newaxis],axis=0)
        return xyz + dispvec

    ### rewrite on set_cell ???
    def scale_cell(self, scale):
        ''' scales the cell by a given fraction (0.1 ^= 10%)
        :Parameters:
            - scale: either single float or list (3,) of floats for x,y,z'''
        if not type(scale) == types.ListType:
            scale = 3*[scale]
        self.cellparams *= np.array(scale+[1,1,1])
        frac_xyz = self.get_frac_xyz()
        self.cell *= np.array(scale)[:,np.newaxis]
        self.images_cellvec = np.dot(images, self.cell)
        self.set_xyz_from_frac(frac_xyz)
        return

    ###  system manipulations ##########################################

    def copy(self):
        ''' returns a copy of the whole mol object'''
        return copy.deepcopy(self)

    def add_mol(self, other, translate=None,rotate=None, scale=None, roteuler=None):
        ''' adds a  nonperiodic mol object to the current one ... self can be both
            :Parameters:
                - other        (mol): an instance of the to-be-inserted mol instance
                - translate=None    : (3,) numpy array as shift vector for the other mol
                - rotate=None       : (3,) rotation triple to apply to the other mol object before insertion
                - scale=None (float): scaling factor for other mol object coodinates
                - roteuler=None     : (3,) euler angles to apply a rotation prior to insertion'''
        if other.periodic:
            if not (self.cell==other.cell).all():
                print "can not add periodic systems with unequal cells!!"
                return
        other_xyz = other.xyz.copy()
        # NOTE: it is important ot keep the order of operations
        #       1) scale
        #       2) rotate by euler angles
        #       3) rotate by orientation triple
        #       4) translate
        if scale    !=None:
            other_xyz *= np.array(scale)
        if roteuler != None:
            other_xyz = rotations.rotate_by_euler(other_xyz, roteuler)
        if rotate   !=None:
            other_xyz = rotations.rotate_by_triple(other_xyz, rotate)
        if translate!=None:
            other_xyz += translate
        if self.natoms==0:
            self.xyz = other_xyz
        else:
            self.xyz = np.array(self.xyz.tolist()+other_xyz.tolist(),"d")
        self.elems += other.elems
        self.atypes+= other.atypes
        for c in other.conn:
            cn = (np.array(c)+self.natoms).tolist()
            self.conn.append(cn)
        self.natoms += other.natoms
        self.fragtypes += other.fragtypes
        self.fragnumbers += other.fragnumbers
        return

    def add_bond(self,a1,a2):
        self.conn[a1].append(a2)
        self.conn[a2].append(a1)
    ###  molecular manipulations #######################################

    def delete_atom(self,bad):
        ''' deletes an atom and its connections and fixes broken indices of all other atoms '''
        new_xyz = []
        new_elems = []
        new_atypes = []
        new_conn = []
        for i in xrange(self.natoms):
            if i != bad:
                new_xyz.append(self.xyz[i].tolist())
                new_elems.append(self.elems[i])
                new_atypes.append(self.atypes[i])
                new_conn.append(self.conn[i])
                for j in xrange(len(new_conn[-1])):
                    if new_conn[-1].count(bad) != 0:
                        new_conn[-1].pop(new_conn[-1].index(bad))
        self.xyz = np.array(new_xyz, "d")
        self.elems = new_elems
        self.natoms = len(self.elems)
        self.atypes = new_atypes
        for i in range(len(new_conn)):
            #try:
                #len(new_conn[i])
            #except:
                #new_conn[i] = [new_conn[i]]
            for j in range(len(new_conn[i])):
                if new_conn[i][j] >= bad:
                    new_conn[i][j]=new_conn[i][j]-1
        self.conn = new_conn
        return

    def translate(self, vec):
        self.xyz += vec
        return

    def translate_frac(self, vec):
        if not self.periodic: return
        self.xyz += np.sum(self.cell*vec, axis=0)
        return

    def rotate_euler(self, euler):
        self.xyz = rotations.rotate_by_euler(self.xyz, euler)
        return

    def rotate_triple(self, triple):
        self.xyz = rotations.rotate_by_triple(self.xyz, triple)
        return

    def center_com(self):
        ''' centers the molsys at the center of mass '''
        if self.periodic: return
        amass = []
        for e in self.elems: amass.append(elements.mass[e])
        amass = np.array(amass)
        center = np.sum((self.xyz*amass[:,np.newaxis]),axis=0)/np.sum(amass)
        self.translate(-center)
        return

    ##### distance measurements #####################

    def get_distvec(self, i, j):
        """ vector from i to j
        This is a tricky bit, because it is needed also for distance detection in the blueprint
        where there can be small cell params wrt to the vertex distances.
        In other words: i can be bonded to j multiple times (each in a different image)
        and i and j could be the same!!
        :Parameters':
            - i,j  : the indices of the atoms for which the distance is to be calculated"""
        ri = self.xyz[i]
        rj = self.xyz[j]
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - ri
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            d_sort = np.argsort(all_d)
            if i == j:
                # if this was requested for i==j then we have to eliminate the shortest
                # distance which NOTE unfinished!!!!!!!!
                pass
            closest = d_sort[0]
            closest=[closest]  # THIS IS A BIT OF A HACK BUT WE MAKE IT ALWAYS A LIST ....
            if (abs(all_d[closest[0]]-all_d[d_sort[1]]) < SMALL_DIST):
                # oops ... there is more then one image atom in the same distance
                #  this means the distance is larger then half the cell width
                # in this case we have to return a list of distances
                for k in d_sort[1:]:
                    if (abs(all_d[d_sort[0]]-all_d[k]) < SMALL_DIST):
                        closest.append(k)
            d = all_d[closest[0]]
            r = all_r[closest[0]]
        else:
            if i == j: return
            r = rj-ri
            d = np.sqrt(np.sum(r*r))
            closest=[0]
        return d, r, closest

    def get_neighb_coords(self, i, ci):
        """ returns coordinates of atom bonded to i which is ci'th in bond list
        :Parameters:
            - i  :  index of the base atom
            - ci :  index of the conn entry of the ith atom"""
        j = self.conn[i][ci]
        rj = self.xyz[j].copy()
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - self.xyz[i]
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            closest = np.argsort(all_d)[0]
            return all_rj[closest]
        return rj

    def get_neighb_dist(self, i, ci):
        """ returns coordinates of atom bonded to i which is ci'th in bond list
        :Parameters:
            - i  :  index of the base atom
            - ci :  index of the conn entry of the ith atom"""
        ri = self.xyz[i]
        j = self.conn[i][ci]
        rj = self.xyz[j].copy()
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - self.xyz[i]
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            closest = np.argsort(all_d)[0]
            return all_rj[closest]
        dr = ri-rj
        d = np.sqrt(np.sum(dr*dr))
        return d

    def get_com(self):
        ''' calculates and returns the center of mass'''
        amass = []
        for e in self.elems: amass.append(elements.mass[e])
        amass = np.array(amass,dtype='float64')
        #print amass, amass[:,np.newaxis].shape,self.xyz.shape
        center = np.sum((amass[:,np.newaxis]*self.xyz),axis=0)/np.sum(amass)
        #center = np.sum((amass[:,np.newaxis]*self.xyz),axis=0)/np.sum(amass)
        return center

    def get_comdist(self,com,i):
        ''' Calculate the distances of an atom i from a given point (e.g. the center of mass)
        :Parameters:
            - com : center of mass
            - i   : index of the atom for which to calculate the distances to the com'''
        ri = self.xyz[i]
        rj = com
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - ri
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            d_sort = np.argsort(all_d)
            closest = d_sort[0]
            closest=[closest]  # THIS IS A BIT OF A HACK BUT WE MAKE IT ALWAYS A LIST ....
            if (abs(all_d[closest[0]]-all_d[d_sort[1]]) < SMALL_DIST):
                # oops ... there is more then one image atom in the same distance
                #  this means the distance is larger then half the cell width
                # in this case we have to return a list of distances
                for k in d_sort[1:]:
                    if (abs(all_d[d_sort[0]]-all_d[k]) < SMALL_DIST):
                        closest.append(k)
            d = all_d[closest[0]]
            r = all_r[closest[0]]
        else:
            r = rj-ri
            d = np.sqrt(np.sum(r*r))
            closest=[0]
        return d, r, closest


    ### get and set methods ###

    def get_natoms(self):
        ''' returns the number of Atoms '''
        return self.natoms

    def get_xyz(self):
        ''' returns the xyz Coordinates '''
        return self.xyz

    def set_xyz(self,xyz):
        ''' set the real xyz coordinates
        :Parameters:
            - xyz: coordinates to be set'''
        assert np.shape(xyz) == (self.natoms,3)
        self.xyz = xyz

    def get_elems(self):
        ''' return the list of element symbols '''
        return self.elems

    def get_elemlist(self):
        ''' Returns a list of unique elements '''
        el = []
        for e in self.elems:
            if not el.count(e): el.append(e)
        return el

    def set_elements(self,elems):
        ''' set the elements
        :Parameters:
            - elems: list of elements to be set'''
        assert len(elems) == self.natoms
        self.elems = elems

    def get_atypes(self):
        ''' return the list of atom types '''
        return self.atypes

    def get_atypelist(self):
        ''' Returns a list of unique atom types '''
        if not self.atypes: return None
        at = []
        for a in self.atypes:
            if not at.count(a): at.append(a)
        return at

    def set_atypes(self,atypes):
        ''' set the atomtypes
        :Parameters:
            - atypes: list of elements to be set'''
        assert len(atypes) == self.natoms
        self.atypes = atypes

    def get_cell(self):
        ''' return unit cell information (a, b, c, alpha, beta, gamma) '''
        return self.cell

    def get_cellparams(self):
        ''' return unit cell information (cell vectors) '''
        return self.cellparams

    def set_cell(self,cell,cell_only = True):
        ''' set unit cell using cell vectors and assign cellparams
        :Parameters:
            - cell: cell vectors (3,3)
            - cell_only (bool)  : if false, also the coordinates are changed
                                  in respect to new cell
            
        '''
        assert np.shape(cell) == (3,3)
        if cell_only == False: frac_xyz = self.get_frac_xyz()
        self.periodic = True
        self.cell = cell
        self.cellparams = unit_cell.abc_from_vectors(self.cell)
        self.images_cellvec = np.dot(images, self.cell)
        if cell_only == False: self.set_xyz_from_frac(frac_xyz)

    def set_cellparams(self,cellparams, cell_only = True):
        ''' set unit cell using cell parameters and assign cell vectors
        :Parameters:
            - cell: cell vectors (3,3)
            - cell_only (bool)  : if false, also the coordinates are changed
                                  in respect to new cell
        '''
        assert len(cellparams) == 6
        if cell_only == False: frac_xyz = self.get_frac_xyz()
        self.periodic = True
        self.cellparams = cellparams
        self.cell = unit_cell.vectors_from_abc(self.cellparams)
        self.images_cellvec = np.dot(images, self.cell)
        if cell_only == False: self.set_xyz_from_frac(frac_xyz)

    def get_fragtypes(self):
        ''' return all fragment types '''
        return self.fragtypes
        
    def get_fragtypes_list(self,count=False):
        ''' return a list of unique fragment types '''
        lset = list(set(self.fragtypes))
        if not count: return lset
        counts = []
        for i,ls in enumerate(lset):
            counts.append(self.fragtypes.count(ls))
        return [lset,counts]

    def set_fragtypes(self,fragtypes):
        ''' set fragment types
        :Parameters:
            - fragtypes: the fragtypes to be set (list of strings)'''
        assert len(fragtypes) == self.natoms
        self.fragtypes = fragtypes

    def get_fragnumbers(self):
        ''' return all fragment numbers, denotes which atom belongs to which fragment '''
        return self.fragnumbers

    def set_fragnumbers(self,fragnumbers):
        ''' set fragment numbers, denotes which atom belongs to which fragment
        :Parameters:
            - fragnumbers: the fragment numbers to be set (list of integers)'''
        assert len(fragnumbers) == self.natoms
        self.fragnumbers = fragnumbers

    def get_conn(self):
        ''' returns the connectivity of the system '''
        return self.conn

    def set_conn(self,conn):
        ''' updates the connectivity of the system
        :Parameters:
            - conn    : List of lists describing the connectivity'''
        self.conn = conn

    def set_nofrags(self):
        ''' in case there are no fragment types and numbers, setup the data structure which is needed in some functions '''
        self.set_fragtypes(['0']*self.natoms)
        self.set_fragnumbers([0]*self.natoms)



