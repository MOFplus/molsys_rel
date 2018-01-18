# -*- coding: utf-8 -*-
### overload print in parallel case (needs to be the first line) [RS] ###
from __future__ import print_function

import string as st
import numpy as np
from scipy.optimize import linear_sum_assignment as hungarian
import types
import copy
import string
import os
import sys
import subprocess
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    mpi_rank = MPI.COMM_WORLD.Get_rank()
    mpi_size = MPI.COMM_WORLD.Get_size()
except ImportError as e:
    mpi_comm = None
    mpi_size = 1
    mpi_rank = 0
    mpi_err = e

from .util import unit_cell
from .util import elems as elements
from .util import rotations
from .util import images
from .fileIO import formats

from . import addon

# set up logging using a logger
# note that this is module level because there is one logger for molsys
# DEBUG/LOG goes to logfile, whereas WARNIGNS/ERRORS go to stdout
#
# NOTE: this needs to be done once only here for the root logger molsys
# any other module can use either this logger or a child logger
# no need to redo this config in the other modules!
# NOTE2: in a parallel run all DEBUG is written by all nodes whereas only the 
#        master node writes INFO to stdout
# TBI: colored logging https://stackoverflow.com/a/384125
import logging
logger    = logging.getLogger("molsys")
logger.setLevel(logging.DEBUG)
if mpi_size > 1:
    logger_file_name = "molsys.%d.log" % mpi_rank
else:
    logger_file_name = "molsys.log"
fhandler  = logging.FileHandler(logger_file_name)
#fhandler.setLevel(logging.DEBUG)
fhandler.setLevel(logging.WARNING)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m-%d %H:%M')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
if mpi_rank == 0:
    shandler  = logging.StreamHandler()
    #shandler.setLevel(logging.INFO)
    shandler.setLevel(logging.WARNING)
    shandler.setFormatter(formatter)
    logger.addHandler(shandler)
    
if mpi_comm == None:
    logger.error("MPI NOT IMPORTED DUE TO ImportError")
    logger.error(mpi_err)

# overload print function in parallel case
try:
    import __builtin__
except ImportError:
    import builtins as __builtin__
def print(*args, **kwargs):
    if mpi_rank == 0:
        return __builtin__.print(*args, **kwargs)
    else:
        return



np.set_printoptions(threshold=20000,precision=5)

SMALL_DIST = 1.0e-3

class mpiobject(object):
    """Basic class to handle parallel attributes
    
    :Parameters:
    - mpi_comm(obj): Intracomm object for parallelization
    - out(str): output file name. If None: stdout
    """
    def __init__(self, mpi_comm = None, out = None):
        if mpi_comm is None:
            try:
                self.mpi_comm = MPI.COMM_WORLD
            except NameError:
                self.mpi_comm = None
        else:
            self.mpi_comm = mpi_comm
        try:
            self.mpi_rank = self.mpi_comm.Get_rank()
            self.mpi_size = self.mpi_comm.Get_size()
        except AttributeError:
            self.mpi_rank = 0
            self.mpi_size = 1
        if out is None:
            self.out = sys.stdout
        elif type(out) == file:
            assert out.mode == 'w'
            self.out = out
        else:
            self.out = open(out, "w")

    def pprint(self, *args, **kwargs):
        """Parallel print function"""
        if self.is_master:
            __builtin__.print(*args, file=self.out, **kwargs)
            self.out.flush()

    @property
    def is_master(self):
        """
        The mpi process with global rank 0 is always the master rank.
        This methods returns True if the current process has rank 0, else
        False
        """
        return self.mpi_rank == 0


    def __getstate__(self):
        """Get state for pickle and pickle-based method (e.g. copy.deepcopy)
        Meant for python3 forward-compatibility.
        Files (and standard output/error as well) are _io.TextIOWrapper
        in python3, and thus they are not pickle-able. A dedicated method for
        files is needed.
        An example: https://docs.python.org/2.0/lib/pickle-example.html

        N.B.: experimental for "out" != sys.stdout
        """
        newone = type(self)()
        newdict = newone.__dict__
        newdict.update(self.__dict__)
        newdict["out.name"] = newdict["out"].name
        newdict["out.mode"] = newdict["out"].mode
        newdict["out.encoding"] = newdict["out"].encoding
        del newdict["out"]
        return newdict

    def __setstate__(self, stored_dict):
        """Set state for pickle and pickle-based method (e.g. copy.deepcopy)
        For python3 forward-compatibility
        Whatever comes out of getstate, goes int setstate.
        https://stackoverflow.com/a/41754104

        N.B.: experimental for "out" != sys.stdout
        """
        if stored_dict["out.name"] == '<stdout>':
            stored_dict["out"] = sys.stdout
        self.__dict__ = stored_dict

class mol(mpiobject):
    """mol class, the basis for any atomistic (or atomistic-like,
    e.g. topo) representation."""

    def __init__(self, mpi_comm = None, out = None):
        super(mol,self).__init__(mpi_comm, out)
        self.natoms=0
        self.cell=None
        self.cellparams=None
        self.images_cellvec=None
        self.bcond = 0
        self.xyz=None
        self.elems=[]
        self.atypes=[]
        self.conn=[]
        self.ctab=[]
        self.fragtypes=[]
        self.fragnumbers=[]
        self.nfrags = 0
        self.periodic=False
        self.is_bb=False
        self.weight=1
        self.loaded_addons =  []
        self.set_logger_level()
        return

    #####  I/O stuff ############################

    def set_logger_level(self,level='INFO'):
        if level=='INFO':
            logger.setLevel(logging.INFO)
        if level=='WARNING':
            logger.setLevel(logging.WARNING)
        if level=='ERROR':
            logger.setLevel(logging.ERROR)
        if level=='DEBUG':
            logger.setLevel(logging.DEBUG)
        return

    def read(self, fname, ftype=None, **kwargs):
        ''' generic reader for the mol class
        :Parameters:
            - fname        : the filename to be read
            - ftype="mfpx" : the parser type that is used to read the file
            - **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        if ftype is None:
            fsplit = fname.rsplit('.',1)[-1]
            if fsplit != fname: #there is an extension
                ftype = fsplit #ftype is inferred from extension
            else: #there is no extension
                ftype = 'mfpx' #default
        logger.info("reading file %s in %s format" % (fname, ftype))
        try:
            f = open(fname, "r")
        except IOError:
            fname += "."+ftype
            f = open(fname, "r")
        if ftype in formats.read:
            formats.read[ftype](self,f,**kwargs)
        else:
            logger.error("unsupported format: %s" % ftype)
            raise IOError("Unsupported format")
        return
    
    @classmethod
    def fromFile(cls, fname, ftype=None, **kwargs):
        ''' reader for the mol class, reading from a file
        :Parameters:
            - fname(str): path to the file (filename included)
            - ftype=None (or str): the parser type that is used to read the file
                if None: assigned by read as mfpx (default)
            - **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        m = cls()
        m.read(fname, ftype, **kwargs)
        return m
    
    @classmethod
    def fromString(cls, istring, ftype='mfpx', **kwargs):
        ''' generic reader for the mol class, reading from a string
        :Parameters:
            - string       : the string to be read
            - ftype="mfpx" : the parser type that is used to read the file
            - **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        m = cls()
        logger.info("reading string as %s" % str(ftype))
        f = StringIO(istring)
        if ftype in formats.read:
            formats.read[ftype](m,f,**kwargs)
        else:
            logger.error("unsupported format: %s" % ftype)
            raise IOError("Unsupported format")
        return m

    @classmethod
    def fromAbinit(cls, elems, xyz, cell, frac = False):
        m = cls()
        logger.info('reading basic data provided by any AbInitio programm')
        m.natoms = len(elems)
        m.set_elems(elems)
        m.set_atypes(elems)
        m.set_cell(cell)
        if frac:
            m.set_xyz_from_frac(xyz)
        else:
            m.set_xyz(xyz)
        m.set_nofrags()
        m.set_empty_conn()
        m.detect_conn()
        return m

    @classmethod
    def fromPymatgen(cls, structure):
        m = cls()
        logger.info('creating mol object from a pymatgen structure object')
        cell = structure.lattice.matrix
        fracs = [] 
        elems = []
        for j, site in enumerate(structure.sites):
            elems.append(site.specie.symbol.lower())
            fracs.append([site.frac_coords[0],site.frac_coords[1], site.frac_coords[2]])
        fracs = np.array(fracs)
        m.natoms=len(elems)
        m.set_elems(elems)
        m.set_atypes(elems)
        m.set_cell(cell)
        m.set_xyz_from_frac(fracs)
        m.set_nofrags()
        m.set_empty_conn()
        m.detect_conn()
        return m
        
    
    @classmethod
    def fromArray(cls, arr, **kwargs):
        ''' generic reader for the mol class, reading from a Nx3 array
        :Parameters:
            - arr         : the array to be read
            - **kwargs    : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        m = cls()
        logger.info("reading array")
        assert arr.shape[1] == 3, "Wrong array dimension (second must be 3): %s" % (a.shape,)
        formats.read['array'](m,arr,**kwargs)
        return m

    @classmethod
    def fromNestedList(cls, nestl, **kwargs):
        ''' generic reader for the mol class, reading from a Nx3 array
        :Parameters:
            - arr         : the array to be read
            - **kwargs    : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        logger.info("reading nested lists")
        for nl in nestl:
            assert len(nl) == 3, "Wrong nested list lenght (must be 3): %s" % (a.shape,)
        arr = np.array(nestl)
        return cls.fromArray(arr, **kwargs)

    def write(self, fname, ftype=None, **kwargs):
        ''' generic writer for the mol class
        :Parameters:
            - fname        : the filename to be written
            - ftype="mfpx" : the parser type that is used to writen the file
            - **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        if self.mpi_rank == 0:
            if ftype is None:
                fsplit = fname.rsplit('.',1)[-1]
                if fsplit != fname: #there is an extension
                    ftype = fsplit #ftype is inferred from extension
                else: #there is no extension
                    ftype = 'mfpx' #default
            logger.info("writing file "+str(fname)+' in '+str(ftype)+' format')
            if ftype in formats.read:
                formats.write[ftype](self,fname,**kwargs)
            else:
                logger.error("unsupported format: %s" % ftype)
                raise IOError("Unsupported format")
        return

    def view(self, program='moldenx', fmt='mfpx', **kwargs):
        ''' launch graphics visualisation tool, i.e. moldenx.
        Debugging purpose.'''
        if self.mpi_rank == 0:
            logger.info("invoking %s as visualisation tool" % (program,))
            pid = str(os.getpid())
            _tmpfname = "_tmpfname_%s.%s" % (pid, fmt)
            self.write(_tmpfname)
            try:
                ret = subprocess.call([program, _tmpfname])
            except KeyboardInterrupt:
                pass
            finally:
                os.remove(_tmpfname)
                logger.info("temporary file "+_tmpfname+" removed")
        return
    
    def molden(self, **kwargs):
        self.view(program='moldenx', fmt='mfpx', **kwargs)
    def pymol(self, **kwargs):
        self.view(program='pymol', fmt='mfpx', **kwargs)

    ##### addons ##################################

    def addon(self, addmod, **kwargs):
        """
        add an addon module to this object

        :Parameters:

            - addmol: string name of the addon module
        """
        if addmod in self.loaded_addons:
            return
        if addmod == "graph":
            if addon.graph != None:
                # ok, it is present and imported ...
                self.graph = addon.graph(self)
            else:
                logger.error("graph_tool is not installed! This addon can not be used")
                return
        elif addmod  == "fragments":
            self.fragments = addon.fragments(self)
        elif addmod  == "bb":
            self.bb = addon.bb(self)
        elif addmod  == "zmat":
            if addon.zmat != None:
                self.zmat = addon.zmat(self)
            else:
                logger.error("pandas/chemcoord is not available! This addon can not be used")
                return
        elif addmod  == "spg":
            if addon.spg != None:
                self.spg = addon.spg(self)
            else:
                logger.error("spglib is not available! This addon can not be used")
                return
        elif addmod == "ric":
            if addon.ric != None:
                self.ric = addon.ric(self)
            else:
                logger.error("ric is not available! This addon can not be used")
                return
        elif addmod == "ff":
            self.ff = addon.ff(self, **kwargs)
        elif addmod == "molecules":
            self.molecules = addon.molecules(self)
        else:
            logger.error("the addon %s is unknown" % addmod)
            return
        self.loaded_addons.append(addmod)
        return

    ##### connectivity ########################

    @staticmethod
    def check_conn(conn):
        """
        checks if connectivity is not broken

        :Parameters:
            - conn (list): list of lists holding the connectivity
        """
        natoms = len(conn)
        for i, c in enumerate(conn):
            for j in c:
                if i not in conn[j]: return False
        return True

    def detect_conn(self, tresh = 0.1,remove_duplicates = False):
        """
        detects the connectivity of the system, based on covalent radii.

        :Parameters:
            - tresh  (float): additive treshhold
            - remove_duplicates  (bool): flag for the detection of duplicates
        """

        logger.info("detecting connectivity by distances ... ")

        xyz = self.xyz
        elems = self.elems
        natoms = self.natoms
        conn = []
        duplicates = []
        for i in range(natoms):
            a = xyz - xyz[i]
            if self.periodic:
                if self.bcond <= 2:
                    cell_abc = self.cellparams[:3]
                    a -= cell_abc * np.around(a/cell_abc)
                elif self.bcond == 3:
                    frac = np.dot(a, self.inv_cell)
                    frac -= np.around(frac)
                    a = np.dot(frac, self.cell)
            dist = ((a**2).sum(axis=1))**0.5 # distances from i to all other atoms
            conn_local = []
            if remove_duplicates == True:
                for j in range(i,natoms):
                    if i != j and dist[j] < tresh:
                        logger.warning("atom %i is duplicate of atom %i" % (j,i))
                        duplicates.append(j)
            else:
                for j in range(natoms):
                    if i != j and dist[j] <= elements.get_covdistance([elems[i],elems[j]])+tresh:
                        conn_local.append(j)
            if remove_duplicates == False: conn.append(conn_local)
        if remove_duplicates:
            if len(duplicates)>0:
                logger.warning("Found %d duplicates" % len(duplicates))
                self.natoms -= len(duplicates)
                self.set_xyz(np.delete(xyz, duplicates,0))
                self.set_elems(np.delete(elems, duplicates))
                self.set_atypes(np.delete(self.atypes,duplicates))
                self.set_fragtypes(np.delete(self.fragtypes,duplicates))
                self.set_fragnumbers(np.delete(self.fragnumbers,duplicates))
            self.detect_conn(tresh = tresh)
        else:
            self.set_conn(conn)
        return

    def report_conn(self):
        ''' Print information on current connectivity, coordination number
            and the respective atomic distances '''

        logger.info("reporting connectivity ... ")
        for i in range(self.natoms):
            conn = self.conn[i]
            self.pprint("atom %3d   %2s coordination number: %3d" % (i, self.elems[i], len(conn)))
            for j in range(len(conn)):
                d = self.get_neighb_dist(i,j)
                self.pprint("   -> %3d %2s : dist %10.5f " % (conn[j], self.elems[conn[j]], d))
        return

    ###  periodic systems .. cell manipulation ############

    def make_supercell(self,supercell):
        ''' Extends the periodic system in all directions by the factors given in the
            supercell upon preserving the connectivity of the initial system
            :Parameters:
                - supercell: List of integers, e.g. [3,2,1] extends the cell three times in x and two times in y'''
        self.supercell = tuple(supercell)
        ntot = np.prod(self.supercell)
        conn =  [copy.deepcopy(self.conn) for i in range(ntot)]
        xyz =   [copy.deepcopy(self.xyz) for i in range(ntot)]
        if sum(self.supercell) == 3:
            logger.warning('Generating %i x %i x %i supercell? No need to do that!' % self.supercell)
            return xyz,conn
        logger.info('Generating %i x %i x %i supercell' % self.supercell)
        img = [np.array(i) for i in images.tolist()]
        nat = copy.deepcopy(self.natoms)
        nx, ny, nz = self.supercell
        #pconn = [copy.deepcopy(self.pconn) for i in range(ntot)]
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
                                print(self.get_distvec(cc,conn[i][cc][c]))
                                print(c,conn[i][cc][c])
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
        cell = self.cell * np.array(self.supercell)[:,np.newaxis]
        self.set_cell(cell)
        self.inv_cell = np.linalg.inv(self.cell)
        self.elems *= ntot
        self.atypes*=ntot
        self.fragtypes*=ntot
        mfn = max(self.fragnumbers)+1
        nfragnumbers = []
        for i in range(ntot):
            nfragnumbers += list(np.array(self.fragnumbers)+i*mfn)
        self.fragnumbers=nfragnumbers
        self.images_cellvec = np.dot(images, self.cell)
        return xyz,conn

    def wrap_in_box(self, thresh=SMALL_DIST):
        ''' In case atoms are outside the box defined by the cell,
            this routine finds and shifts them into the box'''
        if not self.periodic: return
        # puts all atoms into the box again
        frac_xyz = self.get_frac_xyz()
        # now add 1 where the frac coord is negative and subtract where it is larger then 1
        frac_xyz += np.where(np.less(frac_xyz, 0.0), 1.0, 0.0)
        frac_xyz -= np.where(np.greater_equal(frac_xyz, 1.0), 1.0, 0.0)
        # convert back
        self.set_xyz_from_frac(frac_xyz)
        if self.__class__.__name__=='topo':
            self.add_pconn()
        return

    def unwrap_fragments(self, fragments):
        if not self.periodic: return
        for f in fragments:
            xyz = self.xyz[f]
            xyz = self.pbc(xyz,0)
            self.xyz[f] = xyz
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
        self.xyz = np.dot(frac_xyz,self.cell)#

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
        if scale is None:
            scale = [1,1,1]
        if not hasattr(scale, '__iter__'):
            scale = [scale,scale,scale]
        self.cellparams *= np.hstack([scale,[1,1,1]])
        frac_xyz = self.get_frac_xyz()
        self.cell *= np.array(scale)[:,np.newaxis]
        self.images_cellvec = np.dot(images, self.cell)
        self.set_xyz_from_frac(frac_xyz)
        self.inv_cell = np.linalg.inv(self.cell)
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
                raise ValueError("can not add periodic systems with unequal cells!!")
                return
        other_xyz = other.xyz.copy()
        # NOTE: it is important ot keep the order of operations
        #       1) scale
        #       2) rotate by euler angles
        #       3) rotate by orientation triple
        #       4) translate
        if scale    is not None:
            other_xyz *= np.array(scale)
        if roteuler is not None:
            other_xyz = rotations.rotate_by_euler(other_xyz, roteuler)
        if rotate is not None:
            other_xyz = rotations.rotate_by_triple(other_xyz, rotate)
        if translate is not None:
            other_xyz += translate
        if self.natoms==0:
            self.xyz = other_xyz
        else:
            self.xyz = np.concatenate((self.xyz, other_xyz))
        self.elems += other.elems
        self.atypes+= other.atypes
        for c in other.conn:
            cn = (np.array(c)+self.natoms).tolist()
            self.conn.append(cn)
        self.natoms += other.natoms
        if len(other.fragtypes) == 0: other.set_nofrags()
        self.add_fragtypes(other.fragtypes)
        self.add_fragnumbers(other.fragnumbers)
        #self.fragtypes += other.fragtypes
        #start_fragnumber = sorted(self.fragnumbers)[-1]+1
        #self.fragnumbers += list(np.array(other.fragnumbers)+start_fragnumber)
        return

    def add_bond(self,a1,a2):
        """One-to-one connectivity: sets 1 bond between atom a1 and atom a2. Connectivity of both atoms
        is cross-updated by appending. (no sorting)
        :Parameter:
            -a1(int): index of atom1, python-like (starts with 0)
            -a2(int): index of atom2, python-like (starts with 0)"""
        if hasattr(a1,"__iter__"): a1=a1[0] #in case a singleton is passed
        if hasattr(a2,"__iter__"): a2=a2[0] #in case a singleton is passed
        self.conn[a1].append(a2)
        self.conn[a2].append(a1)
        return

    def add_bonds(self,lista1,lista2):
        """Many-to-many connectivity: Sets NxM  bonds, where N and M is the number of atoms per each list.
        Each atom of list 1 is connected to each atom of list 2.
        This is rarely wanted unless (at least) one of the lists has got only one atom.
        In that case, sets Nx1=N bonds, where N is the number of atoms of the "long" list.
        Each atom of the "long" list is connected to the atom of the "short" one.
        If lists have got just one atom per each, sets 1 bond (gracefully collapses to add_bond)
        between atom of list 1 and atom of list 2.
        :Paramters:
            -lista1(iterable of int): iterable 1 of atom indices
            -lista2(iterable of int): iterable 2 of atom indices"""
        if not hasattr(lista1,'__iter__'): lista1 = [lista1]
        if not hasattr(lista2,'__iter__'): lista2 = [lista2]
        for a1 in lista1:
            for a2 in lista2:
                self.add_bond(a1,a2)
        return

    def add_naive_hungarian_bonds(self,lista1,lista2):
        """Valid only in the 2x2 case, four times faster than standard hungarian method"""
        assert len(lista1) == len(lista2) == 2,\
            "only for 2x2 case, here: %dx%d case" % (len(lista1), len(lista2))
        a11, a12 = lista1
        a21, a22 = lista2
        d0 = self.get_distvec(a11,a21)
        d1 = self.get_distvec(a11,a22)
        if d1 > d0: #straight
            self.add_bond(a11,a21)
            self.add_bond(a12,a22)
        else: #cross
            self.add_bond(a11,a22)
            self.add_bond(a12,a21)
        return

    def add_standard_hungarian_bonds(self,lista1,lista2):
        dim = len(lista1)
        assert dim == len(lista2),\
            "only for NxN case (same number of atoms), here: %dx%d case" % (len(lista1), len(lista2))
        dmat = np.zeros([dim,dim])
        for e1,a1 in enumerate(lista1):
            for e2,a2 in enumerate(lista2):
                dmat[e1,e2] = self.get_distvec(a1,a2)[0]
        a1which, a2which = hungarian(dmat)
        for i in range(dim):
            self.add_bond(lista1[a1which[i]], lista2[a2which[i]])
        return

    def add_advanced_hungarian_bonds(self,lista1,lista2):
        raise NotImplementedError

    ###  molecular manipulations #######################################

    def delete_atoms(self,bads):
        ''' deletes an atom and its connections and fixes broken indices of all other atoms '''
        if hasattr(bads, '__iter__'):
            if len(bads) >= 2:
                self.bads = bads
                self.bads.sort()
                self.goods = [i for i in range(self.natoms) if i not in self.bads]
                self.offset = np.zeros(self.natoms, 'int')
                for i in range(self.natoms):
                    if i in self.bads:
                        self.offset[i:] += 1
                self.atypes = np.take(self.atypes,self.goods)
                self.elems  = np.take(self.elems, self.goods)
                self.conn   = np.take(self.conn,  self.goods)
                self.natoms = len(self.elems)
                self.conn =[ [j-self.offset[j] for j in self.conn[i] if j not in bads] for i in range(self.natoms) ]
                self.xyz    = self.xyz[self.goods]
                return
            else:
                if len(bads) != 0:
                    self.delete_atom(bads[0])
        else:
            self.delete_atom(bads)

    def delete_atom(self,bad):
        ''' deletes an atom and its connections and fixes broken indices of all other atoms '''
        new_xyz = []
        new_elems = []
        new_atypes = []
        new_conn = []
        for i in range(self.natoms):
            if i != bad:
                new_xyz.append(self.xyz[i].tolist())
                new_elems.append(self.elems[i])
                new_atypes.append(self.atypes[i])
                new_conn.append(self.conn[i])
                for j in range(len(new_conn[-1])):
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

    def remove_dummies(self,labels=['x','xx']):
        ''' removes atoms by atom labels
        :Parameters:
            - labels (list): atom labels to be removed'''
        badlist = []
        for i,e in enumerate(self.elems):
            if labels.count(e) != 0:
                badlist.append(i)
        logger.info('removing '+ str(badlist[::-1]))
        self.delete_atoms(badlist)
        #for i in badlist[::-1]: self.delete_atom(i)
        return
    
    def randomize_coordinates(self,maxdr=1.0):
        xyz = self.get_xyz()
        xyz += np.random.uniform(-maxdr,maxdr,xyz.shape)
        self.set_xyz(self.pbc(xyz))        

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
        center = self.get_com()
        self.translate(-center)
        return

    def get_com(self, idx = None, xyz = None):
        """
        returns the center of mass of the mol object.

        :Parameters:
            - idx  (list): list of atomindices to calculate the center of mass of a subset of atoms
        """
        if hasattr(self,'masstype') == False: self.set_real_mass()
        #if self.masstype == 'unit': logger.info('Unit mass is used for COM calculation')
        #if self.masstype == 'real': logger.info('Real mass is used for COM calculation')
        if xyz is not None:
            amass = np.array(self.amass)[idx]
        elif idx is None:
            if self.periodic: return None
            xyz = self.get_xyz()
            amass = np.array(self.amass)
        else:
            xyz = self.get_xyz()[idx]
            amass = np.array(self.amass)[idx]
        xyz = self.pbc(xyz, 0)
#        if self.periodic:
#            fix = xyz[0,:]
#            a = xyz[1:,:] - fix
#            if self.bcond <= 2:
#                cell_abc = self.cellparams[:3]
#                xyz[1:,:] -= cell_abc*np.around(a/cell_abc)
#            elif self.bcond == 3:
#                frac = np.dot(a, self.inv_cell)
#                xyz[1:,:] -= np.dot(np.around(frac),self.cell)
        center = np.sum(xyz*amass[:,np.newaxis], axis =0)/np.sum(amass)
        return center

    def map2image(self,xyz):
        if self.periodic == False: return xyz
        fix = xyz[0]
        a = xyz[1:,:] - fix
        if self.bcond <= 2:
            cell_abc = self.cellparams[:3]
            xyz[1:,:] -= cell_abc*np.around(a/cell_abc)
        elif self.bcond == 3:
            frac = np.dot(a, self.inv_cell)
            xyz[1:,:] -= np.dot(np.around(frac),self.cell)
        return xyz

    def pbc(self, xyz, fixidx = 0):
        """
        Compute periodic boundary conditions in an arbitrary (triclinic) cell
        """
        if self.periodic:
            fix = xyz[fixidx,:]
            a = xyz[:,:] - fix
            if self.bcond <= 2:
                cell_abc = self.cellparams[:3]
                xyz[:,:] -= cell_abc*np.around(a/cell_abc)
            elif self.bcond == 3:
                frac = np.dot(a, self.inv_cell)
                xyz[:,:] -= np.dot(np.around(frac),self.cell)
        return xyz

    def new_mol_by_index(self, idx):
        """
        Creates a new mol object which consists of the atoms specified in the argument.
        """
        m = mol()
        m.set_natoms(len(idx))
        d = {}
        elems = []
        xyz = []
        atypes = []
        for n,i in enumerate(idx):
            d[i] = n
            elems.append(self.elems[i])
            xyz.append(self.xyz[i,:])
            atypes.append(self.atypes[i])
        m.set_elems(elems)
        m.set_xyz(np.array(xyz))
        m.set_atypes(atypes)
        conn = []
        for i in idx:
            this_conn = []
            for j in self.conn[i]:
                try:
                    this_conn.append(d[j])
                except KeyError:
                    pass
            conn.append(this_conn)
        m.set_conn(conn)
        # handle periodic boundary conditions
        if type(self.cell) != type(None):
            m.set_cell(self.cell)
            m.periodic = True
            """ ###SOURCE OF BUG, YET NOT STUDIED
            stop = False
            while not stop:
                stop = True
                for i, conns in enumerate(m.conn):
                    for j in conns:
                        d, r, imgi = m.get_distvec(i, j)
                        if imgi != [13]:
                            stop = False
                            for ik, k in enumerate(self.cell):
                                m.xyz[j] += k * images[imgi][0][ik]
                            break
            """
            ### it SEEMS to work now without the while loop, NO WARRANTY (RA+MD)
            for i, conns in enumerate(m.conn):
                for j in conns:
                    d, r, imgi = m.get_distvec(i, j)
                    if imgi != [13]:
                        for ik, k in enumerate(self.cell):
                            m.xyz[j] += k * images[imgi][0][ik]
                        break
            m.cell = None
            m.periodic = False
        return m


    ##### distance measurements #####################

    def get_distvec(self, i, j, thresh=SMALL_DIST):
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
                # distance
                d_sort = d_sort[1:]
            closest = d_sort[0]
            closest=[closest]  # THIS IS A BIT OF A HACK BUT WE MAKE IT ALWAYS A LIST ....
            if (abs(all_d[closest[0]]-all_d[d_sort[1]]) < thresh):
                # oops ... there is more then one image atom in the same distance
                #  this means the distance is larger then half the cell width
                # in this case we have to return a list of distances
                for k in d_sort[1:]:
                    if (abs(all_d[d_sort[0]]-all_d[k]) < thresh):
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
    def add_atom(self, elem, atype, xyz):
        assert type(elem) == str
        assert type(atype)== str
        assert np.shape(xyz) == (3,)
        self.natoms += 1
        self.elems.append(elem)
        self.atypes.append(atype)
        xyz.shape = (1,3)
        if isinstance(self.xyz, np.ndarray):
            self.xyz = np.concatenate((self.xyz, xyz))
        else:
            self.xyz = xyz
        self.conn.append([])
        return self.natoms -1

    def add_conn(self, anum1, anum2):
        ''' add a bond between two atoms
            BEWARE, does not do any checks '''
        self.conn[anum1].append(anum2)
        self.conn[anum2].append(anum1)
        return

    def get_natoms(self):
        ''' returns the number of Atoms '''
        return self.natoms

    def set_natoms(self, natoms):
        """ sets the number of atoms for a new moltype """
        #assert self.natoms == 0
        self.natoms = natoms
        return

    def get_xyz(self):
        ''' returns the xyz Coordinates '''
        return self.xyz

    def set_xyz(self,xyz):
        ''' set the real xyz coordinates
        :Parameters:
            - xyz: coordinates to be set'''
        assert np.shape(xyz) == (self.natoms,3)
        self.xyz = xyz
        return

    def get_sumformula(self):
        """
        returns the sumformula of the mol object
        """
        fsum = ''
        unielems = sorted(list(set(self.elems)))
        elemscount = [self.elems.count(i) for i in unielems]
        for i,e in enumerate(unielems):
            fe = string.upper(e[0])+e[1:]
            fsum += fe
            fsum += str(elemscount[i])
        return fsum

    def get_elems(self):
        ''' return the list of element symbols '''
        return self.elems

    def get_elems_number(self):
        ''' return a list of atomic numbers '''
        return [elements.number[i] for i in self.elems]

    def get_elemlist(self):
        ''' Returns a list of unique elements '''
        el = []
        for e in self.elems:
            if not el.count(e): el.append(e)
        return el

    def set_elems(self,elems):
        ''' set the elements
        :Parameters:
            - elems: list of elements to be set'''
        assert len(elems) == self.natoms
        self.elems = elems

    def set_elems_number(self, elems_number):
        """ set the elemsnts from a list of atomic numbers ""
        :Parameters:
            - elem_number: list of atomic numbers
        """
        assert len(elems_number) == self.natoms
        self.elems = [elements.number.keys()[i] for i in elems_number]
        return

    def get_atypes(self):
        ''' return the list of atom types '''
        return self.atypes

    # just to make compatibel with pydlpoly standard API
    def get_atomtypes(self):
        return self.atypes

    def get_atypelist(self):
        ''' Returns a list of unique atom types '''
        if not self.atypes: return None
        return list(set(self.get_atypes()))

    def set_atypes(self,atypes):
        ''' set the atomtypes
        :Parameters:
            - atypes: list of elements to be set'''
        assert len(atypes) == self.natoms
        self.atypes = atypes

    def get_cell(self):
        ''' return unit cell information (cell vectors) '''
        return self.cell

    def get_cellparams(self):
        ''' return unit cell information (a, b, c, alpha, beta, gamma) '''
        return self.cellparams

    def set_bcond(self):
        """
        sets the boundary conditions. 2 for cubic and orthorombic systems,
        3 for triclinic systems
        """
        if list(self.cellparams[3:]) == [90.0,90.0,90.0]:
            self.bcond = 2
            if self.cellparams[0] == self.cellparams[1] == self.cellparams[2]:
                self.bcond = 1
        else:
            self.bcond = 3
        return

    def get_bcond(self):
        """
        returns the boundary conditions
        """
        return self.bcond

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
        self.inv_cell = np.linalg.inv(self.cell)
        self.images_cellvec = np.dot(images, self.cell)
        self.set_bcond()
        if cell_only == False: self.set_xyz_from_frac(frac_xyz)
        if not hasattr(self, "supercell"): self.supercell = [1,1,1]

    def set_cellparams(self,cellparams, cell_only = True):
        ''' set unit cell using cell parameters and assign cell vectors
        :Parameters:
            - cellparams: vector (6)
            - cell_only (bool)  : if false, also the coordinates are changed
                                  in respect to new cell
        '''
        assert len(list(cellparams)) == 6
        if cell_only == False: frac_xyz = self.get_frac_xyz()
        self.periodic = True
        self.cellparams = cellparams
        self.cell = unit_cell.vectors_from_abc(self.cellparams)
        self.inv_cell = np.linalg.inv(self.cell)
        self.images_cellvec = np.dot(images, self.cell)
        self.set_bcond()
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
        self.nfrags = sorted(self.fragnumbers)[-1]+1

    def get_nfrags(self):
        """
        returns the number of fragments in the actual system
        """
        return self.nfrags

    def add_fragnumbers(self,fragnumbers):
        """
        adds a set of fragnumbers to the actual system
        :Parameters:
            - fragnumbers: the fragment numbers to be set (list of integers)'''
        """
        self.fragnumbers += list(np.array(fragnumbers)+self.get_nfrags())
        self.nfrags = sorted(self.fragnumbers)[-1]  # changed this so they start at 1!
        return

    def add_fragtypes(self,fragtypes):
        """
        adds a set of fragntpyes to the actual system
        :Parameters:
            - fragtypes: the fragtypes to be set (list of strings)'''
        """
        self.fragtypes += fragtypes
        return

    def get_conn(self):
        ''' returns the connectivity of the system '''
        return self.conn

    def set_conn(self, conn, ctab_flag=False):
        ''' updates the connectivity of the system
        :Parameters:
            - conn    : List of lists describing the connectivity'''
        self.conn = conn
        if ctab_flag: self.ctab = self.get_conn_as_tab()

    def get_ctab(self):
        ''' returns the connectivity table (nbonds, 2)'''
        return self.ctab

    def set_ctab(self, ctab, conn_flag=False):
        ''' updates the connectivity table
        :Parameters:
            - ctab  : List of couples describing the connectivity'''
        self.ctab = ctab
        if conn_flag: self.set_conn_from_tab(ctab)

    def set_empty_conn(self):
        """
        sets an empty list of lists for the connectivity
        """
        self.conn = []
        for i in range(self.natoms):
            self.conn.append([])
        return
        
    def get_conn_as_tab(self, pconn_flag=None):
        """
        gets the connectivity as a table of bonds with shape (nbonds, 2)
        N.B.: can return ctab AND ptab if self.use_pconn == True
        """
        if pconn_flag is None: pconn_flag = getattr(self,"use_pconn",False)
        ctab = []
        ptab = []
        if pconn_flag:
            for i in range(self.natoms):
                ci = self.conn[i]
                pi = self.pimages[i]
                for j, pj in zip(ci,pi):
                    if j > i or (j==i and pj <= 13):
                        ctab.append([i,j])
                        ptab.append(pj)
            return ctab, ptab
        else:
            for i, ci in enumerate(self.conn):
                for j in ci:
                    if j > i:
                        ctab.append([i,j])
        return ctab
        
    def set_ctab_from_conn(self, pconn_flag=None):
        if pconn_flag is None: pconn_flag = getattr(self,"use_pconn",False)
        if pconn_flag:
            self.ctab, self.ptab = self.get_conn_as_tab(pconn_flag=True)
        else:
            self.ctab = self.get_conn_as_tab(pconn_flag=False)

    def set_conn_from_tab(self, ctab):
        """
        sets the connectivity froma table of bonds
        :Parameters:
            - ctab   : list of bonds (nbonds, 2)
        """
        self.set_empty_conn()
        for c in ctab:
            i,j = c
            self.conn[i].append(j)
            self.conn[j].append(i)
        return
        
    def set_unit_mass(self):
        """
        sets the mass for every atom to one
        """
        self.masstype = 'unit'
        self.amass = []
        for i in range(self.natoms):
            self.amass.append(1.0)
        return

    def set_real_mass(self):
        """
        sets the physical mass for every atom
        """
        self.masstype = 'real'
        self.amass = []
        for i in self.elems:
            self.amass.append(elements.mass[i])
        return

    def get_mass(self, return_masstype=False):
        """
        returns the mass for every atom as list
        """
        if return_masstype:
            return self.amass, self.masstype
        else:
            return self.amass

    def set_mass(self, mass, masstype='real'):
        """
        returns the mass for every atom as list
        """
        self.amass = mass
        self.masstype = masstype
        return

    def set_nofrags(self):
        ''' in case there are no fragment types and numbers, setup the data structure which is needed in some functions '''
        self.set_fragtypes(['-1']*self.natoms)
        self.set_fragnumbers([-1]*self.natoms)

    def get_comm(self):
        """ dummy call ... returns None ... for compatibility with pydlpoly system objects """
        return None

    def set_weight(self, weight):
        ''' sets the weight of the system
        :Parameters:
            - weight    : int/float'''
        self.weight = weight
        return

    def get_weight(self):
        ''' gets the weight of the system. Default: 1.'''
        return self.weight
