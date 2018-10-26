# -*- coding: utf-8 -*-
### overload print in parallel case (needs to be the first line) [RS] ###
from __future__ import print_function

import numpy as np
#from scipy.optimize import linear_sum_assignment as hungarian
#import types
import copy
import os
#import sys
import subprocess
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from .util import unit_cell
from .util.units import *
from .util import elems as elements
from .util import rotations
from .util import images
from .util.images import arr2idx, idx2arr, idx2revidx
from .util.misc import argsorted

from .fileIO import formats

from . import mpiobject
from . import molsys_mpi
from . import addon
from .prop import Property

import random

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
if molsys_mpi.size > 1:
    logger_file_name = "molsys.%d.log" % molsys_mpi.rank
else:
    logger_file_name = "molsys.log"
fhandler  = logging.FileHandler(logger_file_name)
fhandler.setLevel(logging.DEBUG)
#fhandler.setLevel(logging.WARNING)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m-%d %H:%M')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
if molsys_mpi.rank == 0:
    shandler  = logging.StreamHandler()
    shandler.setLevel(logging.INFO)
    #shandler.setLevel(logging.WARNING)
    shandler.setFormatter(formatter)
    logger.addHandler(shandler)

if molsys_mpi.wcomm is None:
    logger.error("MPI NOT IMPORTED DUE TO ImportError")
    logger.error(molsys_mpi.err)

# overload print function in parallel case
try:
    import __builtin__
except ImportError:
    import builtins as __builtin__
def print(*args, **kwargs):
    if molsys_mpi.rank == 0:
        return __builtin__.print(*args, **kwargs)
    else:
        return


np.set_printoptions(threshold=20000,precision=5)
SMALL_DIST = 1.0e-3


class mol(mpiobject):
    """mol class, the basis for any atomistic (or atomistic-like,
    e.g. topo) representation."""

    def __init__(self, mpi_comm = None, out = None):
        super(mol,self).__init__(mpi_comm, out)
        self.natoms=0
        self.nbonds=0
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
        # defaults
        self.is_topo = False # this flag replaces the old topo object derived from mol
        self.use_pconn = False # extra flag .. we could have toper that do not need pconn
        self.pconn = []
        self.pimages = []
        self.ptab  = []
        self.supercell=[1,1,1]
        self.aprops = {}
        self.bprops = {}
        self._etab = []
        return

    #####  I/O stuff ######################################################################################

    def set_logger_level(self,level='WARNING'):
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
        Parameters:
            fname(str)   : the filename to be read or a generic name
            ftype(str)   : the parser type that is used to read the file (default: "mfpx")
            **kwargs     : all options of the parser are passed by the kwargs
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
            logger.info('the file %s does not exist, trying with extension %s' % (fname,str(ftype)))
            try:
                f = open(fname+'.'+ftype, "r")
            except:
                raise IOError('the file %s does not exist' % (fname,))
        if ftype in formats.read:
            formats.read[ftype](self,f,**kwargs)
        else:
            logger.error("unsupported format: %s" % ftype)
            raise IOError("Unsupported format")
        f.close()
        return

    @classmethod
    def from_file(cls, fname, ftype=None, **kwargs):
        ''' reader for the mol class, reading from a file
        Parameters:
            fname(str): path to the file (filename included)
            ftype=None (or str): the parser type that is used to read the file
                if None: assigned by read as mfpx (default)
            **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        m = cls()
        m.read(fname, ftype, **kwargs)
        return m

    @classmethod
    def from_string(cls, istring, ftype='mfpx', **kwargs):
        ''' generic reader for the mol class, reading from a string
        Parameters:
            string       : the string to be read
            ftype="mfpx" : the parser type that is used to read the file
            **kwargs     : all options of the parser are passed by the kwargs
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
    def from_fileobject(cls, f, ftype='mfpx', **kwargs):
        ''' generic reader for the mol class, reading from a string
        Parameters:
            string       : the string to be read
            ftype="mfpx" : the parser type that is used to read the file
            **kwargs     : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        m = cls()
        logger.info("reading string as %s" % str(ftype))
        if ftype in formats.read:
            formats.read[ftype](m,f,**kwargs)
        else:
            logger.error("unsupported format: %s" % ftype)
            raise IOError("Unsupported format")
        return m

    

    @classmethod
    def from_abinit(cls, elems, xyz, cell, frac = False):
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
    def from_pymatgen(cls, structure):
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
    def from_ff(cls, basename, fit = False):
        m = cls()
        m.read(basename)
        m.addon("ff")
        m.ff.read(basename, fit = fit)
        return m


    @classmethod
    def from_array(cls, arr, **kwargs):
        ''' generic reader for the mol class, reading from a Nx3 array
        Parameters:
            arr         : the array to be read
            **kwargs    : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        m = cls()
        logger.info("reading array")
        assert arr.shape[1] == 3, "Wrong array dimension (second must be 3): %s" % (arr.shape,)
        formats.read['array'](m,arr,**kwargs)
        return m

    @classmethod
    def from_nested_list(cls, nestl, **kwargs):
        ''' generic reader for the mol class, reading from a Nx3 array
        Parameters:
            arr         : the array to be read
            **kwargs    : all options of the parser are passed by the kwargs
                             see molsys.io.* for detailed info'''
        logger.info("reading nested lists")
        for nl in nestl:
            assert len(nl) == 3, "Wrong nested list lenght (must be 3): %s" % (arr.shape,)
        arr = np.array(nestl)
        return cls.fromArray(arr, **kwargs)

    def to_phonopy(self, hessian = None):
        """
            Method to create a phonopy object for lattice dynamic calculations.

        Kwargs:
            hessian (numpy.ndarray, optional): Defaults to None. Hessian matrix of
                shape (3N,3N) in kcal/mol/A**2.
        
        Raises:
            ImportError: Raises Import Error when phonopy is not installed
        
        Returns:
            [Phonopy]: Return the phonopy object.
        """
        try:
            from phonopy import Phonopy
            from phonopy.structure.atoms import PhonopyAtoms
            from phonopy.units import ElkToTHz
        except:
            raise ImportError("Phonopy is not available!")
        assert self.periodic == True; "Requested system is not periodic!"
        unitcell = PhonopyAtoms(symbols = [i.title() for i in self.get_elems()],
            cell = self.get_cell(), scaled_positions = self.get_frac_from_xyz())
        # phonopy is setup by assuming atomic units for the hessian matrix
        phonon = Phonopy(unitcell, [[1,0,0],[0,1,0],[0,0,1]], factor = ElkToTHz)
        if hessian is not None:
            # we convert here the hessian to the phonopy format and  to atomic units
            hessian *= kcalmol/angstrom**2
            h2 = np.zeros((self.natoms,self.natoms,3,3), dtype = "double")
            for i in range(self.natoms):
                for j in range(self.natoms):
                    i3,j3 = 3*i, 3*j
                    h2[i,j,:,:]=hessian[i3:i3+3, j3:j3+3]
            phonon.set_force_constants(h2)
        return phonon

    def write(self, fname, ftype=None, **kwargs):
        ''' generic writer for the mol class
        Parameters:
            fname        : the filename to be written
            ftype="mfpx" : the parser type that is used to writen the file
            **kwargs     : all options of the parser are passed by the kwargs
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
                with open(fname,"w") as f:
                    formats.write[ftype](self,f,**kwargs)
            else:
                logger.error("unsupported format: %s" % ftype)
                raise IOError("Unsupported format")
        return

    def to_string(self, ftype='mfpx', **kwargs):
        """
        Method to output mol object as string in the format
        of the given filetype.

        Kwargs:
            ftype(string): name of the filetype, default to mfpx
        
        Raises:
            IOError
        
        Returns:
            string: mol object as string
        """
        f = StringIO()
        logger.info("writing string as %s" % str(ftype))
        if ftype in formats.read:
            formats.write[ftype](self,f,**kwargs)
        else:
            logger.error("unsupported format: %s" % ftype)
            raise IOError("Unsupported format")
        return f.getvalue()

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

    ##### addons ####################################################################################

    def addon(self, addmod, *args, **kwargs):
        """
        add an addon module to this object
        the adddon will be an instance of the addon class and available as attribute of mol instance

        Args:
            addmod (string): string name of the addon module
            *args          : positional arguments for the addon instantiator
            **kwargs       :    keyword arguments for the addon instantiator
        """
        loaded = False
        if addmod in self.loaded_addons:
            logger.warning("%s addon is already available as attribute of mol instance!" % addmod)
            return
        if addmod in addon.__all__: ### addon is enabled: try to set it
            if getattr(addon, addmod) is not None: ### no error raised during addon/__init__.py import
                try: ### get the addon attribute, initialize it and set as self attribute
                    addinit = getattr(addon, addmod)(self, *args, **kwargs)
                    setattr(self, addmod, addinit)
                except TypeError as e: ### HACK when 'from molsys.addon.addmod import something'
                    # in this case, e.g.: addon.ff is the MODULE, not the CLASS, so that we need TWICE
                    # the 'getattr' to get molsys.addon.ff.ff
                    addinit = getattr(getattr(addon, addmod),addmod)(self, *args, **kwargs)
                    setattr(self, addmod, addinit)
                except Exception as e: ### unexpected error! bugfix needed or addon used improperly
                    import traceback
                    traceback.print_exc()
                    logger.error("%s addon is not available: something unexpectable went wrong!" % addmod)
                else: ### kudos! the addon is now available as self.addmod
                    loaded = True
                    logger.info("%s addon is now available as attribute of mol instance" % addmod)
            else: ### error raised during addon/__init__.py import
                print(addon._errortrace[addmod])
                logger.error("%s addon is not imported: check addon module" % addmod)
        else: ### addon in unknown or disabled in addon.__all__
            logger.error("%s addon is unknown/disabled: check addon.__all__ in addon module" % addmod)
        if loaded:
            ### addmod added to loaded_addons (to prevent further adding)
            self.loaded_addons.append(addmod)
        #assert addmod in self.loaded_addons, "%s not available" % addmod ### KEEP for testing
        return

    ##### connectivity ########################################################################################

    def check_conn(self, conn=None):
        """
        checks if connectivity is not broken

        Args:
            conn (list): list of lists holding the connectivity (default=None, check own )
        """
        if conn is None:
            conn = self.conn
        for i, c in enumerate(conn):
            for j in c:
                logger.debug("%d in conn[%d] == %s? %s" % (i, j, conn[j], i in conn[j]))
                if i not in conn[j]: return False
        return True

    def detect_conn(self, tresh = 0.1,remove_duplicates = False, fixed_dist=False):
        """
        detects the connectivity of the system, based on covalent radii.

        Args:
            tresh (float): additive treshhold
            remove_duplicates (bool): flag for the detection of duplicates
            fixed_dist (bool or float, optional): Defaults to False. If a float is set this distance 
                replaces covalent radii (for blueprints use 1.0)
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
            dist = np.sqrt((a*a).sum(axis=1)) # distances from i to all other atoms
            conn_local = []
            if remove_duplicates == True:
                for j in range(i,natoms):
                    if i != j and dist[j] < tresh:
                        logger.debug("atom %i is duplicate of atom %i" % (j,i))
                        duplicates.append(j)
            else:
                for j in range(natoms):
                    if fixed_dist is False:
                        if i != j and dist[j] <= elements.get_covdistance([elems[i],elems[j]])+tresh:
                            conn_local.append(j)
                    else:
                        if i!= j and dist[j] <= fixed_dist+tresh:
                            conn_local.append(j)
            if remove_duplicates == False: conn.append(conn_local)
        if remove_duplicates:
            if len(duplicates)>0:
                logger.warning("Found and merged %d atom duplicates" % len(duplicates))
                duplicates = list(set(duplicates)) # multiple duplicates are taken once
                self.natoms -= len(duplicates)
                self.set_xyz(np.delete(xyz, duplicates,0))
                self.set_elems(np.delete(elems, duplicates))
                self.set_atypes(np.delete(self.atypes,duplicates))
                self.set_fragtypes(np.delete(self.fragtypes,duplicates))
                self.set_fragnumbers(np.delete(self.fragnumbers,duplicates))
            self.detect_conn(tresh = tresh)
        else:
            self.set_conn(conn)
        if self.use_pconn:
            # we had a pconn and redid the conn --> need to reconstruct the pconn
            self.add_pconn()
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
                print(d, 20*'#')
                self.pprint("   -> %3d %2s : dist %10.5f " % (conn[j], self.elems[conn[j]], d))
        return
    
    def add_pconn(self):
        """ generate the pconn from the exisiting conn
            The pconn contains the image index of the bonded neighbor
            pconn is really needed only for small unit cells (usually topologies) where vertices
            can be bonded to itself (next image) or multiple times to the same vertex in different images.
            """
        # N.B. add_pconn is not bullet-proof!
        # It works pretty always for large frameworks
        # whereas failing quite often with nets
        # (and quite always with the smaller ones)
        # JK+RA proposed FIX [work in progress, needs investigation]
        #for i in self.conn:
        #    i.sort()
        #
        # END FIX
        pimages = []
        pconn = []
        for i,c in enumerate(self.conn):
            atoms_pconn = []
            atoms_image = []
            for ji, j in enumerate(c):
                # If an atom or vertex is connected to another one multiple times (in an image), this
                # will be visible in the self.conn attribute, where the same neighbour will be listed
                # multiple times.
                # Sometimes, the distances are a bit different from each other, and in this case, we
                # have to increase the threshold, until the get_distvec function will find all imgis.
                n_conns = c.count(j)
                t = 0.01
                while True:
                    d,r,imgi = self.get_distvec(i,j,thresh=t)
                    t += 0.01
                    if n_conns == len(imgi):
                        break
                if len(imgi) == 1:
                    # only one neighbor .. all is fine
                    atoms_pconn.append(images[imgi[0]])
                    atoms_image.append(imgi[0])
                else:
                    # we need to assign an image to each connection
                    # if an atom is connected to another atom twice this means it must be another
                    # image
                    for ii in imgi:
                        # test if this image is not used for this atom .. then we can use it
                        if atoms_image.count(ii)==0:
                            atoms_image.append(ii)
                            atoms_pconn.append(images[ii])
                        else:
                            # ok, we have this image already
                            use_it = True
                            for k, iii in enumerate(atoms_image):
                                if (iii == ii) and (c[k] == j): use_it=False
                            if use_it:
                                atoms_image.append(ii)
                                atoms_pconn.append(images[ii])
            pimages.append(atoms_image)
            pconn.append(atoms_pconn)
        self.use_pconn= True
        self.pimages = pimages
        self.pconn = pconn
        return
    
    def check_need_pconn(self):
        """
        check whether pconn is needed or not
        """
        pconn_needed = False
        for i,c in enumerate(self.conn):
            # if atom/vertex bonded to itself we need pconn
            if i in c: pconn_needed = True
            # check if a neigbor appears twice
            for j in c:
                if c.count(j) > 1:
                    pconn_needed = True
            if pconn_needed: break
        return pconn_needed
    
    def omit_pconn(self):
        """
        Omit the pconn (if there is one) if this is acceptable
        """
        if not self.use_pconn: return
        if not self.check_need_pconn():
            # ok we do not need ot and can discard it
            self.pconn = []
            self.use_pconn = False
        return
    
    def make_topo(self, check_flag=True):
        """
        Convert this mol obejct to be a topo object.
        This means a pconn will be generated and written to file as well
        """
        if self.is_topo: return
        self.is_topo = True
        if self.check_need_pconn() and check_flag:
            self.add_pconn()
        return
    
    def unmake_topo(self):
        """
        Convert a topo object back to a "normal" mol object
        """
        if not self.is_topo: return
        self.is_topo = True
        self.omit_pconn()
        return
    

    ###  periodic systems .. cell manipulation ############
    
    #RS  Fix this
    #    currently if pconn is used we just call the old method from topo.py ... a lot of redundant things could 
    #    could be removed and all should be merged into one method at some point
    #

    def make_supercell(self,supercell):
        ''' Extends the periodic system in all directions by the factors given in the
            supercell upon preserving the connectivity of the initial system
            Can be used for systems with and without pconn
            Args:
                supercell: List of integers, e.g. [3,2,1] extends the cell three times in x and two times in y'''
        # HACK
        if self.use_pconn:
            xyz,conn,pconn = self._make_supercell_pconn(supercell)
            return xyz,conn,pconn
        # END HACK
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
        self.elems = list(self.elems)*ntot
        self.atypes=list(self.atypes)*ntot
        self.fragtypes=list(self.fragtypes)*ntot
        mfn = max(self.fragnumbers)+1
        nfragnumbers = []
        for i in range(ntot):
            nfragnumbers += list(np.array(self.fragnumbers)+i*mfn)
        self.fragnumbers=nfragnumbers
        self.images_cellvec = np.dot(images, self.cell)
        return xyz,conn


    def _make_supercell_pconn(self, supercell):
        """ old make_supercell from topo object
        called automatically when pconn exists
        """
        self.supercell = tuple(supercell)
        logger.info('Generating %i x %i x %i supercell' % self.supercell)
        img = [np.array(i) for i in images.tolist()]
        ntot = np.prod(supercell)
        nat = copy.deepcopy(self.natoms)
        nx,ny,nz = self.supercell[0],self.supercell[1],self.supercell[2]
        pconn = [copy.deepcopy(self.pconn) for i in range(ntot)]
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
                            if (img[13] == pconn[i][cc][c]).all():
                                #conn[i][cc][c] += ixyz*nat
                                conn[i][cc][c] = int( conn[i][cc][c] + ixyz*nat )
                                pconn[i][cc][c] = np.array([0,0,0])
                            else:
                                px,py,pz     = pconn[i][cc][c][0],pconn[i][cc][c][1],pconn[i][cc][c][2]
                                #print(px,py,pz)
                                iix,iiy,iiz  = (ix+px)%nx, (iy+py)%ny, (iz+pz)%nz
                                iixyz= iix+nx*iiy+nx*ny*iiz
                                conn[i][cc][c] = int( conn[i][cc][c] + iixyz*nat )
                                pconn[i][cc][c] = np.array([0,0,0])
                                if ((px == -1) and (left.count(ixyz)  != 0)): pconn[i][cc][c][0] = -1
                                if ((px ==  1) and (right.count(ixyz) != 0)): pconn[i][cc][c][0] =  1
                                if ((py == -1) and (bot.count(ixyz)   != 0)): pconn[i][cc][c][1] = -1
                                if ((py ==  1) and (top.count(ixyz)   != 0)): pconn[i][cc][c][1] =  1
                                if ((pz == -1) and (front.count(ixyz) != 0)): pconn[i][cc][c][2] = -1
                                if ((pz ==  1) and (back.count(ixyz)  != 0)): pconn[i][cc][c][2] =  1
                                #print(px,py,pz)
        self.natoms = nat*ntot
        self.conn, self.pconn, self.xyz = [],[],[]
        for cc in conn:
            for c in cc:
                self.conn.append(c)
        for pp in pconn:
            for p in pp:
                self.pconn.append(p)
        self.xyz = np.array(xyz).reshape(nat*ntot,3)
        self.cellparams[0:3] *= np.array(self.supercell)
        self.cell *= np.array(self.supercell)[:,np.newaxis]
        self.inv_cell = np.linalg.inv(self.cell)
        self.elems *= ntot
        self.atypes*=ntot
        self.images_cellvec = np.dot(images, self.cell)
        self.set_ctab_from_conn(pconn_flag=True)
        self.set_etab_from_tabs(sort_flag=True)
        return xyz,conn,pconn

    def apply_pbc(self, xyz=None, fixidx=0):
        ''' 
        apply pbc to the atoms of the system or some external positions
        Note: If pconn is used it is ivalid after this operation and will be reconstructed.
        
        Args:
            xyz (numpy array) : external positions, if None then self.xyz is wrapped into the box
            fixidx (int) : for an external system the origin can be defined (all atoms in one image). default=0 which means atom0 is reference, if fixidx=-1 all atoms will be wrapped
            
        Returns:
            xyz, in case xyz is not None (wrapped coordinates are returned) otherwise None is returned
        '''
        if not self.periodic:
            return xyz
        if xyz is None:
            # apply to structure itself (does not return anything)
            if self.bcond <= 2:
                cell_abc = self.cellparams[:3]
                self.xyz[:,:] -= cell_abc*np.around(self.xyz/cell_abc)
            elif self.bcond == 3:
                frac = self.get_frac_xyz()
                self.xyz[:,:] -= np.dot(np.around(frac),self.cell)
            if self.use_pconn:
                # we need to reconstruct pconn in this case
                self.add_pconn()
            return
        else:
            # apply to xyz
            assert xyz.ndim == 2, "number of dimensions must be 2"
            if fixidx != -1:
                a = xyz[:,:] - xyz[fixidx,:]
            else:
                a = xyz[:,:]
            if self.bcond <= 2:
                cell_abc = self.cellparams[:3]
                xyz[:,:] -= cell_abc*np.around(a/cell_abc)
            elif self.bcond == 3:
                frac = np.dot(a, self.inv_cell)
                xyz[:,:] -= np.dot(np.around(frac),self.cell)
        if self.use_pconn:
            self.add_pconn()
        return xyz

    # legacy name just to keep compat
    def wrap_in_box(self):
        """
        legacy method maps on apply_pbc
        """
        self.apply_pbc()
        return
    
    def get_cell(self):
        """get cell vectors
        
        Get the cell vectors as a 3x3 matrix, where the rows are the individual cell vectors
        
        Returns:
            numpy.ndarray: the cell matrix cell[0] or cell[0,:] is the first cell vector
        """
        return self.cell

    def get_cellparams(self):
        ''' return unit cell information (a, b, c, alpha, beta, gamma) '''
        return self.cellparams

    def get_volume(self):
        """returns volume of the cell
        
        Computes the Volume and returns it in cubic Angstroms
        
        Returns:
            float: Volume
        """
        cx = self.get_cell()
        return np.abs(np.dot(cx[0], np.cross(cx[1],cx[2])))


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
        Parameters:
            cell: cell vectors (3,3)
            cell_only (bool)  : if False, also the coordinates are changed
                                  in respect to new cell

        '''
        assert np.shape(cell) == (3,3)
        if cell_only is False: 
            frac_xyz = self.get_frac_from_xyz()
        self.periodic = True
        self.cell = cell
        self.cellparams = unit_cell.abc_from_vectors(self.cell)
        self.inv_cell = np.linalg.inv(self.cell)
        self.images_cellvec = np.dot(images, self.cell)
        self.set_bcond()
        if cell_only is False:
            self.set_xyz_from_frac(frac_xyz)
        return

    def set_cellparams(self,cellparams, cell_only = True):
        ''' set unit cell using cell parameters and assign cell vectors
        Parameters:
            cellparams: vector (6)
            cell_only (bool)  : if false, also the coordinates are changed
                                  in respect to new cell
        '''
        assert len(list(cellparams)) == 6
        cell = unit_cell.vectors_from_abc(cellparams)
        self.set_cell(cell, cell_only=cell_only)
        return
    
    def set_empty_cell(self):
        ''' set empty cell and related attributes'''
        self.periodic = False
        self.cell = None
        self.cellparams = None
        self.images_cellvec = None

    ### rewrite on set_cell ???
    def scale_cell(self, scale, cell_only=False):
        ''' scales the cell by a given factor
        
        Parameters:
            scale: either single float or an array of len 3'''
            
        cell = self.get_cell().copy()
        cell *= scale
        self.set_cell(cell, cell_only=cell_only)
        return

    def get_frac_xyz(self,xyz=None):
        return self.get_frac_from_xyz(xyz=xyz)
    
    def get_frac_from_xyz(self, xyz=None):
        ''' Returns the fractional atomic coordinates
        
        Parameters:
            xyz=None (array): optional external coordinates
        '''
        if not self.periodic: return None
        if xyz is None:
            xyz = self.xyz
        cell_inv = np.linalg.inv(self.cell)
        return np.dot(xyz, cell_inv)

    def get_xyz_from_frac(self,frac_xyz):
        ''' returns real coordinates from an array of fractional coordinates using the current cell info 
        
        Args:
            frac_xyz (array): fractional coords to be converted to xyz
        '''
        return np.dot(np.array(frac_xyz),self.cell)

    def set_xyz_from_frac(self, frac_xyz):
        ''' Sets atomic coordinates based on input fractional coordinates
        
        Args:
            - frac_xyz (array): fractional coords to be converted to xyz
        '''
        if not self.periodic: return
        assert frac_xyz.shape == (self.natoms, 3)
        self.xyz = np.dot(frac_xyz,self.cell)
        return

    def get_image(self,xyz, img):
        ''' returns the xyz coordinates of a set of coordinates in a specific cell
        Parameters:
            xyz   : xyz coordinates for which the image coordinates are to be retrieved
            img   : descriptor of the image, either an "images" integer (see molsys.util.images)
                      or the unit direction vector, e.g. [1,-1,0]'''
        xyz = np.array(xyz)
        try:
            l = len(img)
            dispvec = np.sum(self.cell*np.array(img)[:,np.newaxis],axis=0)
        except TypeError:
            dispvec = np.sum(self.cell*np.array(images[img])[:,np.newaxis],axis=0)
        return xyz + dispvec


    ###  add mol objects and copy ##########################################

    # TODO: this might not work anymore because of the mpi communicator
    #       is it needed?
    #def copy(self):
    #    ''' returns a copy of the whole mol object'''
    #    return copy.deepcopy(self)

    def add_mol(self, other, translate=None,rotate=None, scale=None, roteuler=None,rotmat=None):
        ''' adds a  nonperiodic mol object to the current one ... self can be both
            Parameters:
                other        (mol)      : an instance of the to-be-inserted mol instance
                translate (numpy.ndarry): numpy array as shift vector for the other mol
                rotate (numpy.ndarry)   : rotation triple to apply to the other mol object before insertion
                scale (float)           : scaling factor for other mol object coodinates
                roteuler (numpy.ndarry) : euler angles to apply a rotation prior to insertion'''
        if self.use_pconn:
            logger.info("Add mols with pconn, which may need tinkering")
        if other.periodic:
            if not (self.cell==other.cell).all():
                raise ValueError("can not add periodic systems with unequal cells!!")
                return
        other_xyz = other.xyz.copy()
        # NOTE: it is important ot keep the order of operations
        #       1 ) scale
        #       2 ) rotate by euler angles
        #       2a) rotate by rotmat
        #       3 ) rotate by orientation triple
        #       4 ) translate
        if scale    is not None:
            other_xyz *= np.array(scale)
        if roteuler is not None:
            other_xyz = rotations.rotate_by_euler(other_xyz, roteuler)
        if rotate is not None:
            other_xyz = rotations.rotate_by_triple(other_xyz, rotate)
        if rotmat is not None:
            other_xyz = np.dot(rotmat,other_xyz.T).T
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

    def new_mol_by_index(self, idx):
        """
        Creates a new mol object which consists of the atoms specified in thfe argument.

        Args:
            idx (list) : list of indices to be extracted as a new mol object
        """
        assert not self.use_pconn, "This method can not be used with pconn!"
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
            m.cellparams = None
            m.periodic = False
        return m

    ##### add and delete atoms and bonds ###########################################################
    
    def add_bond(self,idx1,idx2):
        ''' function necessary for legacy reasons! '''
        self.add_bonds(idx1,idx2)
        return

    def add_bonds(self, lista1, lista2, many2many=False):
        """ 
        add bonds/edges/connections to a mol object between exisiting atoms/vertices        

        If lists have got just one atom per each, sets 1 bond (gracefully collapses to add_bond)
        between atom of list 1 and atom of list 2.
        For many2many == False the length of lista1 and lista2 must be equal

        
        For many2many = True a Many-to-many connectivity is used:
        Sets NxM  bonds, where N and M is the number of atoms per each list.
        Each atom of list 1 is connected to each atom of list 2.
        This is rarely wanted unless (at least) one of the lists has got only one atom.
        In that case, sets Nx1=N bonds, where N is the number of atoms of the "long" list.
        Each atom of the "long" list is connected to the atom of the "short" one.
        
        Args:
            lista1(iterable of int): iterable 1 of atom indices
            lista2(iterable of int): iterable 2 of atom indices
            many2many (boolean):     switch to many2many mode            
            """
        if not hasattr(lista1,'__iter__'): lista1 = [lista1]
        if not hasattr(lista2,'__iter__'): lista2 = [lista2]
        if many2many == False:
            assert len(lista1)==len(lista2)
        if many2many:
            for a1 in lista1:
                for a2 in lista2:
                    self.conn[a1].append(a2)
                    self.conn[a2].append(a1)
                    if self.use_pconn:
                        d,v,imgi = self.get_distvec(a1,a2)
                        self.pconn[a1].append(images[imgi])
                        d,v,imgi = self.get_distvec(a2,a1)
                        self.pconn[a2].append(images[imgi])                
        else:
            for a1,a2 in zip(lista1, lista2):
                    self.conn[a1].append(a2)
                    self.conn[a2].append(a1)
                    if self.use_pconn:
                        d,v,imgi = self.get_distvec(a1,a2)
                        self.pconn[a1].append(images[imgi])
                        d,v,imgi = self.get_distvec(a2,a1)
                        self.pconn[a2].append(images[imgi])                
        return

    def add_shortest_bonds(self,lista1,lista2):
        """
        Adds bonds between atoms from list1 and list2 (same length!) to connect
        the shortest pairs

        in the 2x2 case, simple choice is used whereas for larger sets the hungarian method
        is used
        
        Args:
            lista1 (list) : list of atoms
            lista2 (list) : list of atoms
        
        """
        assert not self.use_pconn
        assert len(lista1) == len(lista2), "only for lists of same length: %dx != %d " % (len(lista1), len(lista2))
        if len(lista1) < 3:
            a11, a12 = lista1
            a21, a22 = lista2
            d0 = self.get_distvec(a11,a21)
            d1 = self.get_distvec(a11,a22)
            if d1 > d0: #straight
                self.add_bonds(a11,a21)
                self.add_bonds(a12,a22)
            else: #cross
                self.add_bonds(a11,a22)
                self.add_bonds(a12,a21)
        else:
            from scipy.optimize import linear_sum_assignment as hungarian
            dim = len(lista1)
            dmat = np.zeros([dim,dim])
            for e1,a1 in enumerate(lista1):
                for e2,a2 in enumerate(lista2):
                    dmat[e1,e2] = self.get_distvec(a1,a2)[0]
            a1which, a2which = hungarian(dmat)
            for i in range(dim):
                self.add_bonds(lista1[a1which[i]], lista2[a2which[i]])
        return

    def delete_bond(self, i, j):
        """delete bond between atom i and j

        Args:
            i (int): atom 1
            j (int): atom 2

        """
        idxj = self.conn[i].index(j)
        idxi = self.conn[j].index(i)
        self.conn[i].remove(idxj)
        self.conn[j].remove(idxi)
        if self.use_pconn:
            self.pconn[i].pop(idxj)
            self.pconn[j].pop(idxi)            
        return

    def add_atom(self, elem, atype, xyz, fragtype='-1', fragnumber=-1):
        """
        add a ato/vertex to the system (unconnected)
        
        Args:
            elem (string):    element symbol
            atype (string):   atom type string
            xyz (ndarry [3]): coordinates
            
        TBI: what to do with fragment types etc?
        """
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
        if self.use_pconn: self.pconn.append([])
        if self.fragtypes:
            self.fragtypes.append(fragtype)
            self.fragnumbers.append(fragnumber)
        return self.natoms-1

    def insert_atom(self,elem, atype, i, j, xyz=None):
        """Inserts an atom in a bond netween i and j

        Adds the atom at position between i and j or at xyz if given
        
        Args:
            elem (str): element
            atype (str): atomtype
            i (integer): atom 1
            j (integer): atom 2
            xyz (ndarray[3]): optional xyz position
        """
        if xyz is None:
            d, r, img = self.get_distvec(i,j)
            new_xyz = self.xyz[i] + 0.5* r[0]
        else:
            new_xyz = xyz
        new_atom = self.add_atom(elem, atype, new_xyz)
        self.delete_bond(i, j)
        self.add_bonds([i,j],[new_atom,new_atom])
        return

    def delete_atoms(self, bads, keep_conn=False):
        ''' 
        deletes an atom and its connections and fixes broken indices of all other atoms
        if keep_conn == True: connectivity is kept when atoms are in the middle of two others '''
        if hasattr(bads, '__iter__'):
            if len(bads) >= 2:
                bads.sort()
                goods = [i for i in range(self.natoms) if i not in bads]
                offset = np.zeros(self.natoms, 'int')
                for i in range(self.natoms):
                    if i in bads:
                        offset[i:] += 1
                self.atypes = np.take(self.atypes, goods).tolist()
                self.elems  = np.take(self.elems, goods).tolist()
                if keep_conn:
                    #works ONLY for edges: ERROR for terminal atoms and TRASH for the rest
                    if self.use_pconn: #must go before setting self.conn
                        self.pconn = [
                            [
                                # get image jp in i-th pconn if j-th atom of i-th conn not in bads
                                jp if self.conn[i][j] not in bads else
                                [
                                # else (j-th atom of i-th con in bads) get image kp
                                # in the pconn of j-th atom in ith conn
                                #    and get the first (TBI: all) among atoms associated to each kp different than i
                                    kp for k,kp in enumerate(self.pconn[self.conn[i][j]])
                                        if self.conn[self.conn[i][j]][k] != i
                                ][0] #works only for edges
                                for j,jp in enumerate(self.pconn[i])
                            ]
                            # if i-th atom not in bads
                            for i in range(self.natoms) if i not in bads
                        ]
                    self.conn = [
                        [
                            # subtract the j-th offset to atom j in i-th conn if j not in bads
                            j-offset[j] if j not in bads else
                            # else (j in bads) subtract the k-th offset to atom k in j-th conn
                            #    and get the first (TBI: all) among atoms different than i
                            [
                                k-offset[k] for k in self.conn[j] if k != i
                            ][0] #works only for edges
                            for j in self.conn[i]
                        ]
                        # if atom i not in bads
                        for i in range(self.natoms) if i not in bads
                    ]
                else:
                    if self.use_pconn: #must go before setting self.conn
                        self.pconn = [
                            [
                                # get image jp in i-th pconn if j-th atom of i-th conn not in bads
                                jp for j,jp in enumerate(self.pconn[i]) if self.conn[i][j] not in bads
                            ]
                            # if atom i not in bads
                            for i in range(self.natoms) if i not in bads
                        ]
                    self.conn = [
                        [
                            # subtract the j-th offset to atom j in i-th conn if j not in bads
                            j-offset[j] for j in self.conn[i] if j not in bads
                        ]
                        # if atom i not in bads
                        for i in range(self.natoms) if i not in bads
                    ]
                self.natoms = len(self.elems)
                self.xyz    = self.xyz[goods]
                return
            else:
                if len(bads) != 0:
                    self.delete_atom(bads[0], keep_conn=keep_conn)
        else:
            self.delete_atom(bads)

    def delete_atom(self,bad, keep_conn=False):
        """deletes an atom and its connections and fixes broken indices of all other atoms
        
        Args:
            bad (integer): atom index to remove
        """
        new_xyz = []
        new_elems = []
        new_atypes = []
        new_conn = []
        new_pconn = []
        for i in range(self.natoms):
            if i != bad:
                new_xyz.append(self.xyz[i].tolist())
                new_elems.append(self.elems[i])
                new_atypes.append(self.atypes[i])
                new_conn.append(self.conn[i])
                if self.use_pconn:
                    new_pconn.append(self.pconn[i])
                for j in range(len(new_conn[-1])):
                    if new_conn[-1].count(bad) != 0:
                        idx = new_conn[-1].index(bad)
                        new_conn[-1].pop(idx)
                        if self.use_pconn:
                            new_pconn[-1].pop(idx)
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
        self.pconn = new_pconn
        return

    def remove_dummies(self, labels=['x','xx'], keep_conn=False):
        ''' removes atoms by atom labels
        Args:
            labels (list): atom labels to be removed'''
        badlist = []
        for i,e in enumerate(self.elems):
            if labels.count(e) != 0:
                badlist.append(i)
        logger.info('removing '+ str(badlist[::-1]))
        self.delete_atoms(badlist, keep_conn=keep_conn)
        #for i in badlist[::-1]: self.delete_atom(i)
        return

    def remove_overlapping_atoms(self, thresh=SMALL_DIST):
        """
        remove atoms/vertices which are closer than thresh
        
        Note that it is unpredictable which atom is removed from the overlapping pair.         
        
        Args:
            thresh : distance threshold
        """
        badlist = []
        for i in range(self.natoms):
            for j in range(i+1, self.natoms):
                d,r,imgi=self.get_distvec(i,j)
                if d < thresh:
                    badlist.append(j)
        self.delete_atoms(badlist)
        return

    def merge_atoms(self, sele=None, parent_index=0, molecules_flag=False):
        """
        merge selected atoms
        sele(list of nested lists of int OR list of int): list of atom indices
        parent_index(int): index of parent atom in the selection which
            attributes are taken from (e.g. element, atomtype, etc.)
        molecules_flag(bool): if True: sele is regrouped accoring to the found
            molecules (e.g. if you select the COO of different linkers, each
            COO is merged per se). The same behavior can be reproduced with
            an appropriate nesting of sele, so consider molecules_flag a
            convenience flag.
        """
        if sele is None: # trivial if molecules_flag=False...
            sele = [range(self.natoms)]
        else:
            if not hasattr(sele[0], '__iter__'): # quick and dirt
                sele = [sele]
        assert len(set().union(*sele)) == len(sum(sele,[])),\
            "multiple occurring atom indices are NOT supported!"
        if molecules_flag:
            # atoms are merged per connected components i.e. molecules
            sele_molecules = []
            molidx = self.get_separated_molecules()
            for midx in molidx:
                for sel in sele:
                    msel = [i for i in midx if i in sel]
                    if msel != []:
                        sele_molecules.append(msel)
            sele = sele_molecules
        while True:
            try:
                sel = sele.pop(0)
            except IndexError:
                return
            else:
                xyz = self.xyz[sel].mean(axis=0)
                parent = sel[parent_index]
                elem = self.elems[parent]
                atype = self.atypes[parent]
                if self.fragtypes:
                    fragtype = self.fragtypes[parent]
                    fragnumber = self.fragnumbers[parent]
                    self.add_atom(elem, atype, xyz, fragtype=fragtype, fragnumber=fragnumber)
                else:
                    self.add_atom(elem, atype, xyz)
                conn_all = sum([self.conn[i] for i in sel],[])
                conn = set(conn_all) - set(sel)
                self.conn[-1] = conn
                for i in conn:
                    self.conn[i].append(self.natoms-1)
                if self.use_pconn:
                    raise NotImplementedError, "TBI! [RA]"
                    frac_xyz = self.get_frac_xyz()
                    frac_j = self.frac_xyz[-1]
                    for i in conn:
                        frac_i = self.frac_xyz[i]
                        a = (frac_j - frac_i)%[1,1,1]
                        xyz_i = self.xyz[i]
                        self.pconn[i].append()
                # offset trick (not new: see delete_atoms)
                # trick must be performed BEFORE delete_atoms!
                offset = np.zeros(self.natoms, 'int')
                for i in range(self.natoms):
                    if i in sel:
                        offset[i:] += 1
                # one of the last call, taking care of conn indices!
                # N.B.: offset must be initialized before delete_atoms
                self.delete_atoms(sel)
                # back to the trick
                for i,s in enumerate(sele):
                    sele[i] = [j-offset[j] for j in s]
        return # will never get it, here for clarity

    def get_separated_molecules(self, sele = None):
        """
        get lists of indices of atoms which are connected together inside the
        list and not connected outside the list.
        same as get islands (see toper) with a native graph-tools algorithm

        :Arguments:
        sele(list of int): selection list of atom indices
            if sele is None: find molecules in the whole mol
            else: find molecules just in the selection, countin non-connected
                atoms as separated molecules (e.g. if you select just the
                COO of a paddlewheel you get 4 molecules)

        >>> import molsys
        >>> m = molsys.mol.from_file("molecules.mfpx")
        >>> molecules_idx = m.get_separated_molecules()
        >>> for m_idx in molecules_idx:
        >>>     m.new_mol_by_index(m_idx).view()
        >>> # if in trouble: CTRL+Z and "kill %%"
        """
        try:
            from graph_tool.topology import label_components
        except ImportError:
            raise ImportError, "install graph-tool via 'pip install graph-tool'"
        from molsys.util.toper import molgraph
        from collections import Counter
        if sele is None:
            mg = molgraph(self)
        else:
            m = self.new_mol_by_index(sele)
            mg = molgraph(m)
        labels = label_components(mg.molg)[0].a.tolist()
        unique_labels = Counter(labels).keys()
        if sele is None:
            molidx = [[j for j,ej in enumerate(labels) if ej==i] for i in unique_labels]
        else:
            molidx = [[sele[j] for j,ej in enumerate(labels) if ej==i] for i in unique_labels]
        return molidx

    def shuffle_atoms(self, sele=None):
        """
        shuffle atom indices, debug purpose

        :Arguments:
        sele(list of int): selection list of atom indices
            if sele is None: all the atoms are shuffled

        many methods should be INVARIANT wrt. atom sorting
        N.B.: using numpy array since for readability
        """
        if sele is None:
            sele = range(self.natoms)
        sele_original = sele[:]
        random.shuffle(sele)
        # selection to original dictionary
        sele2sele_original = dict(zip(sele,sele_original))
        # coordinates #
        self.xyz[sele_original] = self.xyz[sele]
        # elements #
        elems = np.array(self.elems)
        elems[sele_original] = elems[sele]
        self.elems = [str(e) for e in elems.tolist()]
        # atomtypes #
        atypes = np.array(self.atypes)
        atypes[sele_original] = atypes[sele]
        self.atypes = [str(e) for e in atypes]
        # connectivity #
        conn = copy.deepcopy(self.conn)
        for i,ic in enumerate(conn):
            ic = [sele2sele_original[j] if j in sele else j for j in ic]
            conn[i] = ic
        conn = np.array(conn)
        conn[sele_original] = conn[sele][:]
        self.set_conn(conn.tolist())
        return sele2sele_original

### property interface #########################################################

    def get_atom_property(self, pname):
        return self.aprops[pname]

    def set_atom_property(self, pname):
        self[pname] = Property(pname, self.natoms, "atom")
        self.aprops[pname] = self[pname]
        return

    def del_atom_property(self, pname):
        del self.aprops[pname]
        del self[pname]
        return

    def list_atom_properties(self):
        if not self.aprops:
            print("No atom property")
            return
        print("Atom properties:")
        for prop in self.aprops:
            print(prop)
        return

    def get_bond_property(self, pname):
        return self.bprops[pname]

    def set_bond_property(self, pname):
        prop = Property(pname, self.nbonds, "bonds")
        setattr(self, pname, prop)
        self.bprops[pname] = getattr(self, pname)
        return

    def del_bond_property(self, pname):
        del self.bprops[pname]
        del self[pname]
        return

    def list_bond_properties(self):
        if not self.bprops:
            print("No bond property")
            return
        print("Bond properties:")
        for prop in self.bprops:
            print(prop)
        return

    def get_property(self, pname, ptype):
        if ptype.lower() == "atom":
            return self.get_atom_property(pname)
        elif ptype.lower() == "bond":
            return self.get_bond_property(pname)
        else:
            raise AttributeError("No \"%s\" property name: please use \"atom\" or \"bond\"" % pname)

    def set_property(self, pname, ptype):
        if ptype.lower() == "atom":
            self.set_atom_property(pname)
        elif ptype.lower() == "bond":
            self.set_bond_property(pname)
        else:
            raise AttributeError("No \"%s\" property name: please use \"atom\" or \"bond\"" % pname)
        return

    def del_property(self, pname, ptype):
        if ptype.lower() == "atom":
            self.del_atom_property(pname)
        elif ptype.lower() == "bond":
            self.del_bond_property(pname)
        else:
            raise AttributeError("No \"%s\" property name: please use \"atom\" or \"bond\"" % pname)
        return

    def list_properties(self):
        print("Properties:")
        self.list_atom_properties()
        self.list_bond_properties()
        return

    ##### manipulate geomtry #######################################################


    def randomize_coordinates(self,maxdr=1.0):
        """randomizes existing  coordinates
            maxdr (float, optional): Defaults to 1.0. maximum displacement
        """
        xyz = self.get_xyz()
        xyz += np.random.uniform(-maxdr,maxdr,xyz.shape)
        self.set_xyz(self.apply_pbc(xyz))

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

    def center_com(self, check_periodic=True, pbc = True):
        ''' centers the molsys at the center of mass '''
        if check_periodic:
            if self.periodic: return
        center = self.get_com(check_periodic=check_periodic, pbc = pbc)
        self.translate(-center)
        return

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

    def get_dist(self, ri, rj, thresh=SMALL_DIST):
        """ vector from i to j
        This is a tricky bit, because it is needed also for distance detection in the blueprint
        where there can be small cell params wrt to the vertex distances.
        In other words: i can be bonded to j multiple times (each in a different image)
        and i and j could be the same!!
        :Parameters':
            - i,j  : the indices of the atoms for which the distance is to be calculated"""
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - ri
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            d_sort = np.argsort(all_d)
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
            if self.use_pconn:
                img = self.pconn[i][ci]
                rj += np.dot(img, self.cell)
            else:
                all_rj = rj + self.images_cellvec
                all_r = all_rj - self.xyz[i]
                all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
                closest = np.argsort(all_d)[0]
                return all_rj[closest]
        return rj

    def get_neighb_coords_(self, i, j, img):
        """ returns coordinates of atom bonded to i which is ci'th in bond list
        TBI: merge get_neighb_coords and get_neighb_coords_
        :Parameters:
            - i   :  index of base atom
            - ci  :  index of bond atom
            - img :  cell image"""
        rj = self.xyz[j].copy()
        if self.periodic:
            if self.use_pconn:
                rj += np.dot(img, self.cell)
            else:
                all_rj = rj + self.images_cellvec
                all_r = all_rj - self.xyz[i]
                all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
                closest = np.argsort(all_d)[0]
                return all_rj[closest]
        return rj

    def get_neighb_dist(self, i, ci):
        """ returns the distance of atom bonded to i which is ci'th in bond list
        :Parameters:
            - i  :  index of the base atom
            - ci :  index of the conn entry of the ith atom"""
        ri = self.xyz[i]
        j = self.conn[i][ci]
        rj = self.xyz[j].copy()
        if self.periodic:
            if self.use_pconn:
                img = self.pconn[i][ci]
                rj += np.dot(img, self.cell)
            else:
                all_rj = rj + self.images_cellvec
                all_r = all_rj - self.xyz[i]
                all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
                closest = np.argsort(all_d)[0]
                return all_d[closest]
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

    def get_com(self, idx = None, xyz = None, check_periodic=True, pbc = True):
        """
        returns the center of mass of the mol object.

        Parameters:
            idx  (list): list of atomindices to calculate the center of mass of a subset of atoms
        """
        if hasattr(self,'masstype') == False: self.set_real_mass()
        #if self.masstype == 'unit': logger.info('Unit mass is used for COM calculation')
        #if self.masstype == 'real': logger.info('Real mass is used for COM calculation')
        if xyz is not None:
            amass = np.array(self.amass)[idx]
        elif idx is None:
            if self.periodic and check_periodic: return None
            xyz = self.get_xyz()
            amass = np.array(self.amass)
        else:
            xyz = self.get_xyz()[idx]
            amass = np.array(self.amass)[idx]
        if pbc: xyz = self.apply_pbc(xyz, 0)
        if np.sum(amass) > 0.0:
            center = np.sum(xyz*amass[:,np.newaxis], axis =0)/np.sum(amass)
        else: #every atom is dummy!
            center = np.sum(xyz,axis=0)
        return center

    def shift_by_com(self, alpha=2, **kwargs):
        """
        shift by center of mass
        alpha is needed otherwise atom distance is lost for excerpt of former
        periodic structures (e.g. a block)
        """
        ralpha = 1./alpha
        com = self.get_com(check_periodic=False, **kwargs)
        if self.periodic:
            shift = np.dot( np.dot(com, self.inv_cell)%ralpha, self.cell)
        else: # N.B.: reverse alpha has a different meaning
            shift = com*ralpha
        self.xyz -= shift
        return


    ###### get/set  core datastructures ###########################

    def get_natoms(self):
        ''' returns the number of Atoms '''
        return self.natoms

    def set_natoms(self, natoms):
        """ sets the number of atoms for a new moltype """
        #assert self.natoms == 0
        self.natoms = natoms
        return

    def get_xyz(self, idx=None):
        ''' returns the xyz Coordinates 
        
        Args:
            idx=None (list): optional list of indices
        '''
        if idx is None:
            return self.xyz
        else:
            return self.xyz[idx]

    def set_xyz(self,xyz, idx=None):
        ''' set the real xyz coordinates
        Args:
            xyz (array): coordinates to be set
            idx=None (list): optional list of indicies
        '''
        if idx is None:
            assert xyz.shape == (self.natoms,3)
            self.xyz = xyz
        else:
            assert xyz.shape == (len(idx), 3)
            self.xyz[idx] = xyz
        return

    def get_sumformula(self):
        """
        returns the sumformula of the mol object
        """
        fsum = ''
        unielems = sorted(list(set(self.elems)))
        elemscount = [self.elems.count(i) for i in unielems]
        for i,e in enumerate(unielems):
            fe = e[0].upper()+e[1:]
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
        """ set the elements from a list of atomic numbers ""
        :Parameters:
            - elem_number: list of atomic numbers
        """
        assert len(elems_number) == self.natoms
        self.elems = [elements.number.keys()[i] for i in elems_number]
        return

    def get_atypes(self):
        ''' return the list of atom types '''
        return self.atypes

    # just to make compatible with pydlpoly standard API
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
        self.nfrags = np.max(self.fragnumbers)+1 

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
        self.nfrags = np.max(self.fragnumbers)+1  
        return

    def add_fragtypes(self,fragtypes):
        """
        adds a set of fragntpyes to the actual system
        :Parameters:
            - fragtypes: the fragtypes to be set (list of strings)'''
        """
        self.fragtypes += fragtypes
        return

    ### CONNECTIVITY ###
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
                pi = self.pconn[i]
                for j, pj in zip(ci,pi):
                    if j > i or (j==i and arr2idx[pj] < 13):
                        ctab.append((i,j))
                        ptab.append(arr2idx[pj])
            return ctab, ptab
        else:
            for i, ci in enumerate(self.conn):
                for j in ci:
                    if j > i:
                        ctab.append((i,j))
        return ctab
        
    def set_ctab_from_conn(self, pconn_flag=None):
        if pconn_flag is None: pconn_flag = getattr(self,"use_pconn",False)
        if pconn_flag:
            self.ctab, self.ptab = self.get_conn_as_tab(pconn_flag=True)
        else:
            self.ctab = self.get_conn_as_tab(pconn_flag=False)

    def set_conn_from_tab(self, ctab):
        """
        sets the connectivity from a table of bonds
        :Parameters:
            - ctab   : list of bonds (nbonds, 2)
        """
        self.set_empty_conn()
        for c in ctab:
            i,j = c
            self.conn[i].append(j)
            self.conn[j].append(i)
        return
        
    ### PERIODIC CONNECTIVITY ###
    def get_pconn(self):
        ''' returns the periodic connectivity of the system '''
        return self.pconn

    def set_pconn(self, pconn, ptab_flag=False):
        ''' updates the periodic connectivity of the system
        :Parameters:
            - pconn    : List of lists describing the periodic connectivity'''
        self.pconn = pconn
        if ptab_flag: self.ptab = self.get_pconn_as_tab()

    def get_ptab(self):
        ''' returns the periodic connectivity table (nbonds, 2)'''
        return self.ptab

    def set_ptab(self, ptab, pconn_flag=False):
        ''' updates the periodic connectivity table
        :Parameters:
            - ptab  : List of couples describing the periodic connectivity'''
        self.ptab = ptab
        if pconn_flag: self.set_pconn_from_tab(ptab)

    def set_empty_pconn(self):
        """
        sets an empty list of lists for the periodic connectivity
        """
        self.pconn = []
        for i in range(self.natoms):
            self.pconn.append([])
        return
        
    def get_pconn_as_tab(self, pconn_flag=None):
        """
        gets the periodic connectivity as a table of bonds with shape (nbonds, 2)
        N.B.: can return ctab AND ptab if self.use_pconn == True
        """
        ### TBI ### [RA]
        raise NotImplementedError("Use get_conn_as_Tab w/ pconn_flag=True")
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
        
    def set_ptab_from_pconn(self, pconn_flag=None):
        raise NotImplementedError("Use set_ctab_from_conn w/ pconn_flag=True")
        # TBI: see acab for a suggested implementation [RA]
        if pconn_flag is None: pconn_flag = getattr(self,"use_pconn",False)
        if pconn_flag:
            self.ctab, self.ptab = self.get_conn_as_tab(pconn_flag=True)
        else:
            self.ctab = self.get_conn_as_tab(pconn_flag=False)

    def set_pconn_from_tab(self, ptab):
        """
        sets the periodic connectivity froma table of bonds
        :Parameters:
            - ptab   : list of peridioc images per bond (nbonds, 2)
        """
        assert hasattr(self, "ctab"), "ctab is needed for the method"
        self.set_empty_pconn()
        pconn = [[None for ic in c] for c in self.conn]
        for k,p in enumerate(ptab):
            i,j = self.ctab[k]
            ij = self.conn[i].index(j)
            ji = self.conn[j].index(i)
            pconn[i][ij] =  idx2arr[p]
            pconn[j][ji] = -idx2arr[p]
        self.pconn = pconn
        return

    def get_pimages(self):
        """
        return the indices of the periodic images of the system
        """
        return self.pimages

    def set_pimages(self, pimages, pconn_flag=False):
        """
        sets periodic image indices
        if pconn_flag: set periodic connectivity (arrays) from indices
        """
        self.pimages = pimages
        if pconn_flag:
            self.pconn =[[idx2arr[j] for j in pimagesi] for pimagesi in pimages]
        return

    def set_pimages_from_pconn(self, pconn=None):
        """
        sets periodic image indices from periodic connectivity (arrays)
        if pconn is None: get pconn from instance
        """
        if pconn is None:
            pconn = self.pconn
        self.pimages = [[arr2idx[j] for j in pconni] for pconni in pconn]
        return
        
    @property
    def etab(self):
        """ edge tab"""
        return self._etab

    @etab.setter
    def etab(self, etab):
        """any time edge tab is set, ctab, (ptab) and etab are sorted"""
        self._etab = etab
        self.nbonds = len(etab)
        self.sort_tabs(etab_flag=False)

    def set_etab_from_tabs(self, ctab=None, ptab=None, conn_flag=False, sort_flag=True):
        """set etab from ctab (and ptab). Both can be given or got from mol.
        if sort_flag: ctab, (ptab) and etab are sorted too."""
        if ctab is None and ptab is None:
            ctab = self.ctab
            ptab = self.ptab
        elif ctab is None or ptab is None:
            raise ValueError("ctab and ptab can't both be None")
        if self.use_pconn:
            etab_T = zip(*ctab)
            etab_T.append(ptab)
            self.etab = zip(*etab_T)
        else:
            self.etab = ctab[:]
        if sort_flag is True: #it sorts ctab, ptab, and etab too
            self.sort_tabs(etab_flag=False)
        if conn_flag is True:
            self.set_conn_from_tab(self.ctab)
            self.set_pconn_from_tab(self.ptab)
        return

    def set_etab(self, ctab=None, ptab=None, set_tabs=False):
        """set etab without sorting (sorting is default for etab.setter)"""
        if ctab is None and ptab is None:
            ctab = self.ctab
            ptab = self.ptab
        elif ctab is None or ptab is None:
            raise ValueError("ctab and ptab can't both be None")
        if self.use_pconn:
            etab = list(zip(*(zip(*ctab)+[ptab]))) # python3 compl.: zip iterator gets exhausted
        else:
            etab = ctab
        self._etab = etab

    def sort_tabs(self, etab_flag=False, conn_flag=False):
        """sort ctab, (ptab) and etab according to given convention
        Convention is the following:
            1)first ctab atom is lower or equal than the second
            2)if i and j are equal and there is pconn:
                revert the image to an index lower than 13 (the current cell)
        N.B. this sorting is stable.
        """
        ctab = self.ctab
        etab = self.etab
        if self.use_pconn:
            ptab = self.ptab
            for ii,(i,j) in enumerate(ctab):
                if i > j:
                    ctab[ii] = ctab[ii][::-1]
                    ptab[ii] = idx2revidx[ptab[ii]]
                if i == j and ptab[ii] > 13:
                    ptab[ii] = idx2revidx[ptab[ii]]
        else:
            for ii,(i,j) in enumerate(ctab):
                if i > j:
                    ctab[ii] = ctab[ii][::-1]
        asorted = argsorted(ctab)
        ctab_ = [ctab[i] for i in asorted]
        self.set_ctab(ctab_, conn_flag=conn_flag)
        if self.use_pconn:
            ptab_ = [ptab[i] for i in asorted]
            self.set_ptab(ptab_, pconn_flag=conn_flag)
        if etab_flag: # it ensures sorted etab and overwrites previous etab
            self.set_etab_from_tabs(sort_flag=False)
        elif etab:
            self._etab = [etab[i] for i in asorted]
        return

    ### UTILS ##################################################################
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
        ''' in case there are no fragment types and numbers, setup the data
        structure which is needed in some functions '''
        self.set_fragtypes(['-1']*self.natoms)
        self.set_fragnumbers([-1]*self.natoms)

    def get_comm(self):
        """ dummy call ... returns None ...
        for compatibility with pydlpoly system objects """
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
