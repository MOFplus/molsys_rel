# -*- coding: utf-8 -*-
# RS .. to overload print in parallel case (needs to be the first line)
from __future__ import print_function
"""
Created on Thu Mar 23 11:25:43 2017

@author: rochus


        addon module FF to implement force field infrastructure to the molsys
        
        contains class ric
        contains class ff

"""

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


import numpy as np
import uuid
import molsys
from molsys.util.timing import timer, Timer
from molsys.util import elems
from molsys.util.ff_descriptors import desc as ff_desc
from molsys.util.aftypes import aftype, aftype_sort
from molsys.util.ffparameter import potentials, varpars, varpar
from molsys.addon import base

import itertools
import copy
import string
try:
    import cPickle as Pickle
except ImportError:
    import _pickle as Pickle

import logging
import pdb
logger = logging.getLogger("molsys.ff")
#logger.setLevel(logging.INFO)
#import pdb; pdb.set_trace()

if mpi_comm is None:
    logger.error("MPI NOT IMPORTED DUE TO ImportError")
    logger.error(mpi_err)

class AssignmentError(Exception):

    def __init__(self, *args, **kwargs):
        Exception.__init__(self,*args,**kwargs)

    def to_dict(self):
        rv = {}
        rv["error"]="AssignmentError"
        rv["message"]="Set of parameters is incomplete"
        return rv

class ic(list):
    """
    Class inherited from flist which accepts attributes.
    Non-existing attributes return None instead of raising an error.
    """

    def __init__(self, *args, **kwargs):
        list.__init__(self,*args)
        for k,v in kwargs.items():
            setattr(self, k, v)
        self.used = False
        return
        
    def __getattr__(self, name):
        """
        if name is not an attribute return None instead of raising an error
        """
        if not name in self.__dict__:
            return None
        else:
            return self.__dict__[name]
            
    def to_string(self, width=None, filt=None, inc=0):
        """
        Method to generate a string representation of the ic.

        Kwargs:
            width (int): Defaults to None. Width of the string representation
            filt (str): Defaults to None. 
            inc (int): Defaults to 0.
        
        Returns:
            str: string representation of the ic
        """
        form = "%d "
        if width: form = "%%%dd " % width
        attrstring = ""
        if filt:
            for k in self.__dict__:
                if filt == "all" or k in filt:
                    if self.__dict__[k] != None:
                        attrstring += " %s=%s" % (k, self.__dict__[k])
        if inc!=0:
            values = tuple(np.array(list(self))+inc)
        else:
            values = tuple(self)
        return ((len(self)*form) % values) + attrstring


class ric:
    """
    Class to detect and keep all the (redundant) internal coordinates of the system.

    Args:
        mol(molsys.mol): Parent mol object in which the rics should be detected.
    """

    def __init__(self, mol):
        self.timer = Timer()        
        self._mol      = mol
        self.conn      = mol.get_conn()
        self.natoms    = mol.get_natoms()
        self.xyz       = mol.xyz
        self.aftypes   = []
        for i, a in enumerate(mol.get_atypes()):
            self.aftypes.append(aftype(a, mol.fragtypes[i]))
        return
        

    def find_rics(self, specials={"linear": [], "sqp":[]}, smallring = False):
        """
        Method to find all rics of the system

        Kwargs:
            specials(dict): Defaults to {"linear": [], "sqp":[]}. Dictionary of special atom types.
                which has to be included in the search.
            smallring(bool): Defaults to False. If true a smallring check for is performed.
        """
        self.bnd    = self.find_bonds()
        self.ang    = self.find_angles()
        self.oop    = self.find_oops()
        self.dih    = self.find_dihedrals(**specials)
        if smallring: self.check_smallrings()
        self.report()
        self.timer.write_logger(logger.info)
        return


    def set_rics(self, bnd, ang, oop, dih, sanity_test=True):
        """
        Method to set the rics from outsided. Args has to be be properly sorted lists of 
        ic objectos as if they are supplied by find_rics.

        Args:
            bnd(list): list of all bonds
            ang(list): list of all angles
            oop(list): list of all out-of planes
            dih(list): list of all dihedrals

        Kwargs:
            sanity_test(bool): Defaults to True. If True a sanity check is performed
                if the supplied RICs belong to the current system.
        """
        self.bnd = bnd
        self.ang = ang
        self.oop = oop
        self.dih = dih
        if sanity_test:
            # find bonds and check if they are equal ... this should be a sufficent test that the rest is the same, too
            
            mol_bnd = self.find_bonds()
            if mol_bnd != bnd:
                raise ValueError("The rics provided do not match the mol object!")
        return


    @timer("find bonds")
    def find_bonds(self):
        """
        Method to find all bonds of a system.
        
        Returns:
            list: list of internal coordinate objects defining all bonds
        """
        bonds=[]
        for a1 in range(self.natoms):
            for a2 in self.conn[a1]:
                if a2 > a1: bonds.append(ic([a1, a2]))
        return bonds

    def sort_angle(self, idx):
        """
        Method used for sorting angles, according to their
        aftypes.
        
        Args:
            idx (list): list of indices defining an angle
        
        Returns:
            list: sorted list of indeces defining the angle
        """
        if str(self.aftypes[idx[0]]) > str(self.aftypes[idx[2]]):
            idx.reverse()
        elif idx[0] > idx[2] and str(self.aftypes[idx[0]]) == str(self.aftypes[idx[2]]):
            idx.reverse()
        return idx

    @timer("find angles")
    def find_angles(self):
        """
        Method to find all angles of a system

        Returns:
            list: list of internal coordinate objects defining all angles
        """
        angles=[]
        for ca in range(self.natoms):
            apex_atoms = self.conn[ca]
            naa = len(apex_atoms)
            for ia in range(naa):
                aa1 = apex_atoms[ia]
                other_apex_atoms = apex_atoms[ia+1:]
                for aa2 in other_apex_atoms:
                    if str(self.aftypes[aa1]) < str(self.aftypes[aa2]):
                        angles.append(ic([aa1, ca, aa2]))
                    elif str(self.aftypes[aa1]) == str(self.aftypes[aa2]):
                        if aa1 < aa2:
                            angles.append(ic([aa1, ca, aa2]))
                        else:
                            angles.append(ic([aa2, ca, aa1]))
                    else:
                        angles.append(ic([aa2, ca, aa1]))
        return angles

    @timer("find oops")
    def find_oops(self):
        """
        Method to find all oops of a system
        
        Returns:
            list: list of internal coordinate objects defining all oops
        """
        oops=[]
        # there are a lot of ways to find oops ...
        # we assume that only atoms with 3 partners can be an oop center
        for ta in range(self.natoms):
            if (len(self.conn[ta]) == 3):
                # ah! we have an oop
                a1, a2, a3 = tuple(self.conn[ta])
                oops.append(ic([ta, a1, a2, a3]))
                #oops.append([ta, a2, a1, a3])
                #oops.append([ta, a3, a1, a2])
        return oops

    @timer("find dihedrals")
    def find_dihedrals(self, linear = [], sqp = []):
        """
        Method to find all dihedrals of a system

        Kwargs:
            linear(list): Defaults to []. List of linear atom types.
            sqp(list): Defaults to []. List of square planar atom types.

        Returns:
            list: list of internal coordinate objects defining all dihedrals
        """
        dihedrals=[]
        for a2 in range(self.natoms):
            for a3 in self.conn[a2]:
                # avoid counting central bonds twice
                if a3 > a2:
                    endatom1 = list(self.conn[a2])
                    endatom4 = list(self.conn[a3])
                    endatom1.remove(a3)
                    endatom4.remove(a2)
                    ### check if a3 or a2 is a linear one
                    lin = False
                    stubb = False
                    while self.aftypes[a2] in linear:
                        assert len(endatom1) == 1
                        lin = True
                        a2old = a2
                        a2 = endatom1[0]
                        endatom1 = list(self.conn[a2])
                        endatom1.remove(a2old)
                        ### we have now to check for stubbs
                        if len(endatom1) == 0:
                            stubb = True
                            break
                    if stubb: continue
                    while self.aftypes[a3] in linear:
                        assert len(endatom4) == 1
                        lin = True
                        a3old = a3
                        a3 = endatom4[0]
                        endatom4 = list(self.conn[a3])
                        endatom4.remove(a3old)
                        ### we have now to check for stubbs
                        if len(endatom1) == 0:
                            stubb = True
                            break
                    if stubb: continue
                    for a1 in endatom1:
                        con1 = list(self.conn[a1])
                        for a4 in endatom4:
                            ring = None
                            if a1 == a4: continue
                            if con1.count(a4):
                                ring = 4
                            else:
                                con4 = list(self.conn[a4])
                                for c1 in con1:
                                    if con4.count(c1):
                                        ring = 5
                                        break
                            d = ic([a1,a2,a3,a4], ring = ring)
                            if lin:
                                ### in case of a dihedral due to dihedral shifts,
                                ### it has to be checked if we have this dihedral already
                                if d not in dihedrals:
                                    dihedrals.append(d)
                            elif self.aftypes[a2] in sqp:
                                # calculate angle a1 a2 a3
                                if abs(self.get_angle([a1,a2,a3])-180.0) > 10.0:
                                    dihedrals.append(d)
                            elif self.aftypes[a3] in sqp:
                                # calculate angle a2 a3 a4
                                if abs(self.get_angle([a2,a3,a4])-180.0) > 10.0:
                                    dihedrals.append(d)
                            else:
                                dihedrals.append(d)
        return dihedrals

    @timer("find smallrings")
    def check_smallrings(self):
        """
        Method needed to check if smallrings are present and to mark the
        corresponding bonds and angles.
        """
        for d in self.dih:
            if d.ring:
                # bonds are sorted via atom indices
                self.bnd[self.bnd.index(sorted(d[0:2]))].ring = d.ring
                self.bnd[self.bnd.index(sorted(d[1:3]))].ring = d.ring
                self.bnd[self.bnd.index(sorted(d[2:4]))].ring = d.ring
                # angles are sorted in respect to aftypes, we have to
                # check both possibilities
                ang = self.sort_angle(d[0:3])              
                self.ang[self.ang.index(ang)].ring = d.ring
                ang = self.sort_angle(d[1:4])
                self.ang[self.ang.index(ang)].ring = d.ring
        return

    def get_distance(self,atoms):
        """
        Computes distance between two atoms.
        
        Parameters:
            atoms (list): list of atomindices

        Returns:
            float: atom dinstance
        """
        xyz = self._mol.apply_pbc(self.xyz[atoms],fixidx = 0)
        apex_1 = xyz[0]
        apex_2 = xyz[1]
        return np.linalg.norm(apex_1-apex_2)

    def get_angle(self,atoms):
        """
        Computes the angle between three atoms
        
        Parameters:
            atoms (list): list of atomindices

        Returns:
            float: angle in degree
        """
        xyz = self._mol.apply_pbc(self.xyz[atoms], fixidx = 0)
        #xyz = self.xyz[atoms]
        apex_1 = xyz[0]
        apex_2 = xyz[2]
        central = xyz[1]
        r1 = apex_1 - central
        r2 = apex_2 - central
        s = np.dot(r1,r2)/(np.linalg.norm(r1)*np.linalg.norm(r2))
        if s < -1.0: s=-1.0
        if s >  1.0: s=1.0
        phi = np.arccos(s)
        return phi * (180.0/np.pi)

    def get_multiplicity(self, n1, n2):
        """
        Routine to estimate the multiplicity of a dihedral from the local topology
        
        Parameters:
            n1 (int): number of connections of central atom 1
            n2 (int): number of connections of central atom 2
        
        Returns:
            int: multiplicity
        """
        assert type(n1) == type(n2) == int
        if   set([n1,n2])==set([5,5]): return 4
        elif set([n1,n2])==set([6,6]): return 4
        elif set([n1,n2])==set([3,6]): return 4 
        elif set([n1,n2])==set([4,4]): return 3
        elif set([n1,n2])==set([2,4]): return 3
        elif set([n1,n2])==set([3,4]): return 3
        elif set([n1,n2])==set([3,3]): return 2
        elif set([n1,n2])==set([2,3]): return 2
        elif set([n1,n2])==set([2,2]): return 1
        else:                          return None



    def get_dihedral(self, atoms,compute_multiplicity=True):
        """
        Computes dihedral angle between four atoms
        
        Args:
            atoms (list): list of atomindices

        Kwargs:
            compute_multiplicity (bool): Defaults to True. If True multiplicity
                is returned together with the angle.

        Returns:
            tuple: tuple of angle in degree and mutliplicity
        """
        xyz = self._mol.apply_pbc(self.xyz[atoms], fixidx = 0)
        apex1 = xyz[0]
        apex2 = xyz[3]
        central1 = xyz[1]
        central2 = xyz[2]
        b0 = -1.0*(central1-apex1)
        b1 = central2-central1
        b2 = apex2-central2
        n1 = np.cross(b0,b1)
        n2 = np.cross(b1,b2)
        arg = -np.dot(n1,n2)/(np.linalg.norm(n1)*np.linalg.norm(n2))
        if abs(1.0-arg) < 10**-14:
            arg = 1.0
        elif abs(1.0+arg) < 10**-14:
            arg = -1.0
        phi = np.arccos(arg)
        ### get multiplicity
        if compute_multiplicity : 
            m = self.get_multiplicity(len(self._mol.conn[atoms[1]]),
                len(self._mol.conn[atoms[2]]))
            return (phi * (180.0/np.pi), m)
        else:
            return (phi * (180.0/np.pi),0)

    def get_oop(self,atoms):
        """
        Dummy function to compute the value of an oop. It returns always 0.0.
        
        Parameters:
            atoms (list): list of atomindices

        Returns:
            float: always 0.0 is returned
        """
        return 0.0

    def compute_rics(self):
        """
        Computes the values of all rics in the system and attaches 
        them to the corresponding ic
        """
        for b in self.bnd: b.value = self.get_distance(list(b))
        for a in self.ang: a.value = self.get_angle(list(a))
        for d in self.dih: d.value = self.get_dihedral(list(d))
        for o in self.oop: o.value = self.get_oop(list(o))
        return

    def report(self):
        """
        Method to reports all found rics by the help of
        the logger object
        """
        logger.info("Reporting RICs")
        logger.info("%7d bonds"     % len(self.bnd))
        logger.info("%7d angles"    % len(self.ang))
        logger.info("%7d oops"      % len(self.oop))
        logger.info("%7d dihedrals" % len(self.dih))
        return



#class ff(molsys.base):
#class ff(base.base):
class ff(base):
    """Class to administrate the FF parameters and to assign them.

    The ff class is used to administrate the assignment of FF parameters
    following our specific assignment approach and to vary the parameters
    for FFgen runs.
    
    Args:
        molsys (molsys.mol): mol object for which parameters should be
            assigned.
    
    Kwargs:
        par(potentials): Defaults to None. External potentials object
            which should be used in this ff class. Only needed for multistructure
            fits with FFgen.

    Attr:
        ric_type()
        par()
        par_ind()
    """

    def __init__(self, mol, par = None):
        super(ff,self).__init__(mol)
#        self._mol = mol
        self.timer = Timer()
        self.ric = ric(mol)
        # defaults
        self.settings =  {
            "radfact" : 1.0,
            "radrule" : "arithmetic", 
            "epsrule" : "geometric",            
            }
        self.pair_potentials_initalized = False
        self.refsysname = None
        self.fit = False
        if par is not None: 
            assert type(par) == potentials
            self.par = par
        logger.debug("generated the ff addon")
        return

    def _init_data(self, cha=None, vdw=None):
        """
        Method to setup the internal data structres.

        Kwargs:
            cha(list): Defaults to None. If provided the
                indices in the list are setup as cha internal
                coordinates.
            vdw(list): Defaults to None. If provided the
                indices in the list are setup as vdw internal
                coordinates.
        """
        # make data structures . call after ric has been filled with data either in assign or after read
        # these are the relevant datastructures that need to be filled by one or the other way.
        if cha is None:
            cha = [ic([i]) for i in range(self._mol.natoms)]
        if vdw is None:
            vdw = [ic([i]) for i in range(self._mol.natoms)]
        self.ric_type = {
                "cha": cha,
                "vdw": vdw, 
                "bnd": self.ric.bnd, 
                "ang": self.ric.ang, 
                "dih": self.ric.dih, 
                "oop": self.ric.oop}
        self.parind = {
                "cha": [None]*self._mol.natoms,
                "vdw": [None]*self._mol.natoms,
                "bnd": [None]*len(self.ric.bnd),
                "ang": [None]*len(self.ric.ang),
                "dih": [None]*len(self.ric.dih),
                "oop": [None]*len(self.ric.oop),
                }
        return

    def _init_pardata(self, FF=None):
        if not hasattr(self,'par'):
            self.par = potentials({
                    "cha": {},
                    "vdw": {},
                    "vdwpr": {},
                    "bnd": {},
                    "ang": {},
                    "dih": {},
                    "oop": {},
                    })
            self.par.FF=FF

    @timer("assign parameter")
    def assign_params_offline(self, ref, key = False):
        loaded_pots = {'bnd':[],
                'ang':[],
                'dih':[],
                'oop':[],
                'vdw':[],
                'cha':[]}
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            for params in self.par[ic].values():
                pot = params[0]
                if pot not in loaded_pots[ic]: loaded_pots[ic].append(pot)
        with self.timer("find rics"):
            self.ric.find_rics(specials = {'linear':[]})
            self._init_data()
        # check here for assignment based on a ric par file
        # or on an old key file
        with self.timer("make atypes"):
            self.aftypes = []
            for i, a in enumerate(self._mol.get_atypes()):
                if not key:
                    self.aftypes.append(aftype(a, self._mol.fragtypes[i]))
                else:
                    self.aftypes.append(aftype(a, "leg"))
        with self.timer("parameter assignement loop"):
            for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
                for i, r in enumerate(self.ric_type[ic]):
                    if self.parind[ic][i] is None:
                        full_parname_list = []
                        aft_list = self.get_parname_sort(r,ic)
                        parname  = tuple(aft_list)
                        sparname = map(str,parname)
                        full_parname_list = []
                        for p in loaded_pots[ic]:
                            full_parname = p+"->("+string.join(sparname,",")+")|"+ref
                            if full_parname in self.par[ic] and full_parname not in full_parname_list:
                                full_parname_list.append(full_parname)
                        if full_parname_list != []:
                            self.parind[ic][i] = full_parname_list
        # check if all rics have potentials
        self.check_consistency()
        # check if there are potentials which are not needed -> remove them
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            # make a set
            badlist = []
            l = []
            for i in self.parind[ic]: l+=i
            s = set(l)
            for j in self.par[ic].keys():
                if j not in s:
                    badlist.append(j)
            # remove stuff from badlist
            for i in badlist: 
                del self.par[ic][i]
                #logger.warning("")
            #import pdb; pdb.set_trace()



    @timer("assign multi parameters")
    def assign_multi_params(self, FFs, refsysname=None, equivs={}, azone = [], special_atypes = {}, smallring = False, generic = None):
        """
        Method to orchestrate the parameter assignment fo the curent system using multiple force fields
        defined in FFs by getting the corresponding data from the webAPI.

        Args:
            FFs (list): list of strings containing the names of the FFs which should be assigned. Priority
                is defined by the ordering of the list

        Keyword Args:
            refsysname (string): if set this is a refsystem leading to special
                treatment of nonidentified params; defaults to None
            equivs (dict): dictionary holding information on equivalences,
                needed for parameterization processes in order to "downgrade" atomtpyes, defaults to
                {}
            azone (list): list of indices defining the active zone, defaults to []
            special_atypes (dict): dict of special atypes, if empty special atypes from
                mofplus are requested, defaults to {}

        """
        for i, ff in enumerate(FFs):
            # first element
            if i == 0: 
                self.assign_params(ff, refsysname = refsysname,equivs = equivs, azone = azone, 
                        special_atypes = special_atypes, consecutive = True, ricdetect = True, smallring=smallring)
            # last element
            elif i == len(FFs)-1:
                self.par.FF = ff
                self.assign_params(ff, refsysname = refsysname, equivs = equivs, azone = azone, 
                        special_atypes = special_atypes, consecutive = False, ricdetect = False, smallring=smallring, generic = generic)
            # in between
            else:
                self.par.FF = ff
                self.assign_params(ff, refsysname = refsysname, equivs = equivs, azone = azone, 
                        special_atypes = special_atypes,consecutive = True, ricdetect = False, smallring=smallring)
        return


                
    @timer("assign parameter")
    def assign_params(self, FF, verbose=0, refsysname=None, equivs = {}, azone = [], special_atypes = {}, 
            plot=False, consecutive=False, ricdetect=True, smallring = False, generic = None):
        """
        Method to orchestrate the parameter assignment for this system using a force field defined with
        FF getting data from the webAPI

        Args:

            - FF        :    [string] name of the force field to be used in the parameter search
            - verbose   :    [integer, optional] print info on assignement process to logger;
            defaults to 0
            - refsysname:    [string, optional] if set this is a refsystem leading to special 
            treatment of nonidentified params; defaults to None
            - equivs    :    [dict, optional] dictionary holding information on equivalences,
            needed for parameterization processes in order to "downgrade" atomtpyes, defaults to
            {}
            - azone     :    [list, optional] list of indices defining the active zone;
            defaults to []
            - special_atypes : [dict, optional] dict of special atypes, if empty special atypes from
            mofplus are requested, defaults to {}
        """
        assert type(equivs) == dict
        assert type(azone) == list
        self.refsysname = refsysname
        self.equivs = equivs
        self.active_zone = azone
        if self.refsysname is None and len(self.equivs.keys()) > 0:
            raise IOError("Equiv feature can only used together with a defined refsysname")
        if self.refsysname is None and len(self.active_zone) > 0:
            raise IOError("Azone feature can only used together with a defined refsysname")
        with self.timer("connect to DB"):
            ### init api
            if self._mol.mpi_rank == 0:
                from mofplus import FF_api
                self.api = FF_api()
                if len(special_atypes) == 0: special_atypes = self.api.list_special_atypes()
            else:
                self.api = None
                special_atypes = None
            if self._mol.mpi_size > 1:
                special_atypes = self._mol.mpi_comm.bcast(special_atypes, root = 0)
        if ricdetect==True:
            with self.timer("find rics"):
                self.ric.find_rics(specials = special_atypes, smallring = smallring)
                self._init_data()
                self._init_pardata(FF)
            # as a first step we need to generate the fragment graph
            self.timer.start("fragment graph")
            self._mol.addon("fragments")
            self.fragments = self._mol.fragments
            self.fragments.make_frag_graph()
            if plot:
                self.fragments.plot_frag_graph(plot, ptype="png", vsize=20, fsize=20, size=1200)
            # create full atomistic graph
            self._mol.graph.make_graph()
            self.timer.stop()
            # now make a private list of atom types including the fragment name
            self.timer.start("make atypes")
            self.aftypes = []
            for i, a in enumerate(self._mol.get_atypes()):
                self.aftypes.append(aftype(a, self._mol.fragtypes[i]))
            self.timer.stop()
        # detect refsystems
        self.find_refsystems(plot=plot)
        with self.timer("parameter assignement loop"):
            for ref in self.scan_ref:
                counter = 0
                logger.info("assigning params for ref system %s" % ref)
                curr_fraglist = self.ref_fraglists[ref]
                curr_atomlist = self.ref_atomlists[ref]
                curr_par = {\
                    "bnd" : self.ref_params[ref]["twobody"]["bnd"],\
                    "bnd5" : self.ref_params[ref]["twobody"]["bnd5"],\
                    "ang" : self.ref_params[ref]["threebody"]["ang"],\
                    "ang5" : self.ref_params[ref]["threebody"]["ang5"],\
                    "dih" : self.ref_params[ref]["fourbody"]["dih"],\
                    "oop" : self.ref_params[ref]["fourbody"]["oop"],
                    "cha" : self.ref_params[ref]["onebody"]["charge"],
                    "vdw" : self.ref_params[ref]["onebody"]["vdw"]
                    }
                curr_equi_par = {}
                for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
                    if verbose>0: logger.info(" ### Params for %s ###" % ic)
                    for i, r in enumerate(self.ric_type[ic]):
                        if self.parind[ic][i] is None:
                            if ((self.atoms_in_subsys(r, curr_fraglist)) and (self.atoms_in_active(r, curr_atomlist))):
                                # no params yet and in current refsystem => check for params
                                full_parname_list = []
                                aft_list = self.get_parname_equiv(r,ic,ref)
                                #aft_list = map(lambda a: self.aftypes[a], r)
                                 # generate list of permuted tuples according to ic and look up params
                                #parname, par_list = self.pick_params(aft_list, ic, curr_par[ic])
                                parname, par_list = self.pick_params(aft_list, ic, r, curr_par)
                                if par_list != None:
                                    if verbose>1 : logger.info(" found parameter for atoms %20s (types %s) -> %s" % (str(r), aft_list, parname))
                                    for par in par_list:
                                        ### check for equivalences
                                        if par[0] == "equiv":
                                            for j, aft in enumerate(aft_list):
                                                aidx = r[j]
                                                if ((str(aft) == par[1][0]) and (aidx not in curr_equi_par)):
                                                    curr_equi_par[aidx] = par[1][1]
                                                    logger.info("  EQIV: atom %d will be converted from %s to %s" % (aidx, aft, par[1][1]))
                                        else:
                                            sparname = map(str, parname)
                                            full_parname = par[0]+"->("+string.join(sparname,",")+")|"+ref
                                            full_parname_list.append(full_parname)
                                            if not full_parname in self.par[ic]:
                                                logger.info("  added parameter to table: %s" % full_parname)
                                                self.par[ic][full_parname] = par
                                else:
                                    if verbose>1 : logger.info(" NO parameter for atoms %20s (types %s) " % (str(r), aft_list))
                                if full_parname_list != []:
                                    counter += 1
                                    self.parind[ic][i] = full_parname_list
                                #else:
                                #    print("DEBUG DEBUG DEBUG %s" % ic)
                                #    print(self.get_parname(r))
                                #    print(self.get_parname_sort(r, ic))
                logger.info("%i parameters assigned for ref system %s" % (counter,ref))
                #EQUIVALENCE
                # now all params for this ref have been assigned ... any equivalnce will be renamed now in aftypes
                for i, a in enumerate(copy.copy(self.aftypes)):
                    if i in curr_equi_par.keys():
                        at, ft = curr_equi_par[i].split("@")
                        self.aftypes[i] = aftype(at,ft)
        if consecutive==False:
            if refsysname:
                self.fixup_refsysparams()
            else:
                self.check_consistency(generic = generic)
        self.timer.write_logger(logger.info)
        return

    def check_consistency(self, generic = None):
        """
        Method to check the consistency of the assigned parameters, if params are
        missing an AssignmentError is raised.

        Kwargs:
            generic(str): name of the FF which should be used as source for
                the generic parameters, defaults to None
        """
        if generic is not None:
            gen_params = self.api.get_params_from_ref(generic, "generic")["fourbody"]
        complete = True
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            unknown_par = []
            for i, p in enumerate(self.ric_type[ic]):
                if self.parind[ic][i] is None:
                    # if possible we apply here a generic parameter
                    match = False
                    if generic is not None and ic == "dih":
                        aft_list = self.get_parname_equiv(p,ic,"generic")
                        full_parname_list = []
                        parname, par_list = self.pick_params(aft_list, ic, p, gen_params)
                        if par_list != None:
                            for par in par_list:
                                sparname = map(str, parname)
                                full_parname = par[0]+"->("+string.join(sparname,",")+")|generic"
                                full_parname_list.append(full_parname)
                                if not full_parname in self.par[ic]:
                                    logger.info("  added parameter to table: %s" % full_parname)
                                    self.par[ic][full_parname] = par
                        if full_parname_list != []:
                            self.parind[ic][i] = full_parname_list
                            match = True
                    # no params found
                    if match == False:
                        parname = self.get_parname_sort(p,ic)
                        if not parname in unknown_par:
                            unknown_par.append(parname)
            if len(unknown_par) > 0:
                complete = False
                for p in unknown_par: logger.error("No params for %3s %s" % (ic, p))
        if complete == False:
            raise AssignmentError("Assignend parameter set incomplete!")
        else:
            logger.info("Parameter assignment successfull")
        return

    def fixup_refsysparams(self, var_ics = ["bnd", "ang", "dih", "oop"], strbnd = False):
        """
        Equivalent method to check consistency in the case that unknown parameters should
        be determined eg. fitted. The Method prepares the internal data structures for
        parameterization or sets default values for unkown params.
        :Parameters:
            - var_ics(list, optional): list of strings of internal coordinate names
            for which the fixup should be done, defaults to ["bnd", "ang", "dih", "oop"]
            - strbnd(bool, optional): switch for a forcing a fixup of strbnd terms, defauls to
            False
        """
        self.ric.compute_rics()
        self.par.attach_variables()
#        self.variables = varpars()
        if hasattr(self, "active_zone") == False:
            self.active_zone = []
        defaults = {
            "bnd" : ("mm3", 2, "b", ["d","r"], [0.0, 8.0]),
            "ang" : ("mm3", 2, "a", ["d","r"], [0.0, 2.0]),
            "dih" : ("cos3", 3, "d", ["d","d","d"], [0.0, 15.0]),
            "oop" : ("harm", 2, "o", ["d",0.0], [0.0, 1.0]),
            "cha" : ("gaussian", 2, "c", ["d","d"], [0.0, 2.0]),
            "vdw" : ("buck6d", 2, "v", ["d","d"], [0.0, 2.0])}
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            count  = 0
            ric = self.ric_type[ic]
            par = self.par[ic]
            parind = self.parind[ic]
            for i, p in enumerate(ric):
                if parind[i] is None:
                    if ic == "cha" and i not in self.active_zone:
                        self.active_zone.append(i)
                    # not sure if we should sort here or not ... maybe not?
                    # HACK : sort all but angles because here we need strbnd 
                    if ic == "ang":
                        parname = self.get_parname(p)
                    else:
                        parname = self.get_parname_sort(p, ic)
                    sparname = map(str, parname)
                    fullparname = defaults[ic][0]+"->("+string.join(sparname,",")+")|"+self.refsysname
                    ### we have to set the variables here now
                    if not fullparname in par:
                        if ic in var_ics:
                            count+=1
                            vnames = map(lambda a: "$%s%i_%i" % (defaults[ic][2],count,a) 
                                if type(defaults[ic][3][a]) == str else defaults[ic][3][a], range(defaults[ic][1]))
                            par[fullparname] = (defaults[ic][0], vnames)
                            for idx,vn in enumerate(vnames):
                                if type(vn) == str:
                                    if defaults[ic][3][idx] == "r":
                                        self.par.variables[vn]=varpar(self.par,name = vn, 
                                                val = p.value, range = [0.9*p.value, 1.1*p.value])
                                    else:
                                        self.par.variables[vn]=varpar(self.par,name = vn, range = defaults[ic][4])
                                    self.par.variables[vn].pos.append((ic, fullparname, idx))
                            # hack for strbnd
                            if ic == "ang" and strbnd == True:
                                fullparname2 = "strbnd->("+string.join(sparname,",")+")|"+self.refsysname
                                count+=1
                                vnames = map(lambda a: "$a%i_%i" % (count, a), range(6))
                                par[fullparname2] = ("strbnd", vnames)
                                for idx,vn in enumerate(vnames):
                                    self.par.variables[vn] = varpar(self.par, name = vn)
                                    self.par.variables[vn].pos.append((ic,fullparname2,idx))
                        else:
                            par[fullparname] = [defaults[ic][0], defaults[ic][1]*[0.0]]
                    if ic == "ang" and strbnd == True:
                        parind[i] = [fullparname, fullparname2]
                    else:
                        parind[i] = [fullparname]
        self.set_def_sig(self.active_zone)
        self.set_def_vdw(self.active_zone)
        self.fix_strbnd()
        return

    def fix_strbnd(self):
        """
        Method to perform the fixup for strbnd potentials
        """
        ### get potentials to fix
        pots = self.par.variables.varpots
        dels = []
        for p in pots:
            pot, ref, aftypes = self.split_parname(p[1])
            if pot == "strbnd":
                # first check if apex atypes are the same
                if aftypes[0] == aftypes[2]:
                    dels.append(self.par["ang"][p[1]][1][1])
                    self.par["ang"][p[1]][1][1] = self.par["ang"][p[1]][1][0]
                # now distribute ref values
                apot  = "mm3->"+p[1].split("->")[-1] 
                spot1 = self.build_parname("bnd", "mm3", self.refsysname, aftypes[:2])
                spot2 = self.build_parname("bnd", "mm3", self.refsysname, aftypes[1:])
                s1 = self.par["bnd"][spot1][1][1]
                s2 = self.par["bnd"][spot2][1][1]
                a  = self.par["ang"][apot][1][1]
                # del variables
                dels.append(self.par["ang"][p[1]][1][3])
                dels.append(self.par["ang"][p[1]][1][4])
                dels.append(self.par["ang"][p[1]][1][5])
                # rename variables
                self.par["ang"][p[1]][1][3] = s1
                self.par["ang"][p[1]][1][4] = s2
                self.par["ang"][p[1]][1][5] = a
                # redistribute pots to self.variables dictionary
                self.par.variables[s1].pos.append(("ang", p[1],3))
                self.par.variables[s2].pos.append(("ang", p[1],4))
                self.par.variables[a].pos.append(("ang", p[1],5))
        for i in dels: del(self.par.variables[i])




    def set_def_vdw(self,ind):
        """
        Method to set default vdw parameters for a given atom index
        :Parameters:
            - ind(int): atom index
        """
        elements = self._mol.get_elems()
        atypes   = self._mol.get_atypes()
        truncs   = [i.split("_")[0] for i in atypes]
        for i in ind:
            elem  = elements[i]
            at    = atypes[i]
            trunc = truncs[i]
            try:
                prm = elems.vdw_prm[at]
            except:
                try:
                    prm = elems.vdw_prm[trunc]
                except:
                    prm = elems.vdw_prm[elem]
            self.par["vdw"][self.parind["vdw"][i][0]][1] = prm
        return


    def set_def_sig(self,ind):
        """
        Method to set default parameters for gassuain width for the the charges
        :Parameters:
            - ind(int): atom index
        """
        parind = self.parind["cha"]
        fitdat = {"fixed":{}, "equivs": {}, "parnames": []}
        ### first gather already assigned charges and dump them to fitdat["fixed"]
        for i,p in enumerate(parind):
            fitdat["parnames"].append(p[0])
            if i not in ind:
                fitdat["fixed"][i] = self.par["cha"][p[0]][1][0]
        elements = self._mol.get_elems()
        parnames = {}
        for i in ind:
            elem = elements[i]
            parname = parind[i][0]
            if parname in parnames.keys():
                parnames[parname].append(i)
            else:
                parnames[parname]=[i]
            try:
                sig = elems.sigmas[elem]
            except:
                sig = 0.0            
            self.par["cha"][self.parind["cha"][i][0]][1][1] = sig
        ### move stuff from parnames to fitdat["equivs"]
        for k,v in parnames.items():
            for i in v[1:]:
                fitdat["equivs"][i] = v[0]
        ### dump fitdat to json file
        with open("espfit.pickle", "wb") as f: Pickle.dump(fitdat, f)
        return

    def varnames2par(self):
        """
        Forces the paramters in the variables dictionary to be wriiten in the internal
        data structures
        """
        if hasattr(self, 'do_not_varnames2par'): return
        self.par.variables(self.par.variables.keys())
        return

    def remove_pars(self,identifier=[]):
        """Remove Variables to write the numbers to the fpar file
        
        Based on any identifier, variables can be deleted from the dictionary
        Identifiers can be regular expressions
        Usage:
            ff.remove_pars(['d'])  -- removes all dihedral variables
            ff.remove_pars(['a10'])  -- removes all entries of angle 10, e.g. a10_0 & a10_1
            ff.remove_pars(['b*_0'])  -- removes all bond entries zero e.g. b1_0 & b5_0, but keeps b1_1 and b5_1

        Keyword Arguments:
            identifier {list of strings} -- string identifiers, can be regexp (default: {[]})
        """
        import re
        for ident in identifier:
            re_ident = re.compile(ident)
            for k in self.par.variables.keys():
                if k.count(ident) != 0:
                    del self.par.variables[k]
                    continue
                # check regexp
                re_result = re.search(re_ident,k)
                if re_result is None: continue
                if re_result.span()[-1] != 0:  #span[1] is zero if the regexp was not found
                    del self.par.variables[k]
                    continue
        return

    def setup_pair_potentials(self):
        """
        Method to setup the pair potentials based on the per atom type assigned parameters
        :Parameters:
            - radfact (int): factor to be multiplied during radius generation, default to 1.0
            - radrule (str): radiusrule, default to arithmetic
            - epsrule (str): epsilonrule, default to geometric
        """
        self.vdwdata = {}
        self.types2numbers = {} #equivalent to self.dlp_types
        types = self.par["vdw"].keys()
        #print(types)
        for i, t in enumerate(types):
            if t not in self.types2numbers.keys():
                self.types2numbers[t]=str(i)
        ntypes = len(types)
        #print (ntypes)
        for i in range(ntypes):
            for j in range(i, ntypes):
                #TODO check availability of an explicit paramerter
                #if len(self.par["vdwpr"].keys()) > 0:
                #    _,_,ti = self.split_parname(types[i])
                #    _,_,tj = self.split_parname(types[j])
                #    ti =ti[0]
                #    tj = tj[0]
                #    parname = self.build_parname("vdwpr", "buck6d", "leg", [ti,tj])
                #    if parname in self.par["vdwpr"].keys():
                #        # got it
                #        par_ij = self.par["vdwpr"][parname]
                #        self.vdwdata[types[i]+":"+types[j]] = par_ij
                #        self.vdwdata[types[j]+":"+types[i]] = par_ij
                #        continue
                par_i = self.par["vdw"][types[i]][1]
                par_j = self.par["vdw"][types[j]][1]
                pot_i =  self.par["vdw"][types[i]][0]
                pot_j =  self.par["vdw"][types[j]][0]
                if pot_i == pot_j:
                    pot = pot_i
                else:
                    raise IOError("Can not combine %s and %s" % (pot_i, pot_j))
                if self.settings["radrule"] == "arithmetic":
                    rad = self.settings["radfact"]*(par_i[0]+par_j[0])
                elif self.settings["radrule"] == "geometric":
                    rad = 2.0*self.settings["radfact"]*np.sqrt(par_i[0]*par_j[0])
                else:
                    raise IOError("Unknown radius rule %s specified" % self.settings["radrule"])
                if self.settings["epsrule"] == "arithmetic":
                    eps = 0.5 * (par_i[1]+par_j[1])
                elif self.settings["epsrule"] == "geometric":
                    eps = np.sqrt(par_i[1]*par_j[1])
                else:
                    raise IOError("Unknown radius rule %s specified" % self.settings["radrule"])
                par_ij = (pot,[rad,eps])
                # all combinations are symmetric .. store pairs bith ways
                self.vdwdata[types[i]+":"+types[j]] = par_ij
                self.vdwdata[types[j]+":"+types[i]] = par_ij   
        #import pdb; pdb.set_trace()
        self.pair_potentials_initalized = True
        return


    @timer("find reference systems")
    def find_refsystems(self, plot=None):
        """
        function to detect the reference systems:
            - self.scan_ref      : list of ref names in the order to be searched
            - self.ref_systems   : dictionary of mol objects
            - self.ref_fraglist  : list of fragment indices belonging to this refsystem
            - self.ref_params    : paramtere dictionaries per refsystem (n-body/type)
        """
        self.timer.start("get reference systems")
        scan_ref  = []
        scan_prio = []
        if self._mol.mpi_rank == 0:
            ref_dic = self.api.list_FFrefs(self.par.FF)
        else:
            ref_dic = []
        if self._mol.mpi_size > 1:
            ref_dic = self._mol.mpi_comm.bcast(ref_dic, root=0)
        for refname in ref_dic.keys():
            prio, reffrags, active, upgrades, atfix = ref_dic[refname]
            if len(reffrags) > 0 and all(f in self.fragments.get_fragnames() for f in reffrags):
                scan_ref.append(refname)
                scan_prio.append(prio)
            # check for upgrades
            elif upgrades and len(reffrags) > 0:
                oreffrags = copy.deepcopy(reffrags)
                for d,u in upgrades.items():
                    reffrags = [i.replace(d,u) for i in reffrags]
                    if all(f in self.fragments.get_fragnames() for f in reffrags):
                        scan_ref.append(refname)
                        scan_prio.append(prio)
        # sort to be scanned referecnce systems by their prio
        self.scan_ref = [scan_ref[i] for i in np.argsort(scan_prio)]
        self.scan_ref.reverse()
        self.timer.stop()
        # now get the refsystems and make their fraggraphs and atomistic graphs of their active space
        self.timer.start("make ref frag graphs")
        self.ref_systems = {}
        if self._mol.mpi_rank == 0:
            ref_mol_strs  = self.api.get_FFrefs_graph(self.scan_ref, out ="str")
        else:
            ref_mol_strs = {}
        if self._mol.mpi_size > 1:
            ref_mol_strs = self._mol.mpi_comm.bcast(ref_mol_str, root = 0)
        for ref in self.scan_ref:
#            if self._mol.mpi_rank == 0:
#                ref_mol_str = self.api.get_FFref_graph(ref, out="str")
#            else:
#                ref_mol_str = None
#            if self._mol.mpi_size > 1:
#                ref_mol_str = self._mol.mpi_comm.bcast(ref_mol_str, root=0)
            ref_mol = molsys.mol.from_string(ref_mol_strs[ref])
            ref_mol.addon("fragments")
            ref_mol.fragments.make_frag_graph()
            if plot:
                ref_mol.fragments.plot_frag_graph(ref, ptype="png", size=600, vsize=20, fsize=20)
            # if active space is defined create atomistic graph of active zone
            active = ref_dic[ref][2]
            if active: ref_mol.graph.make_graph(active)
            self.ref_systems[ref] = ref_mol
        self.timer.stop()
        # now search in the fraggraph for the reference systems
        self.timer.start("scan for ref systems")
        logger.info("Searching for reference systems:")
        self.ref_fraglists = {}
        self.ref_atomlists = {}
        for ref in copy.copy(self.scan_ref):
            # TODO: if a ref system has only one fragment we do not need to do a substructure search but
            #       could pick it from self.fragemnts.fraglist
            subs = self._mol.graph.find_subgraph(self.fragments.frag_graph, self.ref_systems[ref].fragments.frag_graph)
            # in the case that an upgrade for a reference system is available, it has also to be searched
            # for the upgraded reference systems
            upgrades = ref_dic[ref][3]
            if upgrades:
                # if upgrades should be applied, also an active zone has to be present
                assert ref_dic[ref][2] != None
                for s,r in upgrades.items():
                    self.ref_systems[ref].fragments.upgrade(s, r)
                    subs += self._mol.graph.find_subgraph(self.fragments.frag_graph, self.ref_systems[ref].fragments.frag_graph)
            logger.info("   -> found %5d occurences of reference system %s" % (len(subs), ref))
            if len(subs) == 0:
                # this ref system does not appear => discard
                self.scan_ref.remove(ref)
                del(self.ref_systems[ref])
            else:
                # join all fragments
                subs_flat = list(set(itertools.chain.from_iterable(subs)))
                self.ref_fraglists[ref] = subs_flat
                # now we have to search for the active space
                # first construct the atomistic graph for the sub in the real system if 
                # an active zone is defined
                if ref_dic[ref][2] != None and type(ref_dic[ref][2]) != str:
                    idx = self.fragments.frags2atoms(subs_flat)
                    self._mol.graph.filter_graph(idx)
                    asubs = self._mol.graph.find_subgraph(self._mol.graph.molg, self.ref_systems[ref].graph.molg)
                    ### check for atfixes and change atype accordingly, the atfix number has to be referred to its index in the azone
                    if ref_dic[ref][4] != None and type(ref_dic[ref][4]) != str:
                        atfix = ref_dic[ref][4]
                        for s in asubs:
                            for idx, at in atfix.items():
                                azone = ref_dic[ref][2]
                                self.aftypes[s[azone.index(int(idx))]].atype = at
                    self._mol.graph.molg.clear_filters()
                    asubs_flat = itertools.chain.from_iterable(asubs)
                    self.ref_atomlists[ref] = list(set(asubs_flat))
                else:
                    self.ref_atomlists[ref] = None
        self.timer.stop()
        # get the parameters
        self.timer.start("get ref parmeter sets")
        self.ref_params = {}
        for ref in self.scan_ref:
            logger.info("Getting params for %s" % ref)
            if self._mol.mpi_rank == 0:
                ref_par = self.api.get_params_from_ref(self.par.FF, ref)
            else:
                ref_par = None
            if self._mol.mpi_size > 1:
                ref_par = self._mol.mpi_comm.bcast(ref_par, root=0)                
            self.ref_params[ref] = ref_par
            #print(("DEBUG DEBUG Ref system %s" % ref))
            #print((self.ref_params[ref]))
        self.timer.stop()
        return


    def atoms_in_subsys(self, alist, fsubsys):
        """
        this helper function checks if all fragments of atoms (indices) in alist
        appear in the list of fragments (indices) in fsubsys
        :Parameters:
            - alist(list): list of atom indices
            - fsubsys(list): list of fragment indices 
        :Returns:
            - True or False
        """
        return all(f in fsubsys for f in map(lambda a: self._mol.fragnumbers[a], alist))

    def atoms_in_active(self, alist, subsys):
        """
        this helper function checks if any atom (indices) in alist
        appear in the list of active atoms (indices) in fsubsys

        Args:
            alist (list): list of atom indices
            subsys (list): list of atom indices
            
        Returns:
            bool: True or False
        """
        if subsys is None: return True
        return any(a in subsys for a in alist)

    def get_parname(self, alist):
        """
        helper function to produce the name string using the self.aftypes
        """
        l = map(lambda a: self.aftypes[a], alist)
        return tuple(l)

    def get_parname_equiv(self,alist, ic, refsys):
        """
        helper function to produce the name string using the self.aftypes
        and the information in the self.equivs dictionary
        :Parameters:
            - alist(list): list of atom indices
            - ic(str): corresponding internal coordinate
            - refsys(str): corresponding name of reference system
        :Returns:
            - parname(str): parname
        """
        assert type(ic) == type(refsys) == str
        ### first check if an atom in r is in the predifined active zone
        insides = []
        for i in alist:
            if self.active_zone.count(i) > 0: insides.append(i)
        ### now perform the actual lookup
        try:
            # check for equivs
            equivs = self.equivs[refsys][ic]
            # ok, got some --> try to apply
            # first check if for all insides an equiv is available
            for i in insides: 
                if i not in equivs.keys(): return None
            if len(insides) > 1: return None
            return map(lambda a: self.aftypes[a] if a not in equivs.keys() else equivs[a], alist)
        except:
            if len(insides) > 0: return None
            return map(lambda a: self.aftypes[a], alist)

    def get_parname_sort(self, alist, ic):
        """
        helper function to produce the name string using the self.aftypes
        """
        l = map(lambda a: self.aftypes[a], alist)
        return tuple(aftype_sort(l,ic))

    def split_parname(self,name):
        """
        Helper function to exploit the necessary information from the parname:
        :Parameters:
            - name(str): parname
        :Returns:
            - pot(str): potential type
            - ref(str): name of the reference system
            - aftypes(list): list of aftype names
        """
        pot = name.split("-")[0]
        ref = name.split("|")[-1]
        aftypes = name.split("(")[1].split(")")[0].split(",")
        return pot, ref, aftypes

    def build_parname(self, ic, pot, ref, aftypes):
        """
        Helper function to build the parname out of the necessary information
        :Parameters:
            - ic(str): name of the internal coordinate
            - pot(str): pot name
            - ref(str): name of the reference system
            - aftypes(list): list of the involved aftypes
        :Returns:
            - parname
        """
        sorted = aftype_sort(aftypes, ic)
        return pot + "->("+string.join(sorted, ",")+")|"+ref
        

    def pick_params(self,aft_list,ic,at_list, pardir):
        """
        Hhelper function to pick params from the dictionary pardir using permutations for the given ic
        if len of aft_list == 1 (ic = vdw or cha) no permutations necessary
        :Parameters:
            -aft_list(list): list of aftypes
            -ic(str): name of internal coordinate
            -pardir(dict): dictionary holding the parameters
        :Returns:
            -parname(str): parname
            -params(list): list of the corresponding parameters
        """
        if aft_list is None: return (), None
        ic_perm = {"bnd": ((0,1), (1,0)),
                   "bnd5": ((0,1), (1,0)),
                   "ang": ((0,1,2), (2,1,0)),
                   "ang5": ((0,1,2), (2,1,0)),
                   "dih": ((0,1,2,3),(3,2,1,0)),
                   "oop": ((0,1,2,3),(0,1,3,2),(0,2,1,3),(0,2,3,1),(0,3,1,2),(0,3,2,1))}
        if len(aft_list) == 1:
            parname = tuple(aft_list)
            if pardir[ic].index(parname, "wild_ft") >= 0:
                return parname, pardir[ic].__getitem__(parname, wildcards = True)
            else:
                return (), None
        else:         
            # we have to check here also for ring specific parameters
            if ic == "bnd" and at_list.ring == 5:
                ics = ["bnd5", "bnd"]
            elif ic == "ang" and at_list.ring == 5:
                ics = ["ang5", "ang"]
            else:
                ics = [ic]
            # now we have to loop over ics
            for cic in ics:
                perm = ic_perm[cic]
                for p in perm:
                    parname = tuple(map(aft_list.__getitem__, p))
                    if pardir[cic].index(parname, "wild_ft") >= 0:
                        # if we found a bnd5 or ang5 we have to modify
                        # the name of the parameter
                        if cic == "bnd5" or cic == "ang5":
                            param = copy.deepcopy(pardir[cic].__getitem__(parname, wildcards = True))
                            lparname = list(parname)
                            lparname.append("r5")
                            return tuple(lparname), param
                        else:
                            return parname, pardir[cic].__getitem__(parname, wildcards = True)
            # if we get to this point all permutations gave no result
            return (), None

    def enumerate_types(self):
        """
        Helper function to assign a number to the type, needed for fileIO
        """
        # dummy dicts to assign a number to the type
        par_types = {}
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            ptyp = {}
            i = 1
            for ind in self.par[ic]:
                # cut off the potential type --- if the rest is the same we use the same number
                rind = ind.split("->")[1]
                if not rind in ptyp: 
                    ptyp[rind] = i
                    i += 1
            par_types[ic] = ptyp
        return par_types
 

    ################# IO methods #################################################

    def write(self, fname):
        """
        write the rics including the referencing types to an ascii file
        called <fname>.ric and the parameters to <fname>.par
        :Parameters:
            - fname(str): fname
        """
        if self.mpi_rank > 0:
            return
        # create a random hash to ensure that ric and par file belong together
        hash = str(uuid.uuid4())
        # dummy dicts to assign a number to the type
        par_types = self.enumerate_types()
        # write the RICs first
        f = open(fname+".ric", "w")
        f.write("HASH: %s\n" % hash)
        logger.info("Writing RIC to file %s.ric" % fname)
        # should we add a name here in the file? the FF goes to par. keep it simple ...
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            filt = None
            if ic == "dih":
                filt = ["ring"]
            ric = self.ric_type[ic]
            parind = self.parind[ic]
            ptyp = par_types[ic]
            f.write("%s %d\n" % (ic, len(ric)))
            for i,r in enumerate(ric):
                # we take only the first index and remove the ptype to lookup in ptyp dictionary
                pi = parind[i][0]
                ipi = ptyp[pi.split("->")[1]]
                f.write("%d %d %s\n" % (i+1, ipi, r.to_string(filt=filt, inc=1)))
            f.write("\n")
        f.close()
        # write the par file
        #if self.refsysname:
        if hasattr(self.par, 'variables'):
            # this is a fixed up refsystem for fitting
            f = open(fname+".fpar", "w")
            vals = self.par.variables.vals
            self.varnames2par()
        else:
            f = open(fname+".par", "w")             
            logger.info("Writing parameter to file %s.par" % fname)
        f.write("HASH: %s\n" % hash)
        f.write("FF %s\n\n" % self.par.FF)
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            ptyp = par_types[ic]
            par = self.par[ic]
            f.write(ff_desc[ic])
            f.write("\n")
            f.write("%3s_type %d\n" % (ic, len(par)))
            ind = par.keys()
            ind.sort(key=lambda k: ptyp[k.split("->")[1]])
            for i in ind:
                ipi = ptyp[i.split("->")[1]]
                ptype, values = par[i]
                formatstr = string.join(map(lambda a: "%15.8f" if type(a) != str else "%+15s", values))
                sval = formatstr % tuple(values)
                #sval = (len(values)*"%15.8f ") % tuple(values)
                f.write("%-5d %20s %s           # %s\n" % (ipi, ptype, sval, i))
            f.write("\n")
        if hasattr(self.par, 'variables'):
            self.par.variables(vals)
            if hasattr(self, 'active_zone'):
                active_zone = np.array(self.active_zone)+1
                f.write(("azone "+len(active_zone)*" %d"+"\n\n") % tuple(active_zone))
            if hasattr(self, 'refsysname'): f.write("refsysname %s\n\n" % self.refsysname)
            f.write("variables %d\n" % len(self.par.variables))
            for k,v in self.par.variables.items():
                f.write("%10s %15.8f %15.8f %15.8f %3s %3s\n" % (v.name, v.val, v.range[0], v.range[1], v.bounds[0], v.bounds[1]))
        f.close()
        return

    def write_key(self,fname,atype_map=False,atype_addendum='',write_default_header=True):
        '''
        author: Julian
        try to write a key file from the data available in the class
        needs then to be merged manually in case there is only a partial system
        #### !!!! EQUIVALENCES ARE MISSING !!!! ####
        '''
        a = atype_addendum
        fkey=open(fname, 'w')
        par = self.par
        if write_default_header == True:
            fkey.write(
'''
version        2.0
parameters     none


bondunit       71.94
angleunit      0.02191418
strbndunit     2.51118
opbendunit     0.02191418
torsionunit    0.5
vdwtype        exp6_damped
vdwdampfact    0.25
radiusrule     arithmetic
radiustype     r-min
radiussize     radius
epsilonrule    geometric
a-expterm      184000.0
b-expterm      12.0
c-expterm      2.25

bondtype       mixmorse_bde
strbndtype     mmff
opbendtype     mmff
chargetype     gaussian\n\n''')
            # still in the indentation of write_default_header
            # we need also atom definitions, get unique atype list from 'cha' dictionary
            for i,k in enumerate(par['cha'].keys()):
                atype = k.split('(')[-1].split(')')[0]
                elem = atype.split('_')[0][0:-1] ## convention: all atypes have .lt. 10 connections! 
                                                  # remove only last digit from emenent string
                fkey.write('atom  %s   %s\n' % (atype, elem))
            import pdb; pdb.set_trace()
        #syntax of keyfile:
        # red name [atypes] [params] 
        parkeys= self.par.keys() # cha ang dih oop vdw bnd
        atypes_set = []        
        # bonds 
        for bond in par['bnd'].keys():
            atype1,atype2 = bond.split('(')[-1].split(')')[0].split(',')
            atype1,atype2 = atype1+a,atype2+a
            partype, vals = par['bnd'][bond]
            if atype_map != False:
                #tbi
                pass
            fkey.write('%15s     %15s   %15s   %18.10f   %18.10f\n' % ('bond',atype1,atype2,vals[0],vals[1]))

        # angles
        for angle in par['ang'].keys():
            atype1,atype2,atype3 = angle.split('(')[-1].split(')')[0].split(',')
            atype1,atype2,atype3 = atype1+a,atype2+a,atype3+a
            partype,vals = par['ang'][angle]
            if partype == 'mm3':
                fkey.write('%15s     %15s   %15s   %15s  %18.10f   %18.10f\n' % ('angle',atype1,atype2,atype3,vals[0],vals[1]))
            elif partype == 'fourier':
                fkey.write('%15s     %15s   %15s   %15s  %18.10f   %18.10f %18.10f %18.10f %18.10f\n' % ('fourier',atype1,atype2,atype3,vals[0],vals[1],vals[2],vals[3],vals[4]))
            elif partype == 'strbnd':
                fkey.write('%15s     %15s   %15s   %15s  %18.10f   %18.10f %18.10f\n' % ('strbnd',atype1,atype2,atype3,vals[0],vals[1],vals[2]))
            else:
                raise IOError('partype %s not yet implemented' % partype)
        
        # torsions
        for torsion in par['dih'].keys():
            atype1,atype2,atype3,atype4 = torsion.split('(')[-1].split(')')[0].split(',')
            atype1,atype2,atype3,atype4 = atype1+a,atype2+a,atype3+a,atype4+a
            partype,vals = par['dih'][torsion]
            if partype == 'cos3':
                fkey.write('%15s     %15s   %15s   %15s   %15s   %18.10f   %18.10f   %18.10f\n' % ('torsion',atype1,atype2,atype3,atype4,vals[0],vals[1],vals[2]))
            elif partype == 'cos4':
                fkey.write('%15s     %15s   %15s   %15s   %15s   %18.10f   %18.10f   %18.10f   %18.10f\n' % ('torsion',atype1,atype2,atype3,atype4,vals[0],vals[1],vals[2],vals[3]))
            else:
                raise IOError('partype %s not yet implemented' % partype)
        
        # oops
        for oop in par['oop'].keys():
            atype1,atype2,atype3,atype4 = oop.split('(')[-1].split(')')[0].split(',')
            atype1,atype2,atype3,atype4 = atype1+a,atype2+a,atype3+a,atype4+a
            partype,vals = par['oop'][oop]
            if partype == 'harm':
                fkey.write('%15s     %15s   %15s   %15s   %15s   %18.10f   %18.10f\n' % ('opbend',atype1,atype2,atype3,atype4,vals[0],vals[1]))
            else:
                raise IOError('partype %s not yet implemented' % partype)

        for vdw in par['vdw'].keys():
            atype1 = vdw.split('(')[-1].split(')')[0]
            atype1 = atype1+a
            partype, vals = par['vdw'][vdw]
            if atype_map != False:
                #tbi
                pass
            fkey.write('%15s     %15s   %18.10f   %18.10f\n' % ('vdw',atype1,vals[0],vals[1]))
        
        for charge in par['cha'].keys():
            atype1 = charge.split('(')[-1].split(')')[0]
            atype1 = atype1+a
            partype, vals = par['cha'][charge]
            if atype_map != False:
                #tbi
                pass
            fkey.write('%15s     %15s   %18.10f   %18.10f\n' % ('charge',atype1,vals[0],vals[1]))
        fkey.close() 
        
        return



    def read(self, fname, fit=False):
        """
        read the ric/par files instead of assigning params
        :Parameters:
            - fname(str) : name of <fname>.par und <fname.ric> file
            - fit(bool, optional): specify if an fpar file should be read in, 
            holding fitting information, defaults to False
        """
        hash = None
        fric = open(fname+".ric", "r")
        ric_type = ["bnd", "ang", "dih", "oop", "cha", "vdw"]
        ric_len  = [2    , 3    , 4    , 4    , 1    , 1    ]
        ric      = {}
        # read in ric first, store the type as an attribute in the first place
        stop = False
        assigned = []
        while not stop:
            line = fric.readline()
            if len(line)==0:
                # end of ric file
                stop = True
            sline = line.split()
            if len(sline)> 0:
                if sline[0] == "HASH:": 
                    hash = sline[1]
                elif sline[0] in ric_type:
                    curric = sline[0]
                    curric_len = ric_len[ric_type.index(curric)]
                    assigned.append(curric)
                    nric = int(sline[1])
                    rlist = []
                    for i in range(nric):
                        sline = fric.readline().split()
                        rtype = int(sline[1])
                        aind  = map(int, sline[2:curric_len+2])
                        aind  = np.array(aind)-1
                        icl = ic(aind, type=rtype)
                        for attr in sline[curric_len+2:]:
                            atn,atv = attr.split("=")
                            icl.__setattr__(atn, int(atv))
                        rlist.append(icl)
                    ric[curric] = rlist    
        fric.close()
        logger.info("read RIC from file %s.ric" % fname)
        # now add data to ric object .. it gets only bnd, angl, oop, dih
        self.ric.set_rics(ric["bnd"], ric["ang"], ric["oop"], ric["dih"])
        # time to init the data structures .. supply vdw and cha here
        self._init_data(cha=ric["cha"], vdw=ric["vdw"])
        self._init_pardata()
        # now open and read in the par file
        if fit:
            nkeys={}
            self.fit=True
            fpar = open(fname+".fpar", "r")
            # in the fit case we first screen for the variables block and read it ini
            self.par.attach_variables()
            #self.variables = varpars()
            line = fpar.readline()
            stop = False
            azone = False
            vars  = False
            refsysname = False
            while not stop:
                sline = line.split()
                if len(sline)>0:
                    if sline[0] == "azone":
                        self.active_zone = (np.array(map(int, sline[1:]))-1).tolist()
                        azone = True
                    elif sline[0] == "variables":
                        nvar = int(sline[1])
                        for i in range(nvar):
                            sline = fpar.readline().split()
#                            nkey = self.par.variables[sline[0]] = varpar(self.par, sline[0], 
#                                          val = float(sline[1]), 
#                                          range = [float(sline[2]), float(sline[3])], 
#                                          bounds = [sline[4], sline[5]])
                            nkey = self.par.variables.__setitem__(sline[0], varpar(self.par, sline[0], 
                                          val = float(sline[1]), 
                                          range = [float(sline[2]), float(sline[3])], 
                                          bounds = [sline[4], sline[5]]))
                            if nkey != sline[0]:nkeys[sline[0]]=nkey
                        vars = True
                    elif sline[0] == "refsysname":
                        self.refsysname = sline[1]
                        refsysname = True
#                    if azone == vars == refsysname == True:
                    if vars == True:
                        fpar.seek(0)
                        stop = True
                        break
                line = fpar.readline()
                if len(line) == 0:
                    raise IOError("Variables block and/or azone in fpar is missing!")
        else:
            fpar = open(fname+".par", "r")
        stop = False
        found_hash = False
        while not stop:         
            line = fpar.readline()
            if len(line) == 0:
                stop = True
            sline = line.split()
            if len(sline)>0:
                if sline[0] == "HASH:":
                    found_hash = True
                    assert sline[1] == hash, "Hashes of ric and par file do not match!"
                if sline[0][0] == "#": continue 
                curric = sline[0].split("_")[0]
                if sline[0]=="FF":
                    self.par.FF = sline[1]
                elif curric in ric_type:
                    par = self.par[curric]
                    t2ident = {} # maps integer type to identifier
                    ntypes = int(sline[1])
                    for i in range(ntypes):
                        sline = fpar.readline().split()
                        if sline[0][0] == "#": continue 
                        # now parse the line 
                        itype = int(sline[0])
                        ptype = sline[1]
                        ident = sline[-1]
                        param = sline[2:-2]
                        if self.fit:
                            newparam = []
                            # if we read a fpar file we need to test if there are variables
                            for paridx,p in enumerate(param):
                                if p[0] == "$":
                                    # check if variable name was overwritten
                                    if p in nkeys: p = nkeys[p]
                                    #import pdb; pdb.set_trace()
                                    if not p in self.par.variables:
                                        raise IOError("Varible %s undefiend in variable block" % p)
                                    # for multistruc fits the $ name are not anymore valid, it has
                                    # to be checked firtst if the variable is already defined under
                                    # a different name
                                    found = False
                                    for vname in self.par.variables.keys():
                                        if (curric, ident, paridx) in self.par.variables[vname].pos:
                                            found = True
                                            p = vname
                                            break
                                    if not found:
                                        self.par.variables[p].pos.append((curric,ident,paridx))
                                    newparam.append(p)
                                else:
                                    newparam.append(float(p))
                            param = newparam
                        else:
                            param = map(float, param)
                        if ident in par:
                            logger.warning('Identifier %s already in par dictionary --> will be overwritten' % ident)
                            raise ValueError("Identifier %s appears twice" % ident)
                        par[ident] = (ptype, param)
                        if itype in t2ident:
                            t2ident[itype].append(ident)
                        else:
                            t2ident[itype] = [ident]
                    # now all types are read in: set up the parind datastructure of the ric
                    parind = self.parind[curric]
                    for i,r in enumerate(self.ric_type[curric]):
                        parind[i] = t2ident[r.type]
        fpar.close()
        # check if both fpar and ric file have hashes
        if hash is not None and found_hash == False:
            raise IOError("ric file has a hash, but par has not")
        logger.info("read parameter from file %s.par" % fname)
        ### replace variable names by the current value
        if self.fit: 
            self.par.variables.cleanup()
            self.par.variables()
        return

    def load_params_from_parfile(self, fname, fit = True):
        if fit:
            fname = '%s.fpar' % fname
        else:
            fname = '%s.par' % fname
        self._init_pardata()
        #ric_type = self.ric_type
        if fit:
            with open(fname, 'r') as fpar:
                #self.fit=True
                # in the fit case we first screen for the variables block and read it in
                self.par.attach_variables()
                #self.variables = varpars()
                line = fpar.readline()
                stop = False
                azone = False
                vars  = False
                refsysname = False
                while not stop:
                    sline = line.split()
                    if len(sline)>0:
                        #if sline[0] == "azone":
                        #    self.active_zone = (np.array(map(int, sline[1:]))-1).tolist()
                        #    azone = True
                        if sline[0] == "variables":
                            nvar = int(sline[1])
                            for i in range(nvar):
                                sline = fpar.readline().split()
                                self.par.variables[sline[0]] = varpar(self.par, sline[0], 
                                              val = float(sline[1]), 
                                              range = [float(sline[2]), float(sline[3])], 
                                              bounds = [sline[4], sline[5]])
                            vars = True
                        #elif sline[0] == "refsysname":
                        #    self.refsysname = sline[1]
                        #    refsysname = True
                        #if azone == vars == refsysname == True:
                        if vars == True:
                            stop = True
                            break
                    line = fpar.readline()
                    if len(line) == 0:
                        raise IOError("Variables block in fpar is missing!")
        with open(fname, 'r') as fpar:
            stop = False
            while not stop:
                line = fpar.readline()
                if len(line) == 0:
                    stop = True
                sline = line.split()
                if len(sline)>0:
                    if sline[0][0] == "#": continue 
                    curric = sline[0].split("_")[0]
                    if sline[0]=="FF":
                        self.par.FF = sline[1]
                    elif curric in ['bnd','ang','dih','oop','cha','vdw']:
                        par = self.par[curric]
#                        t2ident = {} # maps integer type to identifier
                        ntypes = int(sline[1]) # number of potentials for current ric
                        for i in range(ntypes):
                            sline = fpar.readline().split()
                            if sline[0][0] == "#": continue 
                            # now parse the line 
                            itype = int(sline[0])
                            ptype = sline[1]
                            ident = sline[-1]
#                            pot = ident.split('-')
#                            if pot not in self.loaded_pots[curric]: self.loaded_pots[curric].append(pot[0])
                            param = sline[2:-2]
                            if fit:
                                newparam = []
                                # if we read a fpar file we need to test if there are variables
                                for paridx,p in enumerate(param):
                                    if p[0] == "$":
                                        if not p in self.par.variables:
                                            raise IOError("Varible %s undefiend in variable block" % p)
                                        self.par.variables[p].pos.append((curric,ident,paridx))
                                        newparam.append(p)
                                    else:
                                        newparam.append(float(p))
                                param = newparam
                            else:
                                param = map(float, param)
                            #if ident in par: raise ValueError("Identifier %s appears twice" % ident)
                            if ident in par:
                                logger.warning('Identifier %s already in par dictionary --> will be overwritten' % ident)
                            par[ident] = (ptype, param)
#                            if itype in t2ident:
#                                t2ident[itype].append(ident)
#                            else:
#                                t2ident[itype] = [ident]
                        # now all types are read in: set up the parind datastructure of the ric
#                        parind = self.parind[curric]
#                        for i,r in enumerate(self.ric_type[curric]):
#                            parind[i] = t2ident[r.type]
        if fit: 
            self.par.variables.cleanup()
            self.par.variables()
        logger.info("read parameter from file %s" % fname)
        return

    def to_latex(self, refsys, ic, datypes = None):
        """
        Method to publish parameters as latex tables. This method needs pandas installed.

        Args:
            refsys (str): reference system for which the parameters should
                be published.
            ic (str): Internal coordinate for which the parameters should
                be published
            datypes (dict, optional): Defaults to None. Dictionary to 
                convert atomtypes to more readable types.
        
        Returns:
            str: latex code as str
        """
   
        import pandas as pd
        mdyn2kcal=143.88
        units = {
            "cha": {},
            "vdw": {},
            "bnd": {"$k_b [\si{\kcal\per\mole\per\AA\squared}]$":mdyn2kcal},
            "bnd5": {"$k_b [\si{\kcal\per\mole\per\AA\squared}]$":mdyn2kcal},
            "ang": {"$k_a [\si{\kcal\per\mole\per\radian\squared}]$":mdyn2kcal},
            "ang5": {"$k_a [\si{\kcal\per\mole\per\radian\squared}]$"},
            "dih": {},
            "oop": {"$k_o [\si{\kcal\per\mole\per\radian\squared}]$":mdyn2kcal},
        }
        potentials = {"cha": {"gaussian": ["$q [e]$", "$sigma [\si{\AA}]$"]},
                "vdw": {"buck6d": ["\text{$r^0_{ij} [\si{\AA}]$}", "$\epsilon_{ij}$"]},
                "bnd": {"mm3": ["$k_b [\si{\kcal\per\mole\per\AA\squared}]$", "\text{$r_b^0 [\si{\AA}]$}"]},
                "bnd5": {"mm3": ["$k_b [\si{\kcal\per\mole\per\AA\squared}]$", "\text{$r_b^0 [\si{\AA}]$}"]},
                "ang": {"mm3": ["$k_a [\si{\kcal\per\mole\per\radian\squared}]$", "\text{$\phi_a^0 [\si{\degree}]$} "]},
                "ang5": {"mm3": ["$k_a [\si{\kcal\per\mole\per\radian\squared}]$", "\text{$\phi_a^0 [\si{\degree}]$} "]},
                "dih": {"cos3": ["\text{$V_1 [\si{\kcal\per\mole}]$}", "\text{$V_2 [\si{\kcal\per\mole}]$}", "\text{$V_3 [\si{\kcal\per\mole}]$}",], 
                "cos3": ["\text{$V_1 [\si{\kcal\per\mole}]$}", "\text{$V_2 [\si{\kcal\per\mole}]$}", "\text{$V_3 [\si{\kcal\per\mole}]$}","\text{$V_4 [\si{\kcal\per\mole}]$}"]},
                "oop": {"harm": ["$k_o [\si{\kcal\per\mole\per\radian\squared}]$", "\text{$\gamma^0$ [\si{\degree}]}"]}}
        sortlist = ["types", "potential",]
        columnformat = "lcc"
        d = {"types":[], "potential":[]}
        for i,l in potentials[ic].items():
            for j in l: 
                d[j]=[]
                sortlist.append(j)
                columnformat+="S"
        for k,v in self.par[ic].items():
            pot, ref, atypes = self.split_parname(k)
            if pot not in potentials[ic].keys():continue
            if ref == refsys:
                if datypes is not None: atypes = map(lambda a: datypes[a],atypes)
                d["types"].append(string.join(atypes, "-"))
                d["potential"].append(pot)
                for i, p in enumerate(v[1]):
                    ptype = potentials[ic][pot][i]
                    if ptype in units[ic].keys(): p *= units[ic][ptype]
                    d[potentials[ic][pot][i]].append(round(p,3))
        # remove empty columns
        rkeys = []
        for k,v in d.items():
            if len(v)==0: rkeys.append(k)
        for k in rkeys: 
            del d[k]
            sortlist.remove(k)
        df = pd.DataFrame(data=d)
        df = df[sortlist]
        return df.to_latex(escape = False, column_format=columnformat)

    def assign_params_from_key(self, fname):
        """
        Method used to setup an ff addon based on the legacy key file format.
        
        Args:
            fname (str): filename of the key file
        """


        def numberify_list(l, formatcode):
            if len(formatcode) == 0:
                return l
            elif formatcode[0] == "*":
                if formatcode[1] == "i":
                    try:
                        l = map(string.atoi,l)
                    except ValueError:
                        print ("Error converting %s to integer!!" % str(l))
                        raise ValueError
                elif formatcode[1] == "f":
                    l = map(string.atoi,l)
                elif formatcode[1] == "s":
                    pass
                else:
                    raise ValueError, "unknown formatcode %s" % formatcode
            else:
                for i in xrange(len(formatcode)):
                    if formatcode[i] == "i":
                        try:
                            l[i] = string.atoi(l[i])
                        except ValueError:
                            print ("Error converting %s to integer!!" % str(l[i]))
                            raise ValueError
                    if formatcode[i] == "f":
                        try:
                            l[i] = string.atof(l[i])
                        except ValueError:
                            print ("Error converting %s to float!!" % str(l[i]))
                            raise ValueError
                        except IndexError:
                            print ("Formatcode: %s" % formatcode)
                            print ("Data      : %s" % str(l))
                            raise IndexError('Maybe is the last opbend parameter missing? (not sure)')
                    if formatcode[i] == "o":
                        # optional float .. need to check if l is long enough ... o must be at the end!!
                        if (len(l) > i) :
                            l[i] = string.atof(l[i])
                        else:
                            l.append(0.0)
            return l
        # idea is to read the stuff first into self.par dictionary and afterwards
        # assign them
        # available potential keywords, keyword not yet implemented here are 
        # commented out
        pterms = {
            #"atom"        : (1, "so") , \
            "vdw"         : (1, "ff", "vdw", "buck6d")    ,\
            "vdwpr"       : (2, "ff", "vdwpr", "buck6d")    ,\
            "bond"        : (2,"ffo", "bnd", "mm3")   ,\
            #"bondq"       : (2,"ffff")   ,\
            #"bond5"       : (2,"ffo")   ,\
            #"bond4"       : (2,"ff")   ,\
            "angle"       : (3, "ff", "ang", "mm3") ,\
            #"angleq"       : (3, "ffff") ,\
            #"angle5"      : (3, "ffoo") ,\
            #"angle4"      : (3, "ffoo") ,\
            "anglef"      : (3, "ffioo", "ang", "fourier")  ,\
            #"anglef-2"    : (3, "ff")   ,\
            "strbnd"      : (3, "fff", "ang", "strbnd") ,\
            "opbend"      : (4, "fo", "oop", "harm")   ,\
            "torsion"     : (4, "fffo", "dih", "cos3") ,\
            #"torsion5"    : (4, "fffo") ,\
            #"torsion4"    : (4, "fffo") ,\
            "charge"      : (1, "ff", "cha", "gaussian") ,\
            #"chargemod"   : (2, "f") ,\
            #"chargeadd"   : (1, "f") ,\
            #"molname"     : (1, "*s") ,\
            #"virtbond"    : (2, "") ,\
            #"restrain-distance"  : (2, "ff"),\
            #"restrain-angle"     : (3, "ff"),\
            #"restrain-torsion"   : (4, "ff"),\
            }
        not_implemented = ["equivalent", "chargemod", "bond4", "bond5", 
            "angle5", "angle4", "restrain-distance","restrain-angle",
            "restrain-angle",]
        self._init_pardata()
        ptnames = pterms.keys()
        with open(fname, "r") as f:
            for line in f.readlines():
                sline = string.split(line)
                # jump over comments and empty lines
                if ((len(sline)>0) and (sline[0][0] != "#")):
                    keyword = sline[0]
                else:
                    continue
                # check for notes
                if sline.count("!!"):
                    sline = sline[:sline.index("!!")]
                # check if keyword is a potential keyword then read it in
                # ignore setttings at the first place
                if keyword in not_implemented:
                    raise NotImplementedError("%s not supported" % keyword)
                if ptnames.count(keyword):
                    natoms, formatcode, ic, pot = pterms[keyword]
                    # we have to get the atomkey and sort it
                    atomkey = aftype_sort([aftype(i,"leg") for i in sline[1:natoms+1]],ic)
                    # get actual params from numberify_list, stolen from pydlpoly.ff module
                    # we need to perform here some checks in order to distinct for example
                    # mm3 bnd from morse bnd potential
                    params = numberify_list(sline[natoms+1:], formatcode)
                    # we need to create a parname we use "leg" for legacy
                    # as name of the refsys and as framentnames to make it
                    # easy to use already implemented helper methods
                    # for the actual assignment
                    full_parname = pot+"->("+string.join(map(str,atomkey),",")+")|leg"
                    if ic == "bnd" and params[2] > 0.0: pot = "morse"
                    if ic == "dih" and params[3] > 0.0: pot = "cos4"
                    self.par[ic][full_parname] = (pot, numberify_list(sline[natoms+1:], formatcode))
        # now we have to add the missing parameters for strbnd by looping over the bond and angle pots
        for k in self.par["ang"].keys():
            if self.par["ang"][k][0] == "strbnd":
                pot,ref, aftypes = self.split_parname(k)
                apot  = "mm3->"+k.split("->")[-1] 
                spot1 = self.build_parname("bnd", "mm3", "leg", aftypes[:2])
                spot2 = self.build_parname("bnd", "mm3", "leg", aftypes[1:])
                s1 = self.par["bnd"][spot1][1][1]
                s2 = self.par["bnd"][spot2][1][1]
                a  = self.par["ang"][apot][1][1]
                pot = ("strbnd", self.par["ang"][k][1]+[s1,s2,a])
                self.par["ang"][k]=pot                
        # now perform the actual assignment loop
        self.assign_params_offline(ref="leg", key = True)
        return




    
    def upload_params(self, FF, refname, dbrefname = None, azone = True, atfix = None, interactive = True):
        """
        Method to upload interactively the parameters to the already connected db.
        
        :Parameters:
            - FF(str): Name of the FF
            - refname(str): name of the refsystem for which params should be uploaded
            - dbrefname(str, optional): name of the refsystem in the db
            - azone(bool, optional): boolean flag indicating if an active zone entry in 
            the db should be created, defaults to True
            - atfix(list, optional): list of special atypes, which should be created 
            in the db for the reference system, defaults to None
            - interactive(bool, optional): specify if an interactive upload should be 
            done, defaults to True
        """
        assert type(refname) == str
        assert type(FF)      == str
        from mofplus import FF_api
        self.api = FF_api()
        if dbrefname is None: dbrefname = refname
        if atfix is not None:
            fixes = {}
            for i, at in enumerate(self._mol.atypes):
                if i in self.active_zone:
                    if at in atfix: fixes[str(i)]=at
        else: 
            fixes = None
        if azone:
            self.api.create_fit(FF, dbrefname, azone = self.active_zone, atfix = fixes)
        else:
            self.api.create_fit(FF, dbrefname)
        uploads = {
                "cha": {},
                "vdw": {}, 
                "bnd": {}, 
                "ang": {}, 
                "dih": {}, 
                "oop": {}}
        for ic,v in self.parind.items():
            for pl in v:
                for pn in pl: 
                    par = self.par[ic][pn]
                    pot, ref, aftypes = self.split_parname(pn)
                    if ref == refname:
                        if (tuple(aftypes), pot) not in uploads[ic].keys():
                            uploads[ic][(tuple(aftypes), pot)] = par[1]
        for ptype, upls in uploads.items():
            for desc, params in upls.items():
                # TODO: remove inconsitenz in db conserning charge and cha
                if ptype == "cha":
                    if interactive:
                        self.api.set_params_interactive(FF, desc[0], "charge", desc[1], dbrefname, params)
                    else:
                        self.api.set_params(FF, desc[0], "charge", desc[1], dbrefname, params)
                else:
                    if interactive:
                        self.api.set_params_interactive(FF, desc[0], ptype, desc[1], dbrefname, params)
                    else:
                        self.api.set_params(FF, desc[0], ptype, desc[1], dbrefname, params)
        return

