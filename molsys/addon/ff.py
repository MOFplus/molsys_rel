# -*- coding: utf-8 -*-
# RS .. to overload print in parallel case (needs to be the first line)
from __future__ import print_function
"""
Created on Thu Mar 23 11:25:43 2017

@author: rochus


        addon module FF to implement force field infrastructure to the molsys

"""

from mpi4py import MPI
mpi_rank = MPI.COMM_WORLD.Get_rank()
# overload print function in parallel case
import __builtin__
def print(*args, **kwargs):
    if mpi_rank == 0:
        return __builtin__.print(*args, **kwargs)
    else:
        return


import numpy as np
from molsys.util.timing import timer, Timer
from molsys.util import elems
from molsys.util.aftypes import aftype, aftype_sort,afdict

import itertools
import copy
import string
import json

import logging
import pdb
logger = logging.getLogger("molsys.ff")

class ic(list):

    def __init__(self, *args, **kwargs):
        list.__init__(self,*args)
        for k,v in kwargs.items():
            setattr(self, k, v)
        return


class ric:
    """
    class to detect and keep all the (redundant) internal coordinates of the system
    """


    def __init__(self, mol, specials = {"linear": []}):
        """
        find all bonds, angles, oops and dihedrals in the molsys mol object
        """

        self.timer = Timer()
        
        self._mol      = mol
        self.conn      = mol.get_conn()
        self.natoms    = mol.get_natoms()
        self.xyz       = mol.xyz
        self.aftypes   = []
        for i, a in enumerate(mol.get_atypes()):
            self.aftypes.append(aftype(a, mol.fragtypes[i]))

        self.bnd    = self.find_bonds()
        self.ang    = self.find_angles()
        self.oop    = self.find_oops()
        self.dih    = self.find_dihedrals(specials["linear"])
        self.report()
        self.timer.write_logger(logger.info)
        return


    # the follwoing code is adapted from pydlpoly/py/assign_FF.py
    # instead of bond objects we use simpler lists of lists with indices

    @timer("find bonds")
    def find_bonds(self):
        bonds=[]
        for a1 in xrange(self.natoms):
            for a2 in self.conn[a1]:
                if a2 > a1: bonds.append(ic([a1, a2]))
        return bonds

    @timer("find angles")
    def find_angles(self):
        angles=[]
        for ca in xrange(self.natoms):
            apex_atoms = self.conn[ca]
            naa = len(apex_atoms)
            for ia in xrange(naa):
                aa1 = apex_atoms[ia]
                other_apex_atoms = apex_atoms[ia+1:]
                for aa2 in other_apex_atoms:
                    angles.append(ic([aa1, ca, aa2]))
        return angles

    @timer("find oops")
    def find_oops(self):
        oops=[]
        # there are a lot of ways to find oops ...
        # we assume that only atoms with 3 partners can be an oop center
        for ta in xrange(self.natoms):
            if (len(self.conn[ta]) == 3):
                # ah! we have an oop
                a1, a2, a3 = tuple(self.conn[ta])
                oops.append(ic([ta, a1, a2, a3]))
                #oops.append([ta, a2, a1, a3])
                #oops.append([ta, a3, a1, a2])
        return oops

    @timer("find dihedrals")
    def find_dihedrals(self, lin_types = []):
        dihedrals=[]
        for a2 in xrange(self.natoms):
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
                    while self.aftypes[a2] in lin_types:
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
                    while self.aftypes[a3] in lin_types:
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
                            smallring = 0
                            if a1 == a4: continue
                            if con1.count(a4):
                                smallring = 4
                            else:
                                con4 = list(self.conn[a4])
                                for c1 in con1:
                                    if con4.count(c1):
                                        smallring = 5
                                        break
                            d = ic([a1,a2,a3,a4], smallring = smallring)
                            if lin:
                                ### in case of a dihedral due to dihedral shifts,
                                ### it has to be checked if we have this dihedral already
                                if d not in dihedrals:
                                    dihedrals.append(d)
                            else:
                                dihedrals.append(d)
        return dihedrals

    def get_distance(self,atoms):
        """
        Computes distance between two atoms
        :Parameters:
            - atoms (list): list of atomindices
        """
        xyz = self._mol.map2image(self.xyz[atoms])
        apex_1 = xyz[0]
        apex_2 = xyz[1]
        return np.linalg.norm(apex_1-apex_2)

    def get_angle(self,atoms):
        """
        Computes angle between three atoms
        :Parameters:
            - atoms (list): list of atomindices
        """
        xyz = self._mol.map2image(self.xyz[atoms])
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

    def get_dihedral(self, atoms):
        """
        Computes dihedral angle between four atoms
        :Parameters:
            - atoms (list): list of atomindices
        """
        xyz = self._mol.map2image(self.xyz[atoms])
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
        return phi * (180.0/np.pi)

    def get_oop(self,atoms):
        """
        Dummy function to the the value of an oop by default to 0.0
        :Parameters:
            - atoms (list): list of atomindices
        """
        return 0.0

    def compute_rics(self):
        """
        Computes the values of the rics and attaches 
        them to the corresponding ic
        """
        for b in self.bnd: b.value = self.get_distance(b)
        for a in self.ang: a.value = self.get_angle(a)
        for d in self.dih: d.value = self.get_dihedral(d)
        for o in self.oop: o.value = self.get_oop(o)
        return

    def report(self):
        logger.info("Reporting RICs")
        logger.info("%7d bonds"     % len(self.bnd))
        logger.info("%7d angles"    % len(self.ang))
        logger.info("%7d oops"      % len(self.oop))
        logger.info("%7d dihedrals" % len(self.dih))
        return



class ff:

    def __init__(self, mol,source = "mofp"):
        """
        instantiate a ff object which will be attached to the parent mol

        :Parameter:

             - mol : a mol type object (can be a derived type like bb or topo as well)
             - source: [string, default=mofp] where to get the data from
        """

        self.timer = Timer()
        self._mol = mol
        ### init api
        self.timer.start("connect to DB")
        self.source = source
        if self.source == "mofp":
            if self._mol.mpi_rank == 0:
                from mofplus import FF_api
                self.api = FF_api()
                self.special_atypes = self.api.list_special_atypes()
            else:
                self.api = None
                self.special_atypes = None
            if self._mol.mpi_size > 1:
                self.special_atypes = self._mol.mpi_comm.bcast(self.special_atypes, root = 0)
        elif self.source == "file":
            raise ValueError, "to be implemented"
        else:
            raise ValueError, "unknown source for assigning parameters"
        self.timer.stop()
        self.ric = ric(mol, specials = self.special_atypes)
        logger.debug("generated the ff addon")
        return

    @timer("assign parameter")
    def assign_params(self, FF):
        """
        method to orchestrate the parameter assignment for this system using a force field defined with
        FF

        :Parameter:

            - FF :    [string] name of the force field to be used in the parameter search
        """
        def get_equivalences():
            ### check for equivalences
            for aidx in r:
                aft = self.aftypes[aidx]
                if ((str(aft) == par[1][0]) and (aidx not in curr_equi_par)):
                    curr_equi_par[aidx] = par[1][1]
            return

        self.FF = FF
        # as a first step we need to generate the fragment graph
        self.timer.start("fragment graph")
        self._mol.addon("fragments")
        self.fragments = self._mol.fragments
        self.fragments.make_frag_graph()
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
        self.find_refsystems()
        # make data structures
        with self.timer("make data structures"):
            self.ric_type = {"cha":[[i] for i in range(self._mol.natoms)],
                    "vdw":[[i] for i in range(self._mol.natoms)], 
                    "bnd":self.ric.bnd, 
                    "ang":self.ric.ang, 
                    "dih":self.ric.dih, 
                    "oop":self.ric.oop}
            self.parind = {
                "bnd": [None]*len(self.ric.bnd),
                "ang": [None]*len(self.ric.ang),
                "dih": [None]*len(self.ric.dih),
                "oop": [None]*len(self.ric.oop),
                "cha": [None]*self._mol.natoms,
                "vdw": [None]*self._mol.natoms,
                }
            self.par = {
                "bnd": {},
                "ang": {},
                "dih": {},
                "oop": {},
                "vdw": {},
                "cha": {},
                }
        with self.timer("parameter assignement loop"):
            for ref in self.scan_ref:
                counter = 0
                logger.info("assigning params for ref system %s" % ref)
                curr_fraglist = self.ref_fraglists[ref]
                curr_atomlist = self.ref_atomlists[ref]
                curr_par = {\
                    "bnd" : self.ref_params[ref]["twobody"]["bnd"],\
                    "ang" : self.ref_params[ref]["threebody"]["ang"],\
                    "dih" : self.ref_params[ref]["fourbody"]["dih"],\
                    "oop" : self.ref_params[ref]["fourbody"]["oop"],
                    "cha" : self.ref_params[ref]["onebody"]["charge"],
                    "vdw" : self.ref_params[ref]["onebody"]["vdw"]
                    }
                curr_equi_par = {}
                for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
                    for i, r in enumerate(self.ric_type[ic]):
                        if self.parind[ic][i] == None:
                            if ((self.atoms_in_subsys(r, curr_fraglist)) and
                                    (self.atoms_in_active(r, curr_atomlist))):
                                # no params yet and in current refsystem => check for params
                                params = []
                                ptypes = []
                                full_parname_list = []
                                # check unsorted list first
                                parname = self.get_parname(r)
                                if parname in curr_par[ic]:
                                    for par in curr_par[ic][parname]:
                                        ptypes.append(par[0])
                                        ### check for equivalences
                                        if par[0] == "equiv":
                                            get_equivalences()
                                            continue
#                                            ### loop over aftypes of ic
#                                            for aidx in r:
#                                                aft = self.aftypes[aidx]
#                                                if aft == par[1][0]:
#                                                    equivs[aidx] = par[1][1]
                                        ### proceed with standard assignment
                                        full_parname = par[0] + "->" + str(parname) + "|" + ref
                                        full_parname_list.append(full_parname)
                                        if not full_parname in self.par[ic]:
                                            self.par[ic][full_parname] = par
                                # now check sorted list (if already in ptype skip), only for manybody ics
                                if ic not in ["cha", "vdw"]:
                                    parname = self.get_parname_sort(r, ic)
                                    if parname in curr_par[ic]:
                                        for par in curr_par[ic][parname]:
                                            if not par[0] in ptypes:
                                                ptypes.append(par[0])
                                                if par[0] == "equiv":
                                                    get_equivalences()
                                                    continue
                                                full_parname = par[0] + "->" + str(parname) + "|" + ref
                                                full_parname_list.append(full_parname)
                                                if not full_parname in self.par[ic]:
                                                    self.par[ic][full_parname] = par
                                if full_parname_list != []:
                                    counter += 1
                                    self.parind[ic][i] = full_parname_list
                                #else:
                                #    print "DEBUG DEBUG DEBUG %s" % ic
                                #    print self.get_parname(r)
                                #    print self.get_parname_sort(r, ic)
                logger.info("%i parameters assigned for ref system %s" % (counter,ref))
                #EQUIVALENCE
                # now all params for this ref have been assigned ... any equivalnce will be renamed now in aftypes
                for i, a in enumerate(copy.copy(self.aftypes)):
                    if i in curr_equi_par.keys():
                        at, ft = curr_equi_par[i].split("@")
                        self.aftypes[i] = aftype(at,ft)
        #self.check_consistency()
        self.setup_pair_potentials()
        self.timer.write_logger(logger.info)
        return

    def check_consistency(self):
        complete = True
        active_zone = []
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            unknown_par = []
            for i, p in enumerate(self.ric_type[ic]):
                if self.parind[ic][i] == None:
                    if ic == "cha" and i not in active_zone: active_zone.append(i)
                    parname = self.get_parname(p)
                    if not parname in unknown_par:
                        unknown_par.append(parname)
            if len(unknown_par) > 0:
                complete = False
                for p in unknown_par: logger.error("No params for %3s %s" % (ic, p))
        if complete == False:
            raise IOError("Assignend parameter set incomplete!")
        else:
            logger.info("Parameter assignment successfull")
        self.setup_pair_potentials()
        self.timer.write_logger(logger.info)
        return

    def write_params_to_key(self):
        self.ric.compute_rics()
        from ff_gen import tools
        kc = tools.keycreator()
        buffer = ""
        ### gather types to make key file more readable
        buffer += "### type lookup table ###\n"
        typemapper =  {}
        for i, aft in enumerate(self.aftypes):
            if str(aft) not in typemapper.keys():
                typemapper[str(aft)] = str(len(typemapper.keys()))
                buffer += "# %-20s %-5s\n" % (str(aft), typemapper[str(aft)])
        buffer += "\n"
        ### gather nonbonded params
        saftypes = map(str, self.aftypes)
        for aft in typemapper.keys():
            idx = saftypes.index(aft)
            buffer += kc.formatter("atom", [typemapper[aft]], params = [elems.mass[self._mol.elems[idx]]])
            #if self.parind["cha"][idx] != None:
            #    buffer += kc.formatter("atom", aft, [elems.mass[self._mol.elems[idx]]])
            #else:
            #    pass

        ### gather bonded params
        icmapper = {"bnd": "bond", "ang": "angle", "dih": "torsion", 
                "oop": "opbend", "cha":"charge", "vdw":"vdw", "strbnd":"strbnd"}
        ics = [ "bnd", "ang", "dih", "oop"]#, "cha", "vdw"]
        data = {"bnd": {}, "ang": {}, "dih": {}, "strbnd": {}, "oop":{}} #, "oop": {}}#, "cha": {}, "vdw": {}}
        for ic in ics:
            for i,p in enumerate(self.ric_type[ic]):
                parname = string.join(map(lambda a: typemapper[a],map(str,self.get_parname_sort(p,ic))), ":")
                if not parname in data[ic]:
                    if self.parind[ic][i] != None:
                        for j in self.parind[ic][i]:
                            if self.par[ic][j][0] == "strbnd":
                                data["strbnd"][parname] = [self.par[ic][j],[[p.value],[],[]]]
                            else:
                                data[ic][parname] = [self.par[ic][j],[[p.value],[],[]]]
                    else:
                        data[ic][parname] = [[],[[p.value],[],[]]]
                else:
                    data[ic][parname][1][0].append(p.value)
        ### write to key
        defaults = {"bnd": 2.0, "ang": 1.0, "dih": 1.0, "oop": 0.1}
        for i in ["bnd", "ang", "strbnd", "dih", "oop",]:
            for ic, p in data[i].items():
                lic = string.split(ic,":")
                if p[0] != []:
                    buffer+= kc.formatter(icmapper[i], lic, params = p[0][1])
                else:
                    if i == "dih":
                        ### tbi
                        pass
                    else:
                        gp = [defaults[i], np.mean(p[1][0])]
                        buffer+= kc.formatter(icmapper[i], lic, params = gp, var = True)
            buffer += "\n"
        print (buffer)
        return data


    
    def setup_pair_potentials(self, radfact = 1.0, radrule="arithmetic", epsrule="geometric"):
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
        for i, t in enumerate(types):
            if t not in self.types2numbers.keys():
                self.types2numbers[t]=str(i)
        ntypes = len(types)
        for i in xrange(ntypes):
            for j in xrange(i, ntypes):
                pair = types[i]+":"+types[j]
                #TODO check availability of an explicit paramerter
                par_i = self.par["vdw"][types[i]][1]
                par_j = self.par["vdw"][types[j]][1]
                pot_i =  self.par["vdw"][types[i]][0]
                pot_j =  self.par["vdw"][types[j]][0]
                if pot_i == pot_j:
                    pot = pot_i
                else:
                    raise IOError("Can not combine %s and %s" % (pot_i, pot_j))
                if radrule == "arithmetic":
                    rad = radfact*(par_i[0]+par_j[0])
                elif radrule == "geometric":
                    rad = 2.0*radfact*np.sqrt(par_i[0]*par_j[0])
                else:
                    raise IOError("Unknown radius rule %s specified" % radrule)
                if epsrule == "arithmetic":
                    eps = 0.5 * (par_i[1]+par_j[1])
                elif epsrule == "geometric":
                    eps = np.sqrt(par_i[1]*par_j[1])
                else:
                    raise IOError("Unknown radius rule %s specified" % radrule)
                self.vdwdata[pair] = (pot,[rad,eps])
        return


    @timer("find reference systems")
    def find_refsystems(self):
        """
        function to detect the reference systems:
            - self.scan_ref      : list of ref names in the order to be searched
            - self.ref_systems   : dictionary of mol objects
            - self.ref_fraglist  : list of fragment indices belonging to this refsystem
            - self.ref_params    : paramtere dictionaries per refsystem (n-body/type)
        """
        if self.source == "mofp":
            self.timer.start("get reference systems")
            scan_ref  = []
            scan_prio = []
            if self._mol.mpi_rank == 0:
                ref_dic = self.api.list_FFrefs(self.FF)
            else:
                ref_dic = []
            if self._mol.mpi_size > 1:
                ref_dic = self._mol.mpi_comm.bcast(ref_dic, root=0)
            for refname in ref_dic.keys():
                prio, reffrags, active, upgrades = ref_dic[refname]
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
            for ref in self.scan_ref:
                if self._mol.mpi_rank == 0:
                    ref_mol = self.api.get_FFref_graph(ref, mol=True)
                else:
                    ref_mol = None
                if self._mol.mpi_size > 1:
                    ref_mol = self._mol.mpi_comm.bcast(ref_mol, root=0)
                ref_mol.addon("fragments")
                ref_mol.fragments.make_frag_graph()
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
                    if ref_dic[ref][2] != None:
                        idx = self.fragments.frags2atoms(subs_flat)
                        self._mol.graph.filter_graph(idx)
                        asubs = self._mol.graph.find_subgraph(self._mol.graph.molg, self.ref_systems[ref].graph.molg)
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
                    ref_par = self.api.get_params_from_ref(self.FF, ref)
                else:
                    ref_par = None
                if self._mol.mpi_size > 1:
                    ref_par = self._mol.mpi_comm.bcast(ref_par, root=0)                
                self.ref_params[ref] = ref_par
                #print ("DEBUG DEBUG Ref system %s" % ref)
                #print (self.ref_params[ref])
            self.timer.stop()
        elif source == "file":
            raise ValueError, "assigning reference systems from file needs to be implemeted"
        else:
            raise ValueError, "unknown source %s" % source
        return


    def atoms_in_subsys(self, alist, fsubsys):
        """
        this helper function checks if all fragments of atoms (indices) in alist
        appear in the list of fragemnts (indices) in fsubsys
        """
        return all(f in fsubsys for f in map(lambda a: self._mol.fragnumbers[a], alist))

    def atoms_in_active(self, alist, subsys):
        """
        this helper function checks if any  atom (indices) in alist
        appear in the list of active atoms (indices) in fsubsys
        """
        if subsys == None: return True
        return any(a in subsys for a in alist)

    def get_parname(self, alist):
        """
        helper function to produce the name string using the self.aftypes
        """
        l = map(lambda a: self.aftypes[a], alist)
        return tuple(l)

    def get_parname_sort(self, alist, ic):
        """
        helper function to produce the name string using the self.aftypes
        """
        l = map(lambda a: self.aftypes[a], alist)
        return tuple(aftype_sort(l,ic))

    def report_params(self):
        """
        TODO: improve this ..write to file or logger in DEBUG level?
        better sorting and commenting?
        """
        for ic in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
            print ("TYPE: " + ic.upper())
            pstrings = self.par[ic].keys()
            pstrings.sort()
            for s in pstrings:
                print (s)
            print ("\n")
        return
