# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 14:19:27 2018

@author: roberto


        addon module acab to implement colors for molsys

        contains class acab
"""
from __future__ import print_function
from __future__ import absolute_import # no RuntimeWarning import error

### LOGGING ###
import logging
logger = logging.getLogger("molsys.acab")
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

### MOLSYS ###
from molsys.addon import base
from molsys.util.color import make_emol, make_vmol, make_mol, ecolor2elem, vcolor2elem
from molsys.util.images import arr2idx
from molsys.util.misc import normalize_ratio, int2base
from molsys.util.sysmisc import isatty, _makedirs, _checkrundir

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

# overload print function in parallel case (python3 compliant [RA])
try:
    import __builtin__
except ImportError:
    import builtins as __builtin__
def print(*args, **kwargs):
    if mpi_rank == 0:
        return __builtin__.print(*args, **kwargs)
    else:
        return

### SCIP ###
try:
    import pyscipopt
except ImportError:
    Model = None
    quicksum = None
else:
    Model = pyscipopt.Model
    quicksum = pyscipopt.quicksum
    quickprod = pyscipopt.quickprod


### UTIL ###
import sys
import os
import numpy as np
from collections import defaultdict, Counter
from itertools import combinations
from numbers import Number
from math import pi, cos
import atexit


if mpi_comm is None:
    logger.error("MPI NOT IMPORTED DUE TO ImportError")
    logger.error(mpi_err)

class acab(base):
    ############################################################################
    ### TODO
    ###
    ### chromo format read/write

    ############################################################################
    ### INITIALIZERS ###

    def __init__(self, mol):
        """
        acab addon object to be attached to the parent mol instance

        :Parameter:
        - mol (obj): mol or mol-derived (bb, topo, etc.) type object
            it follows the addon implementation in molsys
        :Attributes:
         - structures (list of molsys.mol): colored net solutions stored as molecules
         - constrlabels (list of strings): applied constraints stored as strings
         - Model (class): main PySCIPOpt class to setup a model
         - quicksum (function): sum function for PySCIPOpt Expr and ConsExpr
             faster than built-in sum, here as attribute for debugging convenience
         - quickprod (function): product function for PySCIPOpt Expr and ConsExpr
             faster than numpy scalar product, here as attribute for debugging convenience
         - evars (dict): edge variables accessible as:
         (periodic case)
            i (int): source atom index
            j (int): target atom index
            p (int): source-to-target image index (it fulfills periodic edge table convention)
            c (int): color index (see further)
         (non-periodic case)
            i (int): source atom index
            j (int): target atom index
            c (int): color index (see further)
         - vvars (dict): vertex variables accessible as:
            i (int): atom index
            c (int): color index (see further)
        """
        super(acab, self).__init__(mol)
        logger.info("acab addon generated")
        self.structures = [] # list of found structures
        self.constrlabels = [] # list of constraint labels
        # auxiliary: to be available as instance attribute w/o importing (for convenience)
        self.Model = Model # class for constraint integer programming model
        self.quicksum = quicksum # quick sum for variables
        self.quickprod = quickprod # quick product for variable
        self.evars = {}
        self.vvars = {}
        return

    ############################################################################
    ### SETUP (setters w/o getters) ###

    #def setup_color_connectivity(self, otab=None, oconn=None):
    #    """
    #    setup edge colors from color table or color connectivity
    #    
    #    :Parameters:
    #    - otab  (list of ints                 =None): edge color table
    #    - oconn (nested list of lists of ints=None) : edge color connectivity
    #    Provide either otab or oconn. You cannot provide both since one defines
    #    the other. otab and oconn can also be left as None (default) so that
    #    they are set as empty-list attribute (e.g. [] ).
    #    """
    #    pass
    #    #if otab:
    #    #    self.set_otab(otab)
    #    #elif oconn:
    #    #    self.set_oconn(oconn)
    #    #else: # to be set
    #    #    self.otab = []
    #    #    self.oconn = []

    def setup_model(self, verbose=True, debug=False, ctrlc=True, *args, **kwargs):
        """
        Initialize the model and its utility attributes

        :Parameters:
        - verbose (bool=True): output flag for model instance
        - ctrlc (bool=True): experts' flag to en/disable default
            KeyboardInterrupt handling
        :Additional Parameters:
        - args: positional arguments of Model.__init__ (cf. Model)
        - kwargs:  keyword arguments of Model.__init__ (cf. Model)

        :Caveat:
            ctrlc == False:
                - implement a KeyboardInterrupt handler of your own!
                in case of emergency:
                - open a shell
                - ps aux | grep python
                - find the PID of your script
                - kill <yourscript-PID>

        :Toughts:
        - separate edge model and vertex model! (emodel vs vmodel)
        """
        # logo and farewells
        if verbose:
            print_header()
            atexit.register(print_footer)
        self.verbose = verbose
        self.model = Model(*args, **kwargs)
        if not debug:
            self.model.hideOutput()
        self.debug = debug
        if not ctrlc:
            logger.warning("Key Interrupt DISABLED: " \
                "hope you have your own good reasons")
            self.model.setBoolParam('misc/catchctrlc', ctrlc)
        self.ctrlc = ctrlc
        # symmetry enabled check
        try:
            import spglib # no practical use out of this check
            self.sym_enabled = True
        except ImportError:
            logger.error("spglib not imported: symmetry is DISABLED")
            self.sym_enabled = False
        return

    def setup_colors(self, necolors=0, nvcolors=0, *args, **kwargs):
        """
        setup number of edge colors and/or vertex colors.
        if the number of colors is more than zero, the appropriate edge/vertex
            variable setter is called

        :Parameters:
        - necolors (int): number of max edge colors
        - nvcolors (int): number of max vertex colors

        N.B. necolors and nvcolors cannot be:
        - negative (no sense)
        - BOTH zero (why would you ever need a model at all? just set them with
            setters or (e.g.) with <yourinstance>.necolors = colors
        """
        self.assert_nevcolors_number(necolors, nvcolors)
        if necolors > 0:
            self.setup_ecolors(necolors, *args, **kwargs)
        if nvcolors > 0:
            self.setup_vcolors(necolors, *args, **kwargs)
        return

    def setup_ecolors(self, necolors, *args, **kwargs):
        """
        setup number of edge colors
        if the number of edge colors is more than zero, the appropriate edge
            variable setter is called

        :Parameters:
        - necolors (int): number of max edge colors
        N.B. it must be positive

        """
        self.assert_ncolors_number(necolors, strict=True)
        self.set_edge_vars(necolors, set_necolors=True)
        self.setup_vertex2evars()
        self.setup_vertex2edges()
        if necolors > 1:
            label = "ec%d" % necolors
            self.constrlabels.append(label)
        return

    def setup_vcolors(self, nvcolors, *args, **kwargs):
        """
        setup number of vertex colors
        if the number of vertex colors is more than zero, the appropriate vertex
            variable setter is called

        :Parameters:
        - nvcolors (int): number of max vertex colors
        N.B. it must be positive

        """
        self.assert_ncolors_number(nvcolors, strict=True)
        self.set_vertex_vars(nvcolors, set_nvcolors=True)
        self.setup_edge2vertices()
        self.setup_edge2vvars()
        if nvcolors > 1:
            label = "vc%d" % nvcolors
            self.constrlabels.append(label)
        return

    def setup_constraints(self, ecratio=None, vcratio=None):
        """
        DEPRECATED / TBI

        general interface to setup constraints of the model.
        TBI: additional constraints can be set, see tutorial [TBI, easy]

        :Parameters:
        - ecratio (None or list of ints): overall edge   color ratio
        - vcratio (None or list of ints): overall vertex color ratio
        """
        if not None: self.assert_ecratio(ecratio)
        if not None: self.assert_vcratio(vcratio)
        if ecratio:
            self.setup_ecratio(ecratio, set_ecratio=False)
        if vcratio:
            self.setup_vcratio(vcratio, set_vcratio=False)
        self.ecratio = ecratio
        self.vcratio = vcratio

    def setup_ecratio(self, ecratio, vsele=None, set_ecratio=True):
        """
        :Parameters:
        - ecratio (None or list of ints): overall edge color ratio

        TBI: unset ratio with negative integers (convention: -1) and reserve
        the complement of the set elements
        TBI: select atoms which the setup is applied! [JK feature request]
        """
        if not hasattr(self,"necolors"):
            self.setup_ecolors(len(ecratio))
        if vsele is None:
            vertices = range(self._mol.natoms)
        else:
            vertices = vsele
            try:
                vertices[0]
            except TypeError:
                vertices = [vertices]
        self.assert_ecratio(ecratio)
        evars = self.evars
        necolors = self.necolors
        ### ratio ###
        nevars = len(evars) / necolors
        crange = range(necolors)
        ecratio = normalize_ratio(ecratio, nevars)
        etab = self._mol.etab
        ### loop ###
        if self._mol.use_pconn:
            for c in crange:
                self.model.addCons(
                    quicksum([evars[ei,ej,p,c] for ei,ej,p in etab]) == ecratio[c],
                    name = "NEdgeColors(%d)=%d" % (c,ecratio[c])
                )
        else:
            for c in crange:
                self.model.addCons(
                    quicksum([evars[ei,ej,c] for ei,ej in etab]) == ecratio[c],
                    name = "NEdgeColors(%d)=%d" % (c,ecratio[c])
                )
        if set_ecratio: self.ecratio = ecratio
        if len(ecratio) > 0:
            label = "oec" + len(ecratio)*"%d" % tuple(ecratio)
            if vsele is not None:
                label += "s"
            self.constrlabels.append(label)
        return

    def setup_vcratio(self, vcratio, esele=None, set_vcratio=True):
        """
        :Parameters:
        - vcratio (None or list of ints): overall vertex color ratio

        TBI: unset ratio with negative integers (convention: -1) and reserve
        the complement of the set elements without constraint but their total
        sum.
        """
        if not hasattr(self,"nvcolors"):
            self.setup_vcolors(len(vcratio))
        if esele is not None:
            raise NotImplementedError
        #if sele is None:
        #    sele = range(self._mol.natoms)
        #else:
        #    try:
        #        sele[0]
        #    except TypeError:
        #        sele = [sele]
        self.assert_vcratio(vcratio)
        vvars = self.vvars
        nvcolors = self.nvcolors
        ### ratio ###
        nvvars = len(vvars) / nvcolors
        crange = range(nvcolors)
        vcratio = normalize_ratio(vcratio, nvvars)
        ### loop ###
        rvvars = range(nvvars) ### == range(self._mol.natoms)
        for c in crange:
            self.model.addCons(
                quicksum(vvars[v,c] for v in rvvars) == vcratio[c],
                name = "NVertexColors(%d)=%d" % (c,vcratio[c])
            )
        if set_vcratio: self.vcratio = vcratio
        if len(vcratio) > 0:
            label = "ovc" + len(vcratio)*"%d" % tuple(vcratio)
            if esele is not None:
                label += "s"
            self.constrlabels.append(label)
        return

    def setup_ecratio_per_vertex(self, ecratio, vsele=None, set_ecratios=True):
        ### TBI: assert no conflict btw. global and local cratio
        if not hasattr(self,"necolors"):
            self.setup_ecolors(len(ecratio))
        if vsele is None:
            vertices = range(self._mol.natoms)
        else:
            vertices = vsele
            try:
                vertices[0]
            except TypeError:
                vertices = [vertices]
        self.assert_ecratio(ecratio)
        evars = self.evars
        necolors = self.necolors
        crange = range(necolors)
        vertex2edges = self.vertex2edges
        # loop #
        ecratios = []
        for c in crange:
            for v in vertices:
                v2e = vertex2edges[v]
                necratio = normalize_ratio(ecratio,len(v2e))
                ecratios.append(necratio)
                consname = "NEdgeColorsPerVertex(%d,%d)=%d" % (v,c,necratio[c])
                self.model.addCons(
                    quicksum([evars[ei,ej,p,c] for ei,ej,p in v2e]) == necratio[c],
                    name = consname
                )
        if set_ecratios: self.ecratios = ecratios
        if len(ecratio) > 0:
            label = "dec" + len(ecratio)*"%d" % tuple(ecratio)
            if vsele is not None:
                label += "s"
            self.constrlabels.append(label)
        return

    def setup_vcratio_per_edge(self, vcratio, esele=None, set_vcratios=True):
        """
        N.B.: there are always two vertices per edge
        """
        if not hasattr(self,"nvcolors"):
            self.setup_vcolors(len(vcratio))
        if esele is not None:
            raise NotImplementedError
        #if sele is None:
        #    sele = range(self._mol.natoms)
        #else:
        #    try:
        #        sele[0]
        #    except TypeError:
        #        sele = [sele]
        self.assert_vcratio(vcratio)
        vvars = self.vvars
        nvcolors = self.nvcolors
        crange = range(nvcolors)
        edge2vertices = self.edge2vertices
        vcratios = []
        for c in crange:
            for e in edge2vertices:
                e2v = edge2vertices[e]
                nvcratio = normalize_ratio(vcratio,len(e2v))
                vcratios.append(nvcratio)
                if self._mol.use_pconn:
                    consname = "NVertexColorsPerEdge(%d-%d,%d)=%d" % (e[0],e[1],c,nvcratio[c])
                else:
                    consname = "VertexColorRatioPerEdge(%d-%d.%d,%d)=%d" % (e[0],e[1],e[2],c,nvcratio[c])
                self.model.addCons(
                    quicksum([vvars[v,c] for v in e2v]) == nvcratio[c],
                    name = consname
                )
        if set_vcratios: self.vcratios = vcratios
        if len(vcratio) > 0:
            if len(vcratio) > 1:
                label = "dvc" + len(vcratio)*"%d" % tuple(vcratio)
                if esele is not None:
                    label += "s"
            elif esele is None:
                label = ""
            else:
                label = "dvc1s"
            if label != "":
                self.constrlabels.append(label)
        return

    def setup_angle_btw_edges(self, color, theta, sense="min", eps=1e-3, vsele=None):
        """
        Constraint angle between edges to be min/max/close to theta.
        TBI: non-periodic (w/o pconn) version

        IMPLEMENTATION DETAIL
        Instead of explicitly allowing the edge pairs with (range of) angle, it
        forbids the opposite, i.e. edge pairs that do not form/match that angle
        This is easier wrt. integer programming: if edge j and edge k form a
        forbidden angle, then the sum of the edge variabes must be 1 or lower,
        i.e. their colors can't be the same.

        color(int): constraint color (order based on ratio assignment)
        theta(float): angle btw. 0 and pi
        sense(str): sense of constraint [TBI: custom tolerances]
            min:   allowed product btw. edges must be >  cos(theta)
            max:   allowed product btw. edges must be <  cos(theta)
            close: allowed product btw. edges must be =~ cos(theta)
            N.B. theta from 0 to PI and cosine(theta) are antitone mapped,
                so min(theta) -> max(cos(theta))
        eps(float): tolerance to sense
        sele (list of ints or None): selected vertices (if None: all)

        """
        ### assert color < necolors
        if vsele is None:
            vsele = range(self._mol.natoms)
        else:
            try:
                vsele[0]
            except TypeError:
                vsele = [vsele]
        conn = self._mol.conn
        pconn = self._mol.pconn
        if theta < 0 or theta > pi:
            raise TypeError("angle must be btw. 0 and pi")
        cost = cos(theta)
        v2e = [[] for i in range(self._mol.natoms)]
        # compute vectors of selected vertices
        for i in vsele:
            ic = conn[i]
            ip = pconn[i]
            for j,jp in zip(ic,ip):
                rk = self._mol.xyz[j] - self._mol.xyz[i] + np.dot(jp,self._mol.cell)
                dk = np.linalg.norm(rk)
                rk /= dk
                v2e[i].append(rk)
        # compute vectors of sele atoms
        prods = []
        for v_ in v2e:
            if v_:
                v_ = np.array(v_)
                prod = np.dot(v_,v_.T)
                prods.append(prod)
            else:
                prods.append([])
        if sense not in ("min","max","close"):
            raise NotImplementedError("sense \"%s\" not implemented" % sense.__repr__())
        for i,prod in enumerate(prods):
            if len(prod):
                if sense == "min":
                    ea, eb = np.where(prod - eps < cost)
                if sense == "max":
                    ea, eb = np.where(prod + eps > cost)
                if sense == "close":
                    ea, eb = np.where(abs(prod-cost) < eps)
                ew = np.where(ea < eb)
                pairs = zip(ea[ew], eb[ew])
                range_ = range(prod.shape[0])
                # wrong which pairs
                wpairs = [
                    j for j in combinations(range_,2)
                    if j not in pairs
                ]
                # N.B. trick: pairs are cumulated in only one list
                # Easier to code and to read. They are split later
                wpairs_ = sum(wpairs,())
                # wrong j conn triplet (i,j,pconn) pairs
                jpairs_ = [(i,conn[i][j],pconn[i][j]) for j in wpairs_]
                # wrong edge pairs
                epairs_ = [
                    (j,k,arr2idx[p],color) if j < k
                    else (k,j,arr2idx[-p],color)
                    for j,k,p in jpairs_
                ]
                epairs = zip(epairs_[0::2], epairs_[1::2]) # split
                # edge variable pairs
                evpairs_ = [self.evars[e] for e in epairs_]
                evpairs = zip(evpairs_[0::2], evpairs_[1::2]) # split
                # constraints (finally!)
                for iev,(evj,evk) in enumerate(evpairs):
                    ej,ek = epairs[iev] # here just for naming
                    self.model.addCons(
                        evj+evk <= 1, # i.e. they can't be both of the same color!
                        name="WrongEdgePairs(%s-%s.%s,%s;%s-%s.%s,%s)" % \
                            (ej[0],ej[1],ej[2],ej[3],
                             ek[0],ek[1],ek[2],ek[3])
                    )
        senselabels = {"min": "u", "max": "d", "close": "c"}
        label = "e%sax%02d%s" % (color, 10*theta, senselabels[sense])
        if vsele is not None:
            label += "s"
        self.constrlabels.append(label)
        return


    ############################################################################
    ### MAIN ###

    def optimize_model(self):
        self.model.optimize()
        if self.is_optimal_model():
            self.set_colors()
            #if self.evars:
            #    self.set_otab(self.ecolors)
        return

    def free_model(self):
        self.model.freeTransform()
        return

    def is_optimal_model(self):
        return self.model.getStatus() == 'optimal'

    ############################################################################
    ### CYCLE ###

    def cycle_init(self):
        """
        initialize symmetry

        Any symmetry solution subspace belongs to the symmetry space and the 
        subspaces does not overlap each other. In other words, symmetry 
        solution subspaces make a partition of the symmetry solution space.

        :TODO:
        - disable symmetry
        """
        self.colorings = []
        self.cycle_initdir()
        self.cycle_initsym()
        return

    def cycle_initsym(self):
        ### "grey" mol setup (=> symmetry space and permutations) ###
        if self.evars and self.use_edge and self.vvars and self.use_vertex:
            m = make_mol(self._mol, alpha=self.alpha,
                use_edge=self.use_edge, use_vertex=self.use_vertex)
        elif self.evars and self.use_edge:
            m = make_emol(self._mol, alpha=self.alpha)
        elif self.vvars and self.use_vertex:
            m = make_vmol(self._mol)
        else:
            logger.error("step solutions are not constraint: infinite loop!")
            raise TypeError("unconstraint step solutions: infinite loop!")
        if self.init_sym:
            m.addon("spg")
            self.permutations = m.spg.generate_symperms()
        else:
            self.permutations = None
        m.atypes = [0 for i in m.atypes]
        m.write("%s%ssym.mfpx" % (self.rundir, os.sep))
        m.write("%s%ssym.txyz" % (self.rundir, os.sep), pbc=False)
        return

    def cycle_initdir(self):
        """
        create directory to store structures at the end of the loop
        \"colors\" and \"pretty\" contains graph-like structures
        \"colors\" contains useful structures to weave frameworks later
        \"pretty\" contains just clearer views of these structures (do not use)
        \"plain\" contains clearer views in plain tinker format (do not use)
        """
        self.coldir = "%s%scolors" % (self.rundir, os.sep)
        self.predir = "%s%spretty" % (self.rundir, os.sep)
        self.pladir = "%s%splain" % (self.rundir, os.sep)
        _makedirs(self.coldir) #mfpx, w/  pbc, useful
        _makedirs(self.predir) #txyz, w/o pbc, clearer
        _makedirs(self.pladir) #txyz, w/o pbc, clearer, plain
        return

    def cycle_permute(self, ecolors=None, vcolors=None):
        """
        symmetrize solution according to symmetry permutations

        symmetrize the solution of the optimal model in its symmetry solution
        subspace.
        """
        if self.span_sym:
            if self.permutations is not None:
                if ecolors is not None and vcolors is not None:
                    if self.alpha == 2:
                        colors = np.array(ecolors*2 + vcolors)[self.permutations]
                    else:
                        colors = np.array(ecolors + vcolors)[self.permutations]
                elif ecolors is None:
                    colors = np.array(vcolors)[self.permutations]
                elif vcolors is None:
                    if self.alpha == 2:
                        colors = np.array(ecolors)[self.permutations]
                    else:
                        colors = np.array(ecolors*2)[self.permutations]
            else:
                m = self.make_structure(ecolors=ecolors, vcolors=vcolors, alpha=self.alpha)
                m.addon("spg")
                permutations = m.spg.generate_symperms()
                if ecolors is not None and vcolors is not None:
                    if self.alpha == 2:
                        colors = np.array(ecolors*2 + vcolors)[permutations]
                    else:
                        colors = np.array(ecolors + vcolors)[permutations]
                elif ecolors is None:
                    colors = np.array(vcolors)[permutations]
                elif vcolors is None:
                    if self.alpha == 2:
                        colors = np.array(ecolors)[permutations]
                    else:
                        colors = np.array(ecolors*2)[permutations]
        else:
            if ecolors is not None and vcolors is not None:
                if self.alpha == 2:
                    colors = np.array(ecolors*2 + vcolors)[np.newaxis,:]
                else:
                    colors = np.array(ecolors + vcolors)[np.newaxis,:]
            elif ecolors is None:
                colors = np.array(vcolors)[np.newaxis,:]
            elif vcolors is None:
                if self.alpha == 2:
                    colors = np.array(ecolors)[np.newaxis,:]
                else:
                    colors = np.array(ecolors*2)[np.newaxis,:]
        colors = np.vstack({tuple(row) for row in colors})
        return colors

    def cycle_constraint(self):
        ### set solution as negated constraint for the following step
        self.free_model() ### N.B.: to keep original constraints
        if self.use_edge:
            if self.alpha == 2:
                ecolors = self.ecolors[:]
            else:
                ecolors = self.ecolors[:]*2
        if self.use_vertex:
            vcolors = self.vcolors[:]
        ### get permutations of colors ###
        if self.use_edge and self.use_vertex:
            colors = self.cycle_permute(ecolors, vcolors)
        elif self.use_edge:
            colors = self.cycle_permute(ecolors,    None)
        elif self.use_vertex:
            colors = self.cycle_permute(   None, vcolors)
        ### exclude permutations of colors ###
        if self.evars and self.constr_edge:
            if self.alpha == 2:
                ecolorings = colors[:,:len(ecolors)]
            else:
                ecolorings = colors[:,:len(ecolors)/2]
            for ecolors in ecolorings:
                self.exclude_edge_solution(ecolors)
        if self.vvars and self.constr_vertex:
            vcolorings = colors[:,-len(vcolors):]
            for vcolors in vcolorings:
                self.exclude_vertex_solution(vcolors)

    def cycle_step(self):
        """
        step to cycle the model until infeasibility
        the model is optimized, then the new solution/its subspace is set
        as additional negated constraint to the next cycle
        constraints of the previous set (the "original" constraints )are kept
        "freeing the transformation" of the model (cf. SCIP)
        """
        ### optimize model ###
        self.optimize_model()
        if self.is_optimal_model():
            self.cycle_constraint()
            self.colorings.append([self.ecolors[:], self.vcolors[:]])
            return False
        else:
            return True

    def cycle_loop(self, Nmax=1e4, alpha=3, init_sym=True, span_sym=True,
        use_edge=True, use_vertex=True, constr_edge=True, constr_vertex=True,
        color_equivalence=None, rundir="run", newrundir=True):
        """
        cycle model

        :Parameters:
        - Nmax=10.000 (int): maximum number of iterations
        - alpha=3 (int): edge ``quantile'' order (see Notes)
        - init_sym (bool):
        - span_sym (bool):
        - use_edge (bool):
        - use_vertex (bool):
        - constr_edge (bool):
        - constr_vertex (bool):
        - color_equivalence (None or dict):
        - rundir="rundir" (str):
        - newrundir (bool):
        :Attributes:
        - N (int): number of iterations
        :Returns:
        - N (int): abs(N)=number of solutions
        if N > 0: EXACT number of solutions (exhaustiveness)
        if N < 0: MINIMUM number of solutions (there may be others)

        :Notes:
        alpha: for space group symmetry detection, each edge is represented as
            additional atoms in the between of the source vertex and target
            vertex.
            If alpha == 2:
                ONE edge atom is put in the centroid of the two vertices
            If alpha >= 3:
                TWO edge atoms are put in the between of the two vertices as follows:
                1- Bond segment is divided in ``alpha'' equal parts by alpha-1 points
                    between the vertices
                2- The edge atoms lie on the points closest to the vertices, i.e. the
                    first and the (alpha-1)-th
        """
        logger.info("Run exhaustive search of colorings")
        # assert variables
        self.assert_loop_alpha(alpha)
        self.assert_flag(init_sym)
        self.assert_flag(span_sym)
        self.assert_loop_edge(use_edge, constr_edge)
        self.assert_loop_vertex(use_vertex, constr_vertex)
        # set attributes
        self.Nmax = int(Nmax)
        self.alpha = alpha
        self.init_sym = init_sym and self.sym_enabled
        self.span_sym = span_sym and self.sym_enabled
        self.use_edge = use_edge
        self.use_vertex = use_vertex
        self.constr_edge = constr_edge
        self.constr_vertex = constr_vertex
        self.color_equivalence = color_equivalence
        # set run directory
        if newrundir:
            self.rundir = _checkrundir(rundir)
        else:
            self.rundir = rundir
        # initialize loop
        self.cycle_init()
        # run loop
        try:
            logger.info("Key Interrupt DISABLED: " \
                "loop handles CTRL+C")
            self.model.setBoolParam('misc/catchctrlc', True)
            logger.info("Max cycle iteration n.: %s" % self.Nmax)
            N = 0
            while N < self.Nmax:
                self.report_step(N)
                if self.cycle_step():
                    break
                N += 1
        except KeyboardInterrupt:
            N *= -1 ### convention
        if N == self.Nmax:
            N *= -1
        self.report_cycle(N)
        if self.span_sym:
            logger.warning("Symmetry detection is not implemented here")
        #N = self.filter_cycle(N)
        self.set_structures_from_colorings()
        self.write_cycle(abs(N))
        return N

    def report_step(self, i):
        if isatty():
            sys.stdout.write("\r")
        sys.stdout.write("Run cycle iteration n.: %d" % i)
        if not isatty():
            sys.stdout.write("\n")
        sys.stdout.flush()

    def report_last(self, i=None):
        """
        report last cycle status

        :Parameters:
        - i (int): loop index (if None prints None)
        """
        sys.stdout.write("\n")
        status = self.model.getStatus()
        logger.info("Last status: %s; Number of cycles: %s" % (status,i))
        return

    def report_cycle(self, N):
        if N < 0: ### by convention
            N *= -1
            self.report_last(N)
            if N == self.Nmax:
                logger.warning("FAILURE: Iteration limit, convergence NOT reached")
            else:
                logger.warning("FAILURE: Interrupted process, convergence NOT reached")
        else:
            self.report_last(N)
            status = self.model.getStatus()
            if status == "unknown":
                logger.warning("FAILURE: Convergence NOT reached")
            elif status == "timelimit":
                logger.warning("FAILURE: Time limit has been reached")
            elif status == "infeasible":
                if N == 0:
                    logger.info("SUCCESS: No feasible solution found!")
                else:
                    logger.info("SUCCESS: %s feasible solutions found!" % N)
            else:
                logger.error("DISASTER: the unexpected happened!")
        ### timer?
        #print("Total elapsed time: %.3f s" % (self.model.getTotalTime(),))
        return

    def filter_cycle(self, N):
        if self.color_equivalence is None:
            return N
        colorings = self.colorings
        newcolorings = [colorings[0]]
        for j in range(N):
            for k in range(N):
                if j < k:
                    jcolors = colorings[j]
                    kcolors = colorings[k]
                    if not self.color_equivalence(self, jcolors, kcolors):
                        newcolorings.append(kcolors)
        M = len(newcolorings)
        if M < N:
            N = M
            logger.info("%s unequivalent solutions after filtering" % N)
            self.colorings = newcolorings
        return N

    def set_structures_from_colorings(self):
        """
        Set colored net structures from (edge and vertex) colorings
        Structures are stored in `structures` list attribute as mol instances

        """
        for j,coloring in enumerate(self.colorings):
            ecolors, vcolors = coloring
            # N.B.: alpha=2 irrespective to self.alpha by design
            m = self.make_structure(ecolors=ecolors, vcolors=vcolors, alpha=2)
            # helping attributes
            m.coloring = coloring
            m.ecolors = ecolors
            m.vcolors = vcolors
            # store it
            self.structures.append(m)

    def write_cycle(self, N, naming="dummy", labels=None):
        if labels is None and len(self.constrlabels) > 0:
            sorted_labels = sorted(self.constrlabels)
            labels = "-".join(sorted_labels)
        else:
            labels = ""
        naming = "mofplus" ### debug
        fmtsx_naming, fmtdx_naming = self.get_constant_naming(self._mol, N, naming=naming)
        cnt_naming = fmtsx_naming, labels, fmtdx_naming
        cnt_naming = [c for c in cnt_naming if c]
        cst_naming = '_'.join(cnt_naming)
        # auxiliary indexing for alternative naming schemes ()
        if naming in ["acab","mofplus"]:
            spacegroup_numbers = [] # for mofplus
            for jm in self.structures:
                jm.addon("spg")
                jm.spacegroup_number = jm.spg.get_spacegroup()[1]
                spacegroup_numbers.append(jm.spacegroup_number) # for mofplus
            if naming == "mofplus":
                spacegroup_uniquemax = Counter(spacegroup_numbers)
                spacegroup_unique = {k:0 for k in spacegroup_uniquemax.keys()}
                for j,jm in enumerate(self.structures):
                    jspgn = spacegroup_numbers[j]
                    jspgn_max = spacegroup_uniquemax[jspgn]
                    if jspgn_max == 1:
                        spacegroup_label = "%03d" % jspgn
                    else:
                        letterbase = int2base(spacegroup_unique[jspgn], maximum=jspgn_max)
                        spacegroup_label = "%03d%s" % (jspgn, letterbase)
                        spacegroup_unique[jspgn] += 1
                    jm.spacegroup_label = spacegroup_label
        for j,jm in enumerate(self.structures):
            var_naming = cst_naming % self.get_variable_naming(jm, j, N, naming=naming)
            jm.write("%s%s%s.mfpx" % (self.coldir, os.sep, var_naming))
            jm.write("%s%s%s.txyz" % (self.predir, os.sep, var_naming), pbc=False)
            jm.write("%s%s%s.txyz" % (self.pladir, os.sep, var_naming), pbc=False, plain=True)

    @staticmethod
    def get_constant_naming(m, N=1, naming="dummy"):
        """
        Static method to get constant naming for colored net structures
        This part of naming is the same for any of the found solutions

        :Arguments:
        - m (mol): mol object that stores the coloring
        - N=1 (int): total number of found solutions
        - naming="dummy" (str): naming scheme to be used (see comments below)

        :Returns:
        - fmtsx (str): format string at the left side of the constraint labels
        - fmtdx (str): format string at the right side of the constraint labels
        """
        # dummy naming (default): each solution is named after its order of finding
        if naming == "dummy":
            fmtsx = "%%0%dd" % len("%d" % N)
            fmtdx = ""
        # acab naming (standard): each solution is named after its order of finding,
        #   the net name, the consecutive supercell indices, the symmetry spacegroup number
        if naming == "acab":
            if m.supercell:
                scell = ''.join([str(s) for s in m.supercell])
            else:
                scell = '111' # explicit
            m.addon("spg")
            fmtsx = "%%0%dd_%s_%s" % (len("%d" % N), m.name, scell)
            fmtdx = "%03d"
        # mofplus naming (suggested): each solution is named after the net name,
        #   the supercell indices separated by "x", the symmetry spacegroup number and
        #   disambiguation letter/s for the same symmetry spacegroup number
        if naming == "mofplus":
            if m.supercell:
                scell = '_'+'x'.join([str(s) for s in m.supercell])
            else:
                scell = '' # implicit
            m.addon("spg")
            fmtsx = m.name + scell
            fmtdx = "%s"
        return fmtsx, fmtdx

    @staticmethod
    def get_variable_naming(m, j=0, N=1, naming="dummy"):
        """
        Static method to get variable naming for colored net structures
        This part of naming is ensured to be different for any of the found solutions

        :Arguments:
        - m (mol): mol object that stores the coloring
        - j=0 (int): running solution index tracking the order of finding
        - N=1 (int): total number of found solutions
        - naming="dummy" (str): naming scheme to be used (see comments below)

        :Returns:
        - different returns (see comments below)
        """
        # dummy naming (default): the order of finding
        if naming == "dummy":
            return j
        # acab naming (standard): the order of finding and the spacegroup number
        if naming == "acab":
            return j, m.spagroup_number
        # mofplus naming (suggested): the spacegroup number and
        #   disambiguation letter/s for the same symmetry spacegroup number
        if naming == "mofplus":
            return m.spacegroup_label

    def make_structure(self, ecolors=None, vcolors=None, alpha=2):
        """
        Make structures out of edge colors and vertex colors

        :Arguments:
        - ecolors:
        - vcolors:
        - alpha=2 (int):
        """
        if ecolors and vcolors:
            if self.constr_edge:
                ec2e = None
            else:
                ec2e = ecolor2elem[len(ecolor2elem)/2:]+ecolor2elem[:len(ecolor2elem)/2]
            if self.constr_vertex:
                vc2e = None
            else:
                vc2e = vcolor2elem[len(vcolor2elem)/2:]+vcolor2elem[:len(vcolor2elem)/2]
            m = make_mol(self._mol, alpha, ecolors=ecolors, vcolors=vcolors,
                use_vertex=self.use_vertex, use_edge=self.use_edge,
                vc2e=vc2e, ec2e=ec2e)
        elif ecolors:
            m = make_emol(self._mol, alpha, ecolors=ecolors)
        elif vcolors:
            m = make_vmol(self._mol, vcolors=vcolors)
        else:
            logger.error("step solutions are not constraint: infinite loop!")
            raise TypeError("unconstraint step solutions: infinite loop!")
        return m

    ############################################################################
    ### SOLUTION CONSTRAINTS ###

    def validate_edge_constraints(self, ecolors=None):
        if ecolors is None:
            ecolors = self.ecolors
        self.assert_ecolors(ecolors)
        evars = self.evars
        etab = self._mol.etab
        if self._mol.use_pconn:
            for k,c in enumerate(ecolors):
                ei,ej, p = etab[k]
                self.model.addCons(
                    evars[ei,ej,p,c] == 1,
                    name = "ValidateEdgeColor(%d-%d.%d)=%d" % (ei,ej,p,c)
                )
        else:
            for k,c in enumerate(ecolors):
                ei,ej = etab[k]
                self.model.addCons(
                    evars[ei,ej,c] == 1,
                    name = "ValidateEdgeColor(%d-%d)=%d" % (ei,ej,c)
                )
        return

    def validate_vertex_constraints(self, vcolors=None):
        if vcolors is None:
            vcolors = self.vcolors
        self.assert_vcolors(vcolors)
        vvars = self.vvars
        for k,c in enumerate(vcolors):
            self.model.addCons(
                vvars[k,c] == 1,
                name = "ValidateVertexColor(%d)=%d" % (k,c)
            )
        return

    def negate_edge_constraints(self, ecolors=None):
        if ecolors is None:
            ecolors = self.ecolors
        self.assert_ecolors(ecolors)
        evars = self.evars
        etab = self._mol.etab
        if self._mol.use_pconn:
            for k,c in enumerate(ecolors):
                ei,ej, p = etab[k]
                self.model.addCons(
                    evars[ei,ej,p,c] == 0,
                    name = "NegateEdgeColor(%d-%d.%d)=%d" % (ei,ej,p,c)
                )
        else:
            for k,c in enumerate(ecolors):
                ei,ej = etab[k]
                self.model.addCons(
                    evars[ei,ej,c] == 0,
                    name = "NegateEdgeColor(%d-%d)=%d" % (ei,ej,c)
                )
        return

    def negate_vertex_constraints(self, vcolors=None):
        if vcolors is None:
            vcolors = self.vcolors
        self.assert_vcolors(vcolors)
        vvars = self.vvars
        for k,c in enumerate(vcolors):
            self.model.addCons(
                vvars[k,c] == 0,
                name = "NegateVertexColor(%d)=%d" % (k,c)
            )
        return

    def include_edge_solution(self, ecolors=None):
        if ecolors is None:
            ecolors = self.ecolors
        self.assert_ecolors(ecolors)
        esolution = ''.join([str(c) for c in ecolors])
        evars = self.evars
        etab = self._mol.etab
        if self._mol.use_pconn:
            self.model.addCons(
                quicksum([evars[e[0],e[1],e[2],ecolors[k]] for k,e in enumerate(etab)]) >= 1,
                name = "IncludeEdgeSolution(%s)" % esolution
            )
        else:
            self.model.addCons(
                quicksum([evars[e[0],e[1],ecolors[k]] for k,e in enumerate(etab)]) >= 1,
                name = "IncludeEdgeSolution(%s)" % esolution
            )
        return

    def include_vertex_solution(self, vcolors=None):
        if vcolors is None:
            vcolors = self.vcolors
        self.assert_vcolors(vcolors)
        vsolution = ''.join([str(c) for c in vcolors])
        vvars = self.vvars
        self.model.addCons(
            quicksum([vvars[k,c] for k,c in enumerate(vcolors)]) >= 1,
            name = "IncludeVertexSolution(%s)" % vsolution
        )
        return

    def exclude_edge_solution(self, ecolors=None):
        if ecolors is None:
            ecolors = self.ecolors
        self.assert_ecolors(ecolors)
        esolution = ''.join([str(c) for c in ecolors])
        evars = self.evars
        etab = self._mol.etab
        nevars = len(etab)
        if self._mol.use_pconn:
            self.model.addCons(
                quicksum([evars[e[0],e[1],e[2],ecolors[k]] for k,e in enumerate(etab)]) <= nevars - 1,
                name = "ExcludeEdgeSolution(%s)" % esolution
            )
        else:
            self.model.addCons(
                quicksum([evars[e[0],e[1],ecolors[k]] for k,e in enumerate(etab)]) <= nevars - 1,
                name = "ExcludeEdgeSolution(%s)" % esolution
            )
        return

    def exclude_vertex_solution(self, vcolors=None):
        if vcolors is None:
            vcolors = self.vcolors
        self.assert_vcolors(vcolors)
        vsolution = ''.join([str(c) for c in vcolors])
        vvars = self.vvars
        nvvars = self._mol.natoms
        self.model.addCons(
            quicksum([vvars[k,c] for k,c in enumerate(vcolors)]) <= nvvars - 1,
            name = "ExcludeVertexSolution(%s)" % vsolution
        )
        return


    ############################################################################
    ### GET / SET ###

    ### N.B. otab and oconn are meant to be interfaces! They do not constraint
    ###     variables of the model
    #def set_otab(self, otab, set_oconn=True):
    #    self.assert_otab(otab)
    #    self._mol.otab = self.otab = otab
    #    if set_oconn:
    #        self.set_oconn_from_otab(use_otab=False)

    #def set_oconn(self, oconn, set_otab=True):
    #    self.assert_oconn(oconn)
    #    #self.oconn = oconn
    #    self._mol.oconn = self.oconn = oconn
    #    if set_otab:
    #        self.set_otab_from_oconn(use_oconn=False)

    #def get_otab(self):
    #    return self.otab

    #def get_oconn(self):
    #    return self.oconn

    #def set_otab_from_oconn(self, oconn=None, use_oconn=False):
    #    #check/set oconn
    #    if oconn is None:
    #        oconn = self.oconn
    #    elif use_oconn:
    #        self.set_oconn(oconn)
    #    pass

    #def set_oconn_from_otab(self, otab=None, use_otab=False):
    #    # WISHLIST: if pconn: FASTER! (thoughts: use array for same-lenght connectivity)
    #    #check/set otab
    #    if otab is None:
    #        otab = self.otab
    #    elif use_otab:
    #        self.set_otab(otab)
    #    # assign locally to be at hands
    #    # N.B.: DO NOT modify it! They are python lists
    #    conn = self._mol.conn
    #    ctab = self._mol.ctab
    #    if self._mol.use_pconn:
    #        ptab = self._mol.ptab
    #        pconn = self._mol.pconn
    #        ### image index connectivity (to avoid np.array comparison, arguably faster)
    #        iconn = [[arr2idx[j] for j in ic] for ic in pconn]         
    #    # init oconn
    #    oconn = [[None for j in ic] for ic in conn] # init as nested lists of Nones
    #    ### TBI: use etab! ###
    #    if self._mol.use_pconn: # periodic connectivity, ambiguity, slower and harder
    #        for k,(i,j) in enumerate(ctab):
    #            ### i -> j
    #            kp = ptab[k] ### image index
    #            ji_ = [j_ for j_, ji in enumerate( conn[i]) if ji == j ]
    #            jk_ = [j_ for j_, jk in enumerate(iconn[i]) if jk == kp]
    #            if ji_ == jk_:
    #                j_ = ji_
    #            else:
    #                j_ = set(ji_) & set(jk_) # intersect
    #            assert len(j_) == 1, "set of lenght %d, expected lenght is 1" % len(j_)
    #            j_ = list(j_)[0]
    #            oconn[i][j_] = otab[k]
    #            ### j -> i
    #            rp = idx2revidx[kp] ### reverse image index (by convention)
    #            ij_ = [i_ for i_, ij in enumerate( conn[j]) if ij == i ]
    #            ik_ = [i_ for i_, ik in enumerate(iconn[j]) if ik == rp]
    #            if ij_ == ik_:
    #                i_ = ij_
    #            else:
    #                i_ = set(ij_) & set(ik_) # intersect
    #            assert len(i_) == 1, "set of lenght %d, expected lenght is 1" % len(i_)
    #            i_ = list(i_)[0]
    #            oconn[j][i_] = otab[k]
    #    else: # no periodic connectivity, no ambiguity, faster and easier
    #        for k,(i,j) in enumerate(ctab):
    #            j_ = conn[i].index(j)
    #            oconn[i][j_] = otab[k]
    #            i_ = conn[j].index(i)
    #            oconn[j][i_] = otab[k]
    #    self.set_oconn(oconn, set_otab=False)
    #    return

    def set_edge_vars(self, necolors, set_necolors=True):
        self.assert_ncolors_number(necolors, strict=True)
        crange = range(necolors)
        etab = self._mol.etab
        evars = {} # dict[(3or4)-tuple] = var... inefficient yet clear
        if self._mol.use_pconn: # periodic connectivity, ambiguity
            for ei,ej,p in etab:
                _evars = [] # temp. list for later constraint (q.v.)
                # add edge variables to model
                for c in crange:
                    evar = self.model.addVar(
                        vtype="B",
                        name = "e(%d-%d.%d,%d)" % (ei,ej,p,c)
                    )
                    _evars.append(evar)
                    evars[ei,ej,p,c] = evar
                # set color uniqueness (one color per edge)
                self.model.addCons(
                    quicksum(_evars) == 1,
                    name = "EdgeColorsUniqueness(%d-%d.%d)" % (ei,ej,p)
                )
        else: # no periodic connectivity, no ambiguity
            for ei,ej in etab:
                _evars = [] # temp. list for later constraint (q.v.)
                # add edge variables to model
                for c in crange:
                    evar = self.model.addVar(
                        vtype="B",
                        name="e(%d-%d,%d)" % (ei,ej,c)
                    )
                    _evars.append(evar)
                    evars[ei,ej,c] = evar
                # set color uniqueness (one color per edge)
                self.model.addCons(
                    quicksum(_evars) == 1,
                    name = "EdgeColorUniqueness(%d-%d)" % (ei,ej)
                )
        self.evars = evars
        if set_necolors: self.necolors = necolors
        return

    def get_edge_vars(self):
        return self.evars

    def set_vertex_vars(self, nvcolors, set_nvcolors=True):
        self.assert_ncolors_number(nvcolors, strict=True)
        ratoms = range(self._mol.natoms)
        crange = range(nvcolors)
        vvars = {} # dict(2-tuple->var)inefficient yet clear
        for v in ratoms:
            _vvars = []
            for c in crange:
                vvar = self.model.addVar(
                    vtype="B",
                    name="v(%d,%d)" % (v,c)
                )
                _vvars.append(vvar)
                vvars[v,c] = vvar
            # set color uniqueness (one color per vertex)
            self.model.addCons(
                quicksum(_vvars) == 1,
                name = "VertexColorUniqueness(%d)" % v
            )
        self.vvars = vvars
        if set_nvcolors: self.nvcolors = nvcolors
        return

    def get_vertex_vars(self):
        return self.vvars

    def set_edge_colors(self):
        ecolors = []
        evars = self.evars
        if evars:
            necolors = self.necolors
            crange = range(necolors)
            etab = self._mol.etab
            if self._mol.use_pconn:
                for ei,ej,p in etab:
                    try:
                        vals = [self.model.getVal(evars[ei,ej,p,c]) for c in crange]
                        val = vals.index(1)
                    except ValueError:
                        vals = [round(v) for v in vals]
                        val = vals.index(1)
                    except Warning:
                        val = "u"
                    ecolors.append(val)
            else:
                for ei,ej in etab:
                    try:
                        vals = [self.model.getVal(evars[ei,ej,c]) for c in crange]
                        val = vals.index(1)
                    except ValueError:
                        vals = [round(v) for v in vals]
                        val = vals.index(1)
                    except Warning:
                        val = "u"
                    ecolors.append(val)
        self.ecolors = ecolors
        return

    def get_edge_colors(self):
        return self.ecolors

    def set_vertex_colors(self):
        vcolors = []
        vvars = self.vvars
        if vvars:
            nvcolors = self.nvcolors
            nvvars = len(vvars) / nvcolors
            crange = range(nvcolors)
            rvvars = range(nvvars) ### == range(self._mol.natoms)
            for v in rvvars:
                try:
                    vals = [self.model.getVal(vvars[v,c]) for c in crange]
                    val = vals.index(1)
                except ValueError:
                    vals = [round(v) for v in vals]
                    val = vals.index(1)
                except Warning:
                    val = "u"
                vcolors.append(val)
        self.vcolors = vcolors
        return

    def get_vertex_colors(self):
        return self.vcolors

    def set_colors(self):
        """to be called after model optimization
        """
        if self.is_optimal_model():
            self.set_edge_colors()
            self.set_vertex_colors()
        else:
            raise Warning("This method cannot be called at this stage")

    def get_colors(self):
        return self.ecolors, self.vcolors

    ############################################################################
    ### DOWNLOAD / UPLOAD ###

    def upload_():
        pass

    def download_():
        pass

    def find_():
        pass

    ############################################################################
    ### READ / WRITE ###

    def read_():
        pass

    def write_():
        pass

    ############################################################################
    ### UTILS / MISC ###

    def setup_vertex2edges(self, sort_flag=True):
        """
        setup dictionary from vertices to list of edges (convenience)
        keys are integers and values are list of (i,j,p) tuples
        where i,j,p are atom i, atom j, and periodic image index respectively

        :Parameters:
        - sort_flag (bool=True): if True: sorts lists of variables according to
            names of variables
        """
        # init #
        vertex2edges = defaultdict(list)
        evars = self.evars
        # loop #
        for ke in evars:
            if ke[-1] == 0: # easier than reconstruct again the bonds w/ pconn
                vertex2edges[ke[0]].append(ke[:-1])
                vertex2edges[ke[1]].append(ke[:-1])
        # sort #
        if sort_flag:
            for kv2e in vertex2edges:
                vertex2edges[kv2e].sort()
        self.vertex2edges = dict(vertex2edges)
        return

    def setup_edge2vertices(self, sort_flag=True):
        """
        setup dictionary from edges to list of vertices (convenience)
        keys are (i,j,p) tuples and values are list of integers
        where i,j,p are atom i, atom j, and periodic image index respectively

        :Parameters:
        - sort_flag (bool=True): if True: sorts lists of variables according to
            names of variables
        """
        # init #
        edge2vertices = defaultdict(list)
        vvars = self.vvars
        vcrange = range(self.nvcolors)
        etab = self._mol.etab
        # loop #
        if self._mol.use_pconn:
            for ei,ej,p in etab:
                edge2vertices[ei,ej,p] = [ei,ej]
        else:
            for ei,ej in etab:
                edge2vertices[ei,ej,p] = [ei,ej]
        # sort #
        if sort_flag: # useful only if use_pconn, otherwise sorted by design
            for ke2v in edge2vertices:
                edge2vertices[ke2v].sort()
        self.edge2vertices = dict(edge2vertices)
        return

    def setup_vertex2evars(self, sort_flag=True):
        """
        setup dictionary from vertices to list of edge variables (convenience)
        keys are integers and values are list of SCIP variables

        :Parameters:
        - sort_flag (bool=True): if True: sorts lists of variables according to
            names of variables
        """
        # init #
        vertex2evars = defaultdict(list)
        evars = self.evars
        # loop #
        for ke in evars:
            vertex2evars[ke[0]].append(evars[ke])
            vertex2evars[ke[1]].append(evars[ke])
        # sort #
        if sort_flag:
            for kv2e in vertex2evars:
                vertex2evars[kv2e].sort(key=lambda x: x.name)
        self.vertex2evars = dict(vertex2evars)
        return

    def setup_edge2vvars(self, sort_flag=True):
        """
        setup dictionary from edges to list of vertex variables (convenience)
        keys are (i,j,p) tuples and values are list of SCIP variables
        where i,j,p are atom i, atom j, and periodic image index respectively

        :Parameters:
        - sort_flag (bool=True): if True: sorts lists of variables according to
            names of variables
        """
        # init #
        edge2vvars = defaultdict(list)
        vvars = self.vvars
        evars = self.evars # easiest to get etab # TBI: IN MOLSYS
        vcrange = range(self.nvcolors)
        etab = self._mol.etab
        # loop #
        if self._mol.use_pconn:
            for ei,ej,p in etab:
                edge2vvars[ei,ej,p] += [vvars[ei,c] for c in vcrange]
                edge2vvars[ei,ej,p] += [vvars[ej,c] for c in vcrange]
        else:
            for ei,ej in etab:
                edge2vvars[ei,ej] += [vvars[ei,c] for c in vcrange]
                edge2vvars[ei,ej] += [vvars[ej,c] for c in vcrange]
        # sort #
        if sort_flag: # useful only if use_pconn, otherwise sorted by design
            for ke2v in edge2vvars:
                edge2vvars[ke2v].sort(key=lambda x: x.name)
        self.edge2vvars = dict(edge2vvars)
        return

    ############################################################################
    ### TODO ###

    def __set_evars_str():
        def __evars_str():
            pass
        evars.__str__ = __evars_str

    def _set_vvars_str(self, vvars):
        ###TBI###
        ###BUG: AttributeError: 'dict' object attribute '__str__' is read-only
        def __vvars_str():#*args, **kwargs):
            keys = sorted(vvars.keys())
            _strl = []
            for k in keys:
                var = vvars[k]
                try:
                    val = self.model.getVal(var)
                except Warning: # check getVal for more info
                    val = "u" # stands for unsoled
                _str = "%s=%s" % (var.__str__(), val.__str__())
                _strl.append(_str)
            return ', '.join(_strl)
        vvars.__str__ = __vvars_str

    ############################################################################
    ### REPORT ###

    def report_edge_values(self):
        if not self.is_optimal_model():
            return
        evars = self.evars
        if evars:
            for var in evars.values():
                try:
                    val = int(self.model.getVal(var))
                except Warning: #model may be unsolved
                    val = "?" #stands for unsolved
                except ValueError:
                    try:
                        vals = [round(v) for v in vals]
                        val = vals.index(1)
                    except ValueError:
                        val = "!"
                print(var,val)
        return

    def report_vertex_values(self):
        if not self.is_optimal_model():
            return
        vvars = self.vvars
        if vvars:
            for var in vvars.values():
                try:
                    val = int(self.model.getVal(var))
                except Warning: #model may be unsolved
                    val = "?"
                except ValueError: #1 is not found!
                    try:
                        vals = [round(v) for v in vals]
                        val = vals.index(1)
                    except ValueError:
                        val = "!"
                print(var,val)
        return

    def report_values(self):
        ### TBI: needs formatting improvement! ###
        self.report_edge_values()
        self.report_vertex_values()
        return

    def report_edge_colors(self):
        if not self.is_optimal_model():
            return
        evars = self.evars
        if evars:
            necolors = self.necolors
            crange = range(necolors)
            etab = self._mol.etab
            if self._mol.use_pconn:
                for ei,ej,p in etab:
                    try:
                        vals = [self.model.getVal(evars[ei,ej,p,c]) for c in crange]
                        val = vals.index(1)
                    except Warning: #model may be unsolved
                        val = "?"
                    except ValueError: #1 is not found!
                        try:
                            vals = [round(v) for v in vals]
                            val = vals.index(1)
                        except ValueError:
                            val = "!"
                    print("%d-%d.%d = %s" % (ei,ej,p,val))
            else:
                for ei,ej in etab:
                    try:
                        vals = [self.model.getVal(evars[ei,ej,c]) for c in crange]
                        val = vals.index(1)
                    except Warning: #model may be unsolved
                        val = "?"
                    except ValueError: #1 is not found!
                        val = "!"
                    print("%d-%d = %s" % (ei,ej,val))
        return

    def report_vertex_colors(self):
        if not self.is_optimal_model():
            return
        vvars = self.vvars
        if vvars:
            nvcolors = self.nvcolors
            nvvars = len(vvars) / nvcolors
            crange = range(nvcolors)
            rvvars = range(nvvars) ### == range(self._mol.natoms)
            for v in rvvars:
                try:
                    vals = [self.model.getVal(vvars[v,c]) for c in crange]
                    val = vals.index(1)
                except Warning:
                    val = "?"
                except ValueError: #1 is not found!
                    try:
                        vals = [round(v) for v in vals]
                        val = vals.index(1)
                    except ValueError:
                        val = "!"
                print("%d = %s" % (v,val))
        return

    def report_colors(self):
        ### TBI: needs formatting improvement! ###
        if not self.is_optimal_model():
            return
        self.report_edge_colors()
        self.report_vertex_colors()
        return

    ############################################################################
    ### ASSERT ###
    ### TBI?: assertive decorators? ###

    def assert_flag(self,flag):
        """
        Check if flag is boolean

        flag (object): instance to be checked
        """
        assert isinstance(flag,bool), "flag must be boolean" \
            "you gave %s and it is %s" % (flag, flag.__class__.__name__)
        return

    #def assert_otab_or_oconn(self,otab,oconn):
    #    """check otab/oconn arguments"""
    #    self.assert_otab(otab)
    #    self.assert_oconn(oconn)
    #    assert otab is None or oconn is None, "give colors as either color "\
    #        "table or color connectivity"

    def assert_otab(self,otab):
        """
        Check length of color table # EDGE COLORS

        otab (list of integers): color table to be compared against edge table

        """
        assert len(otab) == len(self._mol.ctab)
        return

    def assert_oconn(self,oconn):
        """
        Check length of color connectivity # VERTEX COLORS

        """
        assert len(oconn) == self._mol.natoms
        return

    def assert_hasmodel(self):
        """
        Check if model is setup in acab

        TBI: when edge model and vertex model will be separated (here for future)
        """
        assert hasattr(self, "model"), "model must be setup before this method"
        return
    
    def assert_ncolors_number(self, ncolors, strict=False):
        """
        Check non-negative number of colors

        ncolors (int): number of colors
        strict=False (bool): strict mode
            if True, ncolors must be positive (0 is not allowed)
            if False, ncolors must be non-negative (0 is allowed)
        """
        if strict:
            assert ncolors > 0,   "number of colors is %s but must be positive" % ncolors
        else:
            assert ncolors >= 0,   "number of colors is %s but cannot be negative" % ncolors
        return

    def assert_nevcolors_number(self, necolors, nvcolors):
        """
        Check non-negative number of edge colors and vertex colors. At least
        one must be positive

        necolors (int): number of edge colors
        nvcolors (int): number of vertex colors

        Used when numbers of both edge colors and vertex colors are set with
            the same method

        """
        self.assert_ncolors_number(necolors)
        self.assert_ncolors_number(nvcolors)
        assert necolors > 0 or nvcolors > 0, \
            "At least one btw. edge and vertex colors must be positive: " \
            "you gave %s and %s" % (necolors, nvcolors)
        return

    def assert_ecratio(self, ecratio):
        """
        Check edge color ratio
        The ratio must be an iterable of numbers with the same length as
            number of edge colors

        ecratio (iterable): edge color ratio

        """
        assert hasattr(ecratio, "__iter__"), "color ratio must be iterable, " \
            "you gave %s" % ecratio.__class__.__name__
        assert len(ecratio) == self.necolors, "given edge ratio %s is " \
            "inconsistent with maximum number of edge colors %s" % (ecratio, self.necolors)
        for weight in ecratio:
            assert isinstance(weight,Number), "color weight in color ratio " \
                "must be number, you gave %s and it is %s" % (weight, weight.__class__.__name__)
        return

    def assert_vcratio(self, vcratio):
        """
        Check vertex color ratio
        The ratio must be an iterable of numbers with the same length as
            number of vertex colors

        vcratio (iterable): vertex color ratio

        """
        assert hasattr(vcratio, "__iter__"), "color ratio must be iterable, " \
            "you gave %s" % vcratio.__class__.__name__
        assert len(vcratio) == self.nvcolors, "given vertex ratio %s is " \
            "inconsistent with maximum number of vertex colors %s" % (vcratio, self.nvcolors)
        for weight in vcratio:
            assert isinstance(weight, Number), "color weight in color ratio " \
                "must be number, you gave %s and it is %s" % (weight, weight.__class__.__name__)
        return

    def assert_ecolors(self, ecolors):
        """
        Check length of edge colors

        ecolors (list of ints): edge colors

        """
        assert len(ecolors) == len(self._mol.etab), "length of edge colors and " \
        "number of edges mismatch: %d vs %d" % (len(ecolors), len(self._mol.etab))
        return

    def assert_vcolors(self, vcolors):
        """
        Check length of vertex colors

        vcolors (list of ints): vertex colors

        """
        assert len(vcolors) == self._mol.natoms,"length of vertex colors and " \
        "number of vertex mismatch: %d vs %d" % (len(vcolors), self._mol.natoms)
        return

    def assert_loop_alpha(self, alpha):
        """
        Check alpha order of edge quantile

        alpha (float): alpha

        """
        assert alpha >= 2, "alpha must be >= 2"
        return

    def assert_loop_edge(self, use_edge, constr_edge):
        self.assert_flag(use_edge)
        self.assert_flag(constr_edge)
        assert use_edge or not constr_edge, "edges cannot be constraint"\
            " if they are not used"
        return

    def assert_loop_vertex(self, use_vertex, constr_vertex):
        self.assert_flag(use_vertex)
        self.assert_flag(constr_vertex)
        assert use_vertex or not constr_vertex, "vertex cannot be constraint"\
            " if they are not used"
        return

    ############################################################################
    ### TEST ###
    ### notes [RA]:
    ### - if constraints are maintained as intended before/after optimization
    ###     before/after freeTransform -> CHECK SAME CONSTRAINTS
    ### - enable none/one/either/both of necolors and nvcolors -> OK, expected?
    ### - integer/fractional/zero ratios -> OK
    ### - give negative ratios -> LEAVE UNSET, KEEP THE RATIO OUT
    ### - check in general if ratios are kept globally and locally
    ### - n.b. take care of fractional, zero and negative ratios
    ### - when a new setup_(e/v)cratio is given, it may conflict with
    ###     previously set n(e/v)colors. Set edge/vertex new variables may be
    ###     initialized again if n(e/v)colors are changed.
    def test_():
        pass

    ############################################################################
    ### DEBUG ###

    def debug_():
        pass

__version__ = "2.1.1"

stars= """
********************************************************************************"""
logo = """
                                                          .k:
                                                          XMk
                                                         xMMk
                                                        :MMMx
                                                       .WMMMd
                                                       KMMMMo
                                                      xMMMMMl
                                                     :MMMMMM:
                                                    .WMMMMMM;
                                                    KMMMMMMM,
                                                   dMMMMMMMM'
                           kko;.                  ,MMMMMMMMM.         .:oko
                          .OOOOOkl.  ...',,;;;:; .NMMMMMMMMM.      .lkOOOOk
                          .OOOOOOOOo.,kOOOOOOOO, OMMMMMMMMMM .Okdl:'..,cxOk
                        '' kOOOOOOOOk..xOxo:,'. lMMMMMMMMMMM 'OOOOOOOOo;..;
                    'cxOOk ,OOOOOOOOOO, .':ldd 'MMMMMMMMMMMW ,OOOOOOOOOOOx'
                 ,oOOOOOOOd ,OOOOOOOOOOxOOOOO. NMMMMMMMMMMMX .ckOOOOOOOOOOOl
              'oOOOOOOOOOOOx..dOOOOOOOOOOOOO; kMMMMMMMMMMMMK '; ;OOOOOOOOOOOo
            ;kOOOOOOOOOOOOOOOl..dOOOOOOOOOOo cMMMMMMMMMMMMM0 :Ox..OOOOOOOOOOO;
          ;kOOOOOOOOOOOOOOOOOOx ,OOOOOOOOOk..WMMMMMMMMMMMMMk cOOx ,OOOOOOOOOOk
        'xOOOOOOOOOOOOOOOOOOOO, kOOOOOOOOO, XMMMMMMMMMMMMMMx lOOOc dOOOOOOOOOO'
       lOOOOOOOOOOOOOOOOOOOOOk..OOOOOOOOOk .0XNWMMMMMMMMMMMd lOOOk ,OOOOOOOOOO;
     .xOOOOOOOOOOOOOOOOOOOOOd  ,OOOOOOOOOOx;.................kOOOO..OOOOOOOOOO:
    .kOOOOOOOOOOOOOOOOOOOOOc   ;OOOOOOOOOOOOOOOOOOOkxxdoolloOOOOOO.'OOOOOOOOOO;
   .OOOOOOOOOOOOOOOOOOOOOOc    'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO ;OOOOOOOOOO.
  .kOOOOOOOOOOOOOOOOOOOOOo      OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx cOOOOOOOOOk
  xOOOOOOOOOOOOOOOOOOOOOx       lOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOl oOOOOOOOOO:
 ;OOOOOOOOOOOOOOOOOOOOOO.       .OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOc xOOOOOOOOk
 kOOOOOOOOOOOOOOOOOOOOOo         ;OOOOOO:     OOOOOOOOOO     lOOO: xOOOOOOOO,
,OOOOOOOOOOOOOOOOOOOOOO'          lOOOOO,      .kOOOOd.      cOOOl oOOOOOOOl
lOOOOOOOOOOOOOOOOOOOOOO            lOOOOk,    'dOOOOOOo.    ;OOOOk :OOOOOOx
xOOOOOOOOOOOOOOOOOOOOOk             cOOOOOOxxkOOOOOOOOOOkxxOOOOOO; .OOOOOO.
OOOOOOOOOOOOOOOOOOOOOOx              ,kOOOOOOOOOOOOOOOOOOOOOOOOk.   kOOOOc
OOOOOOOOOOOOOOOOOOOOOOO               .dOOOOOOOOOOOOOOOOOOOOOOo     ;OOOO,
kOOOOOOOOOOOOOOOOOOOOOO.                oOOOOOOOOOOOOOOOOOOOOc       :OOO'
oOOOOOOOOOOOOOOOOOOOOOOl                 dOOOOOOOOOOOOOOOOOOl         'xO.
,OOOOOOOOOOOOOOOOOOOOOOO.                 xOOOOOOOOOOOOOOOOx
 kOOOOOOOOOOOOOOOOOOOOOOx                 .OOOklcxOOdcoOOOO.
 :OOOOOOOOOOOOOOOOOOOOOOOx                 oOO;   Ox   cOOc
  kOOOOOOOOOOOOOOOOOOOOOOOx.               'OOklcxOOdclOOO.
  .OOOOOOOOOOOOOOOOOOOOOOOOO:               ;kOOOOOOOOOOx'
   ,OOOOOOOOOOOOOOOOOOOOOOOOOk;               ':odxxdl:.
    ,OOOOOOOOOOOOOOOOOOOOOOOOOOkc.                                   ..
     'kOOOOOOOOOOOOOOOOOOOOOOOOOOOx:.                               ,k;
      .dOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOd:'.                        ;xOx
        ;OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxoc;,...          .'cdOOOo
          cOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOkkkkkOOOOOOOO:
           .cOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOd.
              ,dOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOd.
                .;oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOx:.
                    .:okOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOxl,.
                        ..;coxkOOOOOOOOOOOOOOOOOkdoc;'.
                                 .;coxkOdoc;'.
"""
info = """
        .,;;'.              ..oOOOOo;            .';;,.         ..oOOOOOoo.     
    .oKMMMMMMMNx'       ,lxOOOOOOOOOOk:      :kNMMMMMMM0c      xOOOOOOOOOOOxl'  
  .0MMMMMMMMMMMMMO    :kOOOx.       ,OO,   oWMMMMMMMMMMMMN;    OOc       cOOOOl
 oMMMNdc;'..'oMMMMX  dOOOOo         .OO, .XMMMOl;,...:KMMMM;  ;OO:       xOOOOx
:MMMX         dMMMM,lOOOOk          .Ok  NMMM:        .MMMMO  kOO,     .dOOOx:  
OMMM:         lMMMM.kOOOOo           O' ,MMMX          MMMMd ;OOOdllodkOkd:.    
OMMMd,,,,'''''KMMMk xOOOOk           .  ,MMMX,,,,,'''.lMMMW. kOOc      .kOOkc   
xMMWXNNNNWWWWWMMMN. .OOOOOd.            .MMMNXNNNNWWWWMMMMl 'OOO.       kOOOOx  
lMM:          cMMl   'kOOOOOd:'.   .,.   MM0           NMX  ;OOx       cOOOOOd  
lMo            WW      ;dOOOOOOOOOOo.    MX.           dMc   oOO;.  ..lOOOOo,   
ck             od        .:dkOOxl;       X.            .N.    .lxkOOkkxdl;      


                        ACAB = ALL COLORS ARE BEAUTIFUL
             by Roberto Amabile <roberto d0t amabile at rub d0t de>

  Net coloring + advanced Reverse Topological Approach: R. Amabile, R. Schmid
             (C) 2018- Computational Materials Chemistry (Germany).

********************************************************************************
"""

footer="""
********************************************************************************

                   Your ride with ACAB ends here. Houyhnhnm!

********************************************************************************
"""

def print_header():
    print(stars)
    if isatty():
        print(logo)
    print(info)

def print_footer():
    print(footer)

