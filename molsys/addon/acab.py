# -*- coding: utf-8 -*-
from __future__ import print_function
"""
Created on Mon Jun 11 14:19:27 2018

@author: roberto


        addon module acab to implement colors for molsys
        
        contains class acab
"""
from molsys.addon import base

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

### UTIL ###
from molsys.util.images import arr2idx
from molsys.util.images import idx2revidx
from collections import defaultdict
from numbers import Number

### LOGGING ###
import logging
logger = logging.getLogger("molsys.acab")
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)
### HACK ###******
logger.debug = print
logger.info = print
logger.warning = print
logger.error = print
logger.critical = print
#*****************
import sys

if mpi_comm is None:
    logger.error("MPI NOT IMPORTED DUE TO ImportError")
    logger.error(mpi_err)

#ctab: connectivity table
#ptab: periodic image table
#otab: EDGE color table
class acab(base):
    ############################################################################
    ### TODO
    ###
    ### cycle
    ### read/write
    ### testing cases! (w/ & w/o pconn)

    ############################################################################
    ### INITIALIZERS ###

    def __init__(self, mol, otab=None, oconn=None, *args, **kwargs):
        """
        acab addon object to be attached to the parent mol instance
        
        :Parameter:
            - mol (obj): mol or mol-derived (bb, topo, etc.) type object
            - otab (list of couples): color table
            - oconn (list of lists): color connectivity
        Provide either otab or oconn. You cannot provide both since one defines
            the other.
        otab and oconn can be left None so that they are set as empty list
            attribute (e.g. [] ).
        :Attributes:
            - Model (class): main PySCIPOpt class to setup a model
            - quicksum (function): sum function for PySCIPOpt Expr and ConsExpr
                faster than built-in sum
        """
        super(acab, self).__init__(mol)
        ### default settings
        acab.Model = Model
        acab.quicksum = quicksum
        self.set_etab() ### TBI: MUST GO ELSEWHERE... IN MOL INSTANCE
        ### set colors from color table or color connectivity
        # one at least must be None
        assert otab is None or oconn is None, "give colors as either color table or color connectivity"
        if otab:  #former not None, latter None
            self.set_otab(otab)
        elif oconn: #former None, latter not None
            self.set_oconn(oconn)
        else: #both None
            ### to be initialized
            self.otab = []
            self.oconn = []
        ### debug
        logger.debug("acab addon generated")
        return

    ############################################################################
    ### SETUP (setters w/o getters) ###

    def setup_model(self, verbose=False, ctrlc=True, *args, **kwargs):
        """
        initialize the model and its utility attributes

        :Parameters:
        - verbose (bool=True): output flag for model instance
        - ctrlc (bool=True): experts' flag to en/disable default
            KeyboardInterrupt handling
        :Additional Parameters:
        - args: positional arguments of Model.__init__ (cf. Model)
        - kwargs:  keyword arguments of Model.__init__ (cf. Model)
        
        if ctrlc == False: implement a KeyboardInterrupt handler of your own!
        in case of emergency:
        - open a shell
        - ps aux | grep python
        - find the PID of your script
        - kill <yourscript-PID>

        :Toughts:
        - separate edge model and vertex model! (emodel vs vmodel)
        """
        self.model = Model(*args, **kwargs)
        self.assert_flag(ctrlc)
        if not verbose:
            self.model.hideOutput()
        if not ctrlc:
            logger.warning("KeyboardInterrupt is DISABLED: " \
                "hope you have your own good reasons")
            self.model.setBoolParam('misc/catchctrlc', ctrlc)
        self.evars = {}
        self.vvars = {}
    
    def setup_colors(self, necolors=0, nvcolors=0, *args, **kwargs):
        """
        setup number of edge colors and/or vertex colors.
        if a color is more than zero, the appropriate variable setter is called
        
        :Parameters:
        - necolors (int): number of max edge   colors
        - nvcolors (int): number of max vertex colors

        N.B. necolors and nvcolors cannot be:
        - negative (no sense)
        - BOTH zero (why would you ever need a model at all? just set them with
            setters or (e.g.) with <yourinstance>.necolors = colors
        TBI: THEY CAN BE NONE SO THAT IT IS RETRIEVED FROM CRATIO (?!)
        """
        self.assert_hasmodel()
        self.assert_nevcolors_number(necolors, nvcolors)
        if necolors > 0:
            self.set_edge_vars(necolors, set_necolors=True)
            self.setup_vertex2evars()
            self.setup_vertex2edges()
        if nvcolors > 0:
            self.set_vertex_vars(nvcolors, set_nvcolors=True)
            self.setup_edge2vertices()
            self.setup_edge2vvars()

    def setup_constraints(self, ecratio=None, vcratio=None):
        """
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

    def setup_ecratio(self, ecratio, set_ecratio=True):
        """
        :Parameters:
        - ecratio (None or list of ints): overall edge   color ratio

        TBI: unset ratio with negative integers (convention: -1) and reserve
        the complement of the set elements
        """
        if not hasattr(self,"necolors"):
            self.setup_colors(necolors=len(ecratio))
        self.assert_ecratio(ecratio)
        evars = self.evars
        necolors = self.necolors
        ### ratio ###
        nevars = len(evars) / necolors
        crange = range(necolors)
        ecratio = self.normalize_ratio(ecratio, nevars)
        etab = self.etab
        ### loop ###
        if self._mol.use_pconn:
            for c in crange:
                self.model.addCons(
                    quicksum([evars[ei,ej,p,c] for (ei,ej),p in etab]) == ecratio[c],
                    name = "NEdgeColors(%d)=%d" % (c,ecratio[c])
                )
        else:
            for c in crange:
                self.model.addCons(
                    quicksum([evars[ei,ej,c] for ei,ej in etab]) == ecratio[c],
                    name = "NEdgeColors(%d)=%d" % (c,ecratio[c])
                )
        if set_ecratio: self.ecratio = ecratio
        return

    def setup_vcratio(self, vcratio, set_vcratio=True):
        """
        :Parameters:
        - vcratio (None or list of ints): overall vertex color ratio

        TBI: unset ratio with negative integers (convention: -1) and reserve
        the complement of the set elements
        """
        if not hasattr(self,"nvcolors"):
            self.setup_colors(nvcolors=len(vcratio))
        self.assert_vcratio(vcratio)
        vvars = self.vvars
        nvcolors = self.nvcolors
        ### ratio ###
        nvvars = len(vvars) / nvcolors
        crange = range(nvcolors)
        vcratio = self.normalize_ratio(vcratio, nvvars)
        ### loop ###
        rvvars = range(nvvars) ### == range(self._mol.natoms)
        for c in crange:
            self.model.addCons(
                quicksum(vvars[v,c] for v in rvvars) == vcratio[c],
                name = "NVertexColors(%d)=%d" % (c,vcratio[c])
            )
        if set_vcratio: self.vcratio = vcratio
        return

    def setup_ecratio_per_vertex(self, ecratio, set_ecratios=True):
        ### TBI: assert no conflict btw. global and local cratio
        if not hasattr(self,"necolors"):
            self.setup_colors(necolors=len(ecratio))
        self.assert_ecratio(ecratio)
        evars = self.evars
        necolors = self.necolors
        crange = range(necolors)
        vertex2edges = self.vertex2edges
        # loop #
        ecratios = []
        for c in crange:
            for v in vertex2edges:
                v2e = vertex2edges[v]
                necratio = self.normalize_ratio(ecratio,len(v2e))
                ecratios.append(necratio)
                consname = "NEdgeColorsPerVertex(%d,%d)=%d" % (v,c,necratio[c])
                self.model.addCons(
                    quicksum([evars[ei,ej,p,c] for ei,ej,p in v2e]) == necratio[c],
                    name = consname
                )
        if set_ecratios: self.ecratios = ecratios

    def setup_vcratio_per_edge(self, vcratio, set_vcratios=True):
        """
        N.B.: there are always two vertex per edge
        """
        if not hasattr(self,"nvcolors"):
            self.setup_colors(nvcolors=len(vcratio))
        self.assert_vcratio(vcratio)
        vvars = self.vvars
        nvcolors = self.nvcolors
        crange = range(nvcolors)
        edge2vertices = self.edge2vertices
        vcratios = []
        for c in crange:
            for e in edge2vertices:
                e2v = edge2vertices[e]
                nvcratio = self.normalize_ratio(vcratio,len(e2v))
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

    ############################################################################
    ### MAIN ###
    def optimize_model(self):
        self.model.optimize()
        if self.is_optimal_model():
            self.set_colors()
            if self.evars:
                self.set_otab(self.ecolors)
        return

    def free_model(self):
        self.model.freeTransform()
        return

    def is_optimal_model(self):
        return self.model.getStatus() == 'optimal'

    ############################################################################
    ### CYCLE ###

    def cycle_model(self, imax=1e2):
        """
        :TODO: URGENTLY!
        - separate vmodel from emodel is NEEDED here
        """
        try:
            i = 0
            while i < imax:
                self.report_step(i)
                self.optimize_model()
                if not self.is_optimal_model():
                    break
                self.free_model()
                if self.evars:
                    self.exclude_edge_solution()
                if self.vvars:
                    self.exclude_vertex_solution()
                i += 1
        except KeyboardInterrupt:
            i *= -1
        self.report_cycle_status(i)
        return
        """
        # init #
        self.optimizeReference() #???
        if sym:
            self.m.self.symperms = self.cycleInitSym(alpha)
        # loop #
        symMolAll = [] #???
        try:
            while i < imax:
                i += 1
            # self.cycleCheckMaxiter(i) ???
        except KeyboardInterrupt:
            logger.warning("\nINTERRUPTED: Convergence NOT reached!" + self.cycleStatusStr(i))
        ### logger.timer???
        # end #
        #print("Total elapsed time: %.3f s" % (self.model.getTotalTime(),))
        return symMolsAll #???
        """

    ############################################################################
    ### CONSTRAINTS ###

    def validate_edge_constraints(self, ecolors=None):
        if ecolors is None:
            ecolors = self.ecolors
        self.assert_ecolors(ecolors)
        evars = self.evars
        etab = self.etab
        if self._mol.use_pconn:
            for k,c in enumerate(ecolors):
                (ei,ej), p = etab[k]
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
        etab = self.etab
        if self._mol.use_pconn:
            for k,c in enumerate(ecolors):
                (ei,ej), p = etab[k]
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
        etab = self.etab
        if self._mol.use_pconn:
            self.model.addCons(
                quicksum([evars[e[0][0],e[0][1],e[1],ecolors[k]] for k,e in enumerate(etab)]) >= 1,
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
        etab = self.etab
        nevars = len(etab)
        if self._mol.use_pconn:
            self.model.addCons(
                quicksum([evars[e[0][0],e[0][1],e[1],ecolors[k]] for k,e in enumerate(etab)]) <= nevars - 1,
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
    ### SYMMETRY ###

    ############################################################################
    ### GET / SET ###

    ### N.B. otab and oconn are meant to be interfaces! They do not constraint
    ###     variables of the model
    def set_otab(self, otab, set_oconn=True):
        self.assert_otab(otab)
        self._mol.otab = self.otab = otab
        if set_oconn:
            self.set_oconn_from_otab(use_otab=False)

    def set_oconn(self, oconn, set_otab=True):
        self.assert_oconn(oconn)
        #self.oconn = oconn
        self._mol.oconn = self.oconn = oconn
        if set_otab:
            self.set_otab_from_oconn(use_oconn=False)

    def get_otab(self):
        return self.otab

    def get_oconn(self):
        return self.oconn

    def set_otab_from_oconn(self, oconn=None, use_oconn=False):
        #check/set oconn
        if oconn is None:
            oconn = self.oconn
        elif use_oconn:
            self.set_oconn(oconn)
        pass

    def set_oconn_from_otab(self, otab=None, use_otab=False):
        # WISHLIST: if pconn: FASTER! (thoughts: use array for same-lenght connectivity)
        #check/set otab
        if otab is None:
            otab = self.otab
        elif use_otab:
            self.set_otab(otab)
        # assign locally to be at hands
        # N.B.: DO NOT modify it! They are python lists
        conn = self._mol.conn
        ctab = self._mol.ctab
        if self._mol.use_pconn:
            ptab = self._mol.ptab
            pconn = self._mol.pconn
            ### image index connectivity (to avoid np.array comparison, arguably faster)
            iconn = [[arr2idx[j] for j in ic] for ic in pconn]         
        # init oconn
        oconn = [[None for j in ic] for ic in conn] # init as nested lists of Nones
        ### TBI: use etab! ###
        if self._mol.use_pconn: # periodic connectivity, ambiguity, slower and harder
            for k,(i,j) in enumerate(ctab):
                ### i -> j
                kp = ptab[k] ### image index
                ji_ = [j_ for j_, ji in enumerate( conn[i]) if ji == j ]
                jk_ = [j_ for j_, jk in enumerate(iconn[i]) if jk == kp]
                if ji_ == jk_:
                    j_ = ji_
                else:
                    j_ = set(ji_) & set(jk_) # intersect
                assert len(j_) == 1, "set of lenght %d, expected lenght is 1" % len(j_)
                j_ = list(j_)[0]
                oconn[i][j_] = otab[k]
                ### j -> i
                rp = idx2revidx[kp] ### reverse image index (by convention)
                ij_ = [i_ for i_, ij in enumerate( conn[j]) if ij == i ]
                ik_ = [i_ for i_, ik in enumerate(iconn[j]) if ik == rp]
                if ij_ == ik_:
                    i_ = ij_
                else:
                    i_ = set(ij_) & set(ik_) # intersect
                assert len(i_) == 1, "set of lenght %d, expected lenght is 1" % len(i_)
                i_ = list(i_)[0]
                oconn[j][i_] = otab[k]
        else: # no periodic connectivity, no ambiguity, faster and easier
            for k,(i,j) in enumerate(ctab):
                j_ = conn[i].index(j)
                oconn[i][j_] = otab[k]
                i_ = conn[j].index(i)
                oconn[j][i_] = otab[k]
        self.set_oconn(oconn, set_otab=False)
        return

    def set_edge_vars(self, necolors, set_necolors=True):
        ###TBI: test necolors = 0
        self.assert_ncolors_number(necolors)
        crange = range(necolors)
        etab = self.etab
        evars = {} # dict[(3or4)-tuple] = var... inefficient yet clear
        if self._mol.use_pconn: # periodic connectivity, ambiguity
            for (ei,ej),p in etab:
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
        ###TBI: test nvcolors = 0
        self.assert_ncolors_number(nvcolors)
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
            etab = self.etab
            if self._mol.use_pconn:
                for (ei,ej),p in etab:
                    try:
                        vals = [self.model.getVal(evars[ei,ej,p,c]) for c in crange]
                        val = vals.index(1)
                    except Warning:
                        val = "u"
                    ecolors.append(val)
            else:
                for ei,ej in etab:
                    try:
                        vals = [self.model.getVal(evars[ei,ej,c]) for c in crange]
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
    ### CHECK ###
    # WISHLIST: connectivity checker
    # - table against list
    # - i < j for ctab

    def check_():
        pass

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
    
    def set_etab(self):
        #THIS METHOD SHOULD BE EITHER IN CORE MOLSYS OR IN UTIL, DEF.LY NOT HERE
        ctab = self._mol.ctab
        if self._mol.use_pconn:
            ptab = self._mol.ptab
            etab = list(zip(ctab, ptab)) # python3 compl.: zip iterator gets exhausted
        else:
            etab = ctab
        self.etab = etab


    def get_etab(self):
        #THIS METHOD SHOULD BE EITHER IN CORE MOLSYS OR IN UTIL, DEF.LY NOT HERE
        return self.etab
    def normalize_ratio(self, cratio, total):
        """
        return normalized color ratio so that the total number of elements (edges 
        or vertices) is colored and the actual (float) ratio among colors is
        close to the given (int) ratio. It implements D'Hondt's quotients
        internally. Credits: https://github.com/rg3
        N.B.: in case of "ties", the latter color gets the higher ratio of
        elements. (this is consistent AND reproducible!)
        
        :Parameters:
         - cratio (list of ints): color ratio
         - total (int): total number of elements to color with given ratio
        :Returns:
         - norm_cratio (list of ints): color ratio normalized wrt. elements
        """
        def dhondt_quotient(cratio, subtotal):
            return float(cratio) / (subtotal + 1)
        
        # Calculate the quotients matrix (list in this case)
        dict_cratio = dict(enumerate(cratio))
        quot = []
        ret ={} 
        for p in dict_cratio:
            ret[p] = 0
            for s in range(0, total):
                q = dhondt_quotient(dict_cratio[p], s)
                quot.append((q, p))

        # Sort the quotients by value
        quot.sort(reverse=True)

        # Take the highest quotients with the assigned parties
        for s in range(0, total):
            ret[quot[s][1]] += 1
        norm_cratio = ret.values()
        return norm_cratio # already ordered

    def setup_vertex2edges(self, sort_flag=True):
        """
        setup dictionary from vertices to list of edges
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
        setup dictionary from edges to list of vertices
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
        etab = self.etab
        # loop #
        if self._mol.use_pconn:
            for (ei,ej),p in etab:
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
        setup dictionary from vertices to list of edge variables
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
        setup dictionary from edges to list of vertex variables
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
        etab = self.etab
        # loop #
        if self._mol.use_pconn:
            for (ei,ej),p in etab:
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
                    val = "!" #stands for unsolved
                except ValueError:
                    val = "?"
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
            etab = self.etab
            if self._mol.use_pconn:
                for (ei,ej),p in etab:
                    try:
                        vals = [self.model.getVal(evars[ei,ej,p,c]) for c in crange]
                        val = vals.index(1)
                    except Warning: #model may be unsolved
                        val = "?"
                    except ValueError: #1 is not found!
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

    def report_step(self, i):
        sys.stdout.write("cycle n.: %d\n" % i)

    def report_last_status(self, i):
        status = self.model.getStatus()
        logger.info("Last status: %s; Number of cycles: %s" % (status,i))
        return

    def report_cycle_status(self, i):
        if i < 0:
            self.report_last_status(-i)
            logger.warning("INTERRUPTED: Convergence NOT reached")
        else:
            self.report_last_status(i)
            status = self.model.getStatus()
            if status == "unknown":
                logger.warning("FAILURE: Convergence NOT reached")
            elif status == "timelimit":
                logger.warning("FAILURE: Time limit has been reached")
            elif status == "infeasible":
                if i == 0:
                    logger.info("FAILURE: No feasible solution found")
                else:
                    logger.info("SUCCESS: Convergence reached")
            else:
                logger.error("DISASTER: the unexpected happened!")
        return

    ############################################################################
    ### ASSERT ###
    
    def assert_flag(self,flag):
        assert isinstance(flag,bool), "flag must be boolean" \
            "you gave %s and it is %s" % (flag, flag.__class__.__name__)
        return
    
    def assert_otab(self,otab):
        assert len(otab) == len(self._mol.ctab)
        return
    
    def assert_oconn(self,oconn):
        assert len(oconn) == self._mol.natoms
        return
    
    def assert_hasmodel(self):
        assert hasattr(self, "model"), "model must be setup before this method"
        return
    
    def assert_ncolors_number(self, ncolors):
        assert ncolors >= 0,   "number of colors is %s but cannot be negative" % ncolors
        return

    def assert_nevcolors_number(self, necolors, nvcolors):
        self.assert_ncolors_number(necolors)
        self.assert_ncolors_number(nvcolors)
        assert necolors > 0 or nvcolors > 0, \
            "At least one btw. edge and vertex colors must be positive: " \
            "you gave %s and %s" % (necolors, nvcolors)
        return

    def assert_ecratio(self, ecratio):
        assert hasattr(ecratio, "__iter__"), "color ratio must be iterable, " \
            "you gave %s" % ecratio.__class__.__name__
        assert len(ecratio) == self.necolors, "given edge ratio %s is " \
            "inconsistent with maximum number of edge colors %s" % (ecratio, self.necolors)
        for weight in ecratio:
            assert isinstance(weight,Number), "color weight in color ratio " \
                "must be number, you gave %s and it is %s" % (weight, weight.__class__.__name__)
        return

    def assert_vcratio(self, vcratio):
        assert hasattr(vcratio, "__iter__"), "color ratio must be iterable, " \
            "you gave %s" % vcratio.__class__.__name__
        assert len(vcratio) == self.nvcolors, "given vertex ratio %s is " \
            "inconsistent with maximum number of vertex colors %s" % (vcratio, self.nvcolors)
        for weight in vcratio:
            assert isinstance(weight, Number), "color weight in color ratio " \
                "must be number, you gave %s and it is %s" % (weight, weight.__class__.__name__)
        return

    def assert_ecolors(self, ecolors):
        assert len(ecolors) == len(self.etab), "length of edge colors and " \
        "number of edges mismatch: %d vs %d" % (len(ecolors), len(self.etab))

    def assert_vcolors(self, vcolors):
        assert len(vcolors) == self._mol.natoms,"length of vertex colors and " \
        "number of vertex mismatch: %d vs %d" % (len(vcolors), self._mol.natoms)

    ############################################################################
    ### TEST ###
    ### notes [RA]:
    ### - if constraints are maintained as intended before/after optimization
    ###     before/after freeTransform -> CHECK SAME CONSTRAINTS
    ### - enable none/one/either/both of necolors and nvcolors -> OK, expected?
    ### - integer/fractional/zero ratios -> OK
    ### - give negative ratios -> LEAVE UNSET, KEEP THE RATIO OUT
    ### - None, [], else? for otab and oconn -> ERROR?
    ### - give both otab and oconn -> ERROR
    ### - give neither otab nor oconn -> OK
    ### - check if printing is right
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

__version__ = "2.0.0"

header= """
                                 .';ldxO00KKKK0Okdl:'.
                            ':lxO0kdoc:;;      ;:cldxO0koc,.
                       .;ok0xl;..            :xxdc.      :dK0.
                    'lOOo;. .cdd:            OKxd0Kc      cKl               ..
                 ,o0x:.     OKkkKO'         'KK,.kKc     ,Kd              .',.
              .l0k:.       .K0  lKK,        oKKKKk.     .00'            .,;'
            .d0d.          .KK'  dK0.       0KdkKKc     OKxkkl,       .,;,.
          .dKo. ,xd.        dKO.'lKKO      'KO  OK0    dKc,:loxko,  .,;,.
         cKx.  'KKk.        .OKK0xl0Ko     oKc  OK0   cKk:'..':lod:,,;,.
       'OO,    xKK.          .0Kl  '0Kc   .0KdcxK0,  ;K0ooooc,.';,,,;c.
      cKo      ,0Kd.          'KKc  .dc    .,:clc.  .0KNX0xolc;,'',;cxOl
     dK:         lK0l.     .   ,xd   :ldxkkOOkxdlc,'O0ok0XKo;,'',;,.':lkk'
    dK;           .oK0oclx0KO.  ,oOOdc;'         :ok0'  ,c;,'.',:loc..;cdOc
   oK;.,,,'.        .oOK0ko;.,d0d;.          ....      .,,..',:O0dooo,.'coOo
  cKl,KK0OKKKOxl;.         ;OO:      .',;,''......;..','...,,'x0NNkooo:..coOd
 .0k '0Kc     kKKK0kdc.  .x0:     ...';,,,..lo....;,,'...','   oOXWOooo:..coOl
 dK'  .dKkc   kKc   ld. ,0k.     .,;;;,...cOO:...;,'...',,.     cOXWOood;.'cdO,
.Kx     .;o0KKKKc      ;Kd      .','....,xO0O;..;'....,,.        cOXWkooo'.;ckk
lK,           ;lxOKKO.'Kx      ,;'.....dO0XX0Oollc::coddddl;.     d0NXoooc..cdO;
kKd                   xK.       ,;,..,kOKXXXXXKK0000KKKKXXK0k;    ,OKWkooo..:lOo
kKdliclkK0OO00OOkxdollKx      .,'...lO0XXXx..,cok0XXXXXXXXXXOk.   .O0W0ood,.;ckk
       cOoc..oddkNK0l;:,      ,,.,;xOKXXXXk   ,xl; lXXXXXXXXOk.    x0WKood;.,ckO
       cOoc..oookWXO,    ',''',;lkOKXXXXXXXo.  ,;,lXXXXXXXX0k,     k0WKood;.,ckO
       ;Odc..loodNNOc    'cc;cdO0KXXXXXXXXXXXOdxOXXXXXXXXKOd.     .OKW0ood'.;cOx
       .Oxc,.;dooKW0k    ,xOO0KXXXXXXXXXXXXXXXXXXKOkxOOOko,       :OXWxooo..coOc
        dOl:..oooxWXOc     ;xOKXXXXXXXXXXXXXXXXXKOc              .k0WKood;.,cxO.
        'Oxc,.;doo0WKO;      'dO0XXXXXXXXXXXXXXXOx              .x0NXdool..:lOo
         lOoc..:ooo0WKO:       .oO0XXXXXXXXXXXXKO;             .x0NNxooo..;ckk.
          dOl:..cooo0WX0o.       .lO0KXXXXXXXXXOx             ,kKNXdooo'.;cxO,
          .xkl:..:oookNN0k;        .ck0KXXXXXXOk.           .d0XN0oool'.;cxO;
           .dOoc'.,ooodONXKkc.       .:kOKXX0Oo.          ,d0XNKxooo:..;lxO,
             cOxc;..:ooodONNK0x:.       ,dOxc.        .,oOKXNKxoool'.'cokx.
              .xkoc,..:ooookKNNK0Oo:,..          .';lx0KXNXOdoool,..:lxOc
                :kxoc,..;looook0XNNXKK00kxxddxkO00KXXNXKOdoooo:'.';ldOd.
                  :kkoc;..';loooooxkO0KXXNNNNNNXXK0Oxdoooooc,..,:ldOd'
                    ,dkxoc;'..';clooooooooooooooooooool:,...,:ldkkc.
                      .;dOxdl:;'....',;::ccccc::;;,'...',:codkkc.
                         .,cxkxdolc:;,,'''....''',,;:lodxkko;.
                              .;ldkkkxxddddodddddxxkkxo:,.
                                   .,:ldxkOOOOkxoc;.
"""

footer="""Your ride with ACAB ends here. Houyhnhnm!"""

