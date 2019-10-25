# -*- coding: utf-8 -*-

"""

       module to implement an addon feature: graphs using the graph_tool library

       NOTE: this is only imported by __init__.py if graph_tool is present

"""

from graph_tool import Graph, GraphView
from graph_tool.topology import *
import copy
import numpy as np
import molsys
from molsys.util import elems

import logging
logger = logging.getLogger("molsys.graph")

class graph(object):

    def __init__(self, mol):
        """
        instantiate a graph object which will be attached to the parent mol

        :Parameter:

             - mol : a mol type object
        """
        self._mol = mol
        self.molg = Graph(directed=False)
        self.molg.vp.type = self.molg.new_vertex_property("string")
        self.molg.vp.molid = self.molg.new_vertex_property("int")
        # defaults
        self.moldg = None
        self.bbg   = None
        logger.debug("generated the graph addon")
        return

    def make_graph(self, idx = None, hashes = True):
        """
        generate a graph for the mol object (atoms should be typed)
        we use the atomtype name with the "_" and everything after it (rule=2) truncated.
        in other words the vertex property is the element plus the coordination number

        """
        if idx is None: idx = range(self._mol.natoms)
        # now add vertices
        self.vert2atom = [] # this list maps vertex indices to the real atoms becasue we omit the hydrogens in the graph
        ig = 0
        for i in idx:
            if self._mol.elems[i] != "x":
                self.molg.add_vertex()
                self.vert2atom.append(i)
                vtype = self._mol.atypes[i]
                # extract element and coordination number
                if "_" in vtype:
                    vtype = vtype.split("_")[0]
                # if the coordination number is one replace the element by a #
                if hashes:
                    if vtype[-1] == "1":
                        vtype = "#"
                self.molg.vp.type[ig] = vtype
                ig += 1
        self.nvertices = len(self.vert2atom)
        logger.info("generated a graph for a mol object with %d vertices" % self.nvertices)
        # now add edges ... only bonds between vertices
        for i in range(self.nvertices):
            ia = self.vert2atom[i]
            for ja in self._mol.conn[ia]:
                if ja>=ia:   #we need a .le. here for those atoms/vertices connected to itself twice in different boxes
                    if ja in self.vert2atom:
                    # print("bond from %d to %d" % (ia, ja))
                    # print(self._mol.atypes[ia], self._mol.atypes[ja])
                        self.molg.add_edge(self.molg.vertex(i), self.molg.vertex(self.vert2atom.index(ja)))
                        #self.molg.add_edge( self.molg.vertex(self.vert2atom.index(ja)),self.molg.vertex(i))
        return

    def plot_graph(self, fname, g = None, label=None, edge_label=None, size=1000, fsize=16, vsize=8, ptype = "pdf",method='arf'):
        """
        plot the graph (needs more tuning options) [via **kwargs? RA]

        :Parameter:
            - fname  : filename (will write filename.pdf)
            - size   : outputsize will be (size, size) in px [default 800]
            - fsize  : font size [default 10]
            - method : placement method to draw graph, can be one of
                       arf
                       frucht
                       radtree
                       sfdp
                       random
        """
        import graph_tool.draw as gtd
        if g:
            draw_g = g
        else:
            draw_g = self.molg
        if label:
            vlabel = label
        else:
            vlabel = "type"
        g=draw_g
        if method=='arf':
            pos = gtd.arf_layout(draw_g, max_iter=0)
        elif method=='frucht':
            pos = gtd.fruchterman_reingold_layout(draw_g, n_iter=1000)
        elif method=='radtree':
            pos = gtd.radial_tree_layout(draw_g, draw_g.vertex(0))
        elif method=='sfdp':
            pos = gtd.sfdp_layout(draw_g)
        elif method=='random':
            pos = gtd.random_layout(draw_g)
        else:
            pos=None
        gtd.graph_draw(draw_g,pos=pos, vertex_text=draw_g.vp[vlabel], vertex_font_size=fsize, vertex_size=vsize, \
            output_size=(size, size), output=fname+"."+ptype, bg_color=[1,1,1,1])
        return

    @staticmethod
    def find_subgraph(graph, subg):
        """
        use graph_tools subgraph_isomorphism tool to find substructures

        :Parameter:

            - graph : parent graph to be searched
            - subg  : graph to be found

        :Returns:

            a list of lists with the (sorted) vertex indices of the substructure
        """
        maps = subgraph_isomorphism(subg, graph, vertex_label=(subg.vp.type, graph.vp.type))
        subs = []
        subs_check = []
        for m in maps:
            sl = list(m)
            sl_check = copy.deepcopy(sl)
            sl_check.sort()
            if sl_check not in subs_check: 
                subs.append(sl)
                subs_check.append(sl_check)
        return subs


    def find_sub(self, subg):
        """
        use graph_tools subgraph_isomorphism tool to find substructures

        :Parameter:

            - subg : graph object (from another molsys) to be searched

        :Returns:

            a list of lists with the (sorted) vertex indices of the substructure
        """
        subs = self.find_subgraph(self.molg, subg.molg)
        return subs

    def find_fragment(self, frag,add_hydrogen=False):
        """
        find a complete fragment (including the hydrogen atoms not included in the graph)
        Note that the fragment found can be different from the fragment by the number of hydrogen atoms!!

        :Parameter:

            - frag : mol object with graph addon to be found

        :Returns:

            a list of lists with the atom indices of the fragment in the full system
        """
        subs = self.find_sub(frag.graph)
        frags = []
        for s in subs:
            # loop over all vertices
            f = []
            for v in s:
                a = self.vert2atom[v]
                f.append(a)
                # check all atoms connected to this atom if they are hydrogen
                if add_hydrogen:
                    for ca in self._mol.conn[a]:
                        if self._mol.elems[ca] == "h":
                            f.append(ca)
            frags.append(f)
        return frags

    def util_graph(self, vertices, conn):
        """
        generate a graph with vertices and connectivity in conn
        """
        g = Graph(directed=False)
        # now add vertices
        g.vp.type = g.new_vertex_property("string")
        for i, v in enumerate(vertices):
            g.add_vertex()
            g.vp.type[i] = v
        # now add edges ...
        for i, v in enumerate(vertices):
            for j in conn[i]:
                if j>=i:
                    g.add_edge(g.vertex(i), g.vertex(j))
        return g

    def filter_graph(self, idx):
        """
        filters all atoms besides the given out of the graph
        :Parameters:
            - idx (list): indices of atoms to keep
        """
        # TODO use vert2atom
        assert type(idx) == list
        self.molg.clear_filters()
        filter = self.molg.new_vertex_property("bool")
        filter.set_value(False)
        for i in idx:
            filter[self.molg.vertex(i)]=True
        self.molg.set_vertex_filter(filter)
        return

    def get_components(self):
        """Get all the components aka molecules from the atomic graph

        it adds a property map molid to the graph

        """
        label_components(self.molg, vprop=self.molg.vp.molid)
        return

    """
    #############################################################################################
           decomposition (or deconstruction)
    #############################################################################################

    Since the molg graph above is really used for fragment finding with a special way to treat 
    hydrogens we use a seperate graph called moldg (mol decomposition graph) which contains all atoms
    and uses only lower case element symbols as vertex symbols. 
    """

    def decompose(self, mode="ringsize", join_organic=True):
        """decompose the molecular graph into BBs and the underlying net

        There are different strategies to split the graph into parts and also the returned information
        can be tuned.
        In general, a topo mol object of the underlying net and a set of BB mol objects with a mapping 
        to vertices and edges will be returned.
        
        modes:
            ringsize: use clusters of minimum ring size as used by Topos for decomposition

        Args:
            mode (str, optional): mode of splitting into BBs. Defaults to "ringsize".
        """
        # set up the graph
        self.make_decomp_graph()
        # split it
        if mode == "ringsize":
            self.split_ringsize()
            self.make_bb_graph()
            # if we join organic BBs do this now and regenerate bbg afterwards
            if join_organic:
                if self.join_organic():
                    self.make_bb_graph()
            # now test if there have been 2c BBs side by side and remove these splits
            if self.join_2c():
                self.make_bb_graph()
        else:
            print("unknown decomosition mode")
            return
        # now get the BBs and the net and collect the output
        net = self.get_net()
        # generate the BBs as mol objects and theri distribution to the vertices and edges
        #   NOTE: return values are also available as self.decomp_<name>
        vbb, vbb_map, ebb, ebb_map = self.get_bbs()
        # return complete info as a tuple
        return (net, vbb, vbb_map, ebb, ebb_map)

    def make_decomp_graph(self):
        """make a mol graph for decomposition
        """
        self.moldg = Graph(directed=False)
        self.moldg.vp.bb   = self.moldg.new_vertex_property("int")     # bb index
        self.moldg.vp.con  = self.moldg.new_vertex_property("string")  # connector atom (-1 no con, else atom to which is bonded)
        self.moldg.vp.filt = self.moldg.new_vertex_property("bool")    # vertex filter (for endgroups aka hydrogen)
        self.moldg.vp.atype= self.moldg.new_vertex_property("string")  # atomtype
        # generate vertices and edges
        self.moldg.vp.elem = self.moldg.new_vertex_property("string")  # element
        elems = self._mol.get_elems()
        conn  = self._mol.get_conn()
        for i in range(self._mol.natoms):
            self.moldg.add_vertex()
            self.moldg.vp.elem[i] = elems[i]
            for j in conn[i]:
                if j<i:
                    self.moldg.add_edge(self.moldg.vertex(i), self.moldg.vertex(j))
        # add all further properties
        self.moldg.ep.Nk   = self.moldg.new_edge_property("int")      # smallest ringsize
        self.moldg.ep.filt = self.moldg.new_edge_property("bool")     # filter for splitting into BBs
        # init properties
        for v in self.moldg.vertices():
            self.moldg.vp.con[v] = ""
            # determine the atype from the graph
            atype = self.moldg.vp.elem[v]
            atype += "%d_" % v.out_degree()
            nel = [self.moldg.vp.elem[nv] for nv in v.all_neighbors()]
            nels = list(set(nel))
            nels.sort()
            for e in nels:
                atype += "%s%d" % (e, nel.count(e))
            self.moldg.vp.atype[v] = atype
        self.moldg.ep.Nk.a[:] = 0
        for e in self.moldg.edges():
            self.moldg.ep.filt[e] = True
        # set up flags
        self.decomp_split = False
        self.decomp_bbg = False
        self.decomp_bb_exist = False
        return

    def split_ringsize(self):
        """
        split by determining ringsizes follwoing the topos strategy
        """
        assert self.moldg is not None
        assert self.decomp_split is False
        # first step: remove all vertices of degree 1 (hydrogen)
        k = kcore_decomposition(self.moldg).get_array()
        k1 = np.not_equal(k, 1)
        self.moldg.vp.filt.a = k1
        self.moldg.set_vertex_filter(self.moldg.vp.filt)
        # second step: determine minimal ring size to each edge
        self.moldg.set_edge_filter(self.moldg.ep.filt)
        for e in self.moldg.edges():
            self.moldg.ep.filt[e] = False
            dist = shortest_distance(self.moldg, source=e.source(), target= e.target())
            self.moldg.ep.filt[e] = True
            if dist < 2147483647:
                self.moldg.ep.Nk[e] = dist+1
            else:
                self.moldg.ep.Nk[e] = 0
        self.moldg.set_edge_filter(self.moldg.ep.filt)
        # third step: determine where to cut
        Nks = list(set(self.moldg.ep.Nk.get_array()))
        Nks.sort()
        Nks.remove(0)
        # print Nks
        cut_at_Nk = []
        for i in range(len(Nks)-1):
            if Nks[i+1]-Nks[i]>2:
                #cut_at_Nk.append(Nks[i+1])
                cut_at_Nk = Nks[i+1:]
                break
        # print cut_at_Nk
        # fourth step: get all vertices back and set the cuts by filtering the edges
        self.moldg.set_vertex_filter(None)
        for e in self.moldg.edges():
            if self.moldg.ep.Nk[e] in cut_at_Nk:
                # this is a cut -> set edge filter to False
                self.moldg.ep.filt[e] = False
                # print "cut at bond %d %d" % (int(e.source()), int(e.target()))
        self.decomp_split = True
        # DEBUG DEBUG DEBUG
        # self.plot_graph("topo", g=self.moldg, label="elem", method="sfdb")
        return

    def make_bb_graph(self, plot=False):
        assert self.decomp_split is True
        # label the BBs -> this generates the BBid valid for both vertices and edges
        self.moldg.set_edge_filter(self.moldg.ep.filt)
        label_components(self.moldg, vprop=self.moldg.vp.bb)
        self.moldg.set_edge_filter(None)
        # prepare lists with mappings for the vertices to the BBs (which atom is in which vertex BB etc)
        # Note:
        #    the number of vertices nv is not equal to the number of BBs (nbb) since some BBs sit on edges
        #    we thus need a mapping from BBs to atom indices and another mapping of vertices to BBs
        self.decomp_nbb = self.moldg.vp.bb.a.max()+1
        # now we make a temporary graph bbg for the BBs
        self.bbg = Graph(directed=False)
        self.bbg.ep.satom = self.bbg.new_edge_property("int") #  these two edge properties store the connecting atom indices 
        self.bbg.ep.tatom = self.bbg.new_edge_property("int") #  from the original molg 
        for b in range(self.decomp_nbb):
            self.bbg.add_vertex()  # vertices of this graph can be indexed by bbid (it contains also the 2c edges)
        # now check all split edges, register the connected atoms and add edges to the bbg graph
        # reinit the con vertex property og moldg in case we do this again
        for v in self.moldg.vertices():
            self.moldg.vp.con[v] = ""
        for e in self.moldg.edges():
            if self.moldg.ep.filt[e] == False:
                i = e.source() 
                j = e.target()
                # also add the connector info as a string because we could have multiple
                if len(self.moldg.vp.con[i]) > 0:
                    print ("atom %d has already con %s ... connects to %d" % (i, self.moldg.vp.con[i], j)) 
                if len(self.moldg.vp.con[j]) > 0:
                    print ("atom %d has already con %s ... connects to %d" % (j, self.moldg.vp.con[j], i)) 
                self.moldg.vp.con[i] += "%d " % j
                self.moldg.vp.con[j] += "%d " % i
                ibb = self.moldg.vp.bb[i]
                jbb = self.moldg.vp.bb[j]
                # print ("bond %d-%d for bbs %d %d" % (i, j, ibb, jbb))
                # TBI if we have two bonds between BBs we get edges twice here ... we could check if the edge is already there
                ne = self.bbg.add_edge(self.bbg.vertex(ibb),self.bbg.vertex(jbb))
                # we define source(small bb) - > target(large bb)
                if ibb<jbb:
                    self.bbg.ep.satom[ne] = i
                    self.bbg.ep.tatom[ne] = j
                else:
                    self.bbg.ep.satom[ne] = j
                    self.bbg.ep.tatom[ne] = i
        # add vertex and edge properties for bbs and init them (edges are -1 by default => no bb)
        self.bbg.vp.bb = self.bbg.new_vertex_property("int")
        self.bbg.ep.bb = self.bbg.new_edge_property("int")
        self.bbg.vp.bb.a = np.arange(self.bbg.num_vertices())
        self.bbg.ep.bb.a = -1
        self.decomp_map_bb2atoms = []
        self.decomp_map_bb2cons = []
        for i in range(self.decomp_nbb):
            self.decomp_map_bb2atoms.append([])
            self.decomp_map_bb2cons.append([])
        # iterate over all vertices of the moldg 
        for v in self.moldg.vertices():
            vbb = self.moldg.vp.bb[v]
            self.decomp_map_bb2atoms[vbb].append(int(v))
            if len(self.moldg.vp.con[v]) > 0:
                # this is a connector
                self.decomp_map_bb2cons[vbb].append(int(v))
        # set flag that we can proceed
        self.decomp_bbg =True
        if plot:
            self.plot_graph("bbg", g=self.bbg, label="bb")
        return

    def join_2c(self):
        """this method combines all 2c BBs into one, needs bbg but a call to make_bbg_graph has to be repeated

        returns a flag which is True if a merge happend
        """
        assert self.decomp_bbg
        # iterate over edges in the molgraph ... those which are edges of the bbg -> check if both ends are 2c and remove
        bb2v = list(self.bbg.vp.bb.a)
        merge = False
        for e in self.moldg.edges():
            if self.moldg.ep.filt[e] == False:
                ibb = self.moldg.vp.bb[e.source()]
                jbb = self.moldg.vp.bb[e.target()]
                iv = bb2v.index(ibb)
                jv = bb2v.index(jbb)
                if (self.bbg.vertex(iv).out_degree() == 2) and (self.bbg.vertex(jv).out_degree() == 2):
                    # print ("DEBUG: vertices %d and %d are both 2c: merging" % (ibb, jbb))
                    self.moldg.ep.filt[e] = True
                    merge = True
        # now invalidate bbg
        if merge: 
            self.decomp_bbg = False
        return merge

    def join_organic(self):
        """this method combines organic BBs into one (be it edges or vertices)

        returns a flag which is True if a merge happened and make_bbg_graph needs to be rerun
        """
        assert self.decomp_bbg
        # iterate over vertices -> flag all vertices as organic or inorganic
        orgbbs = []
        molelem = self._mol.get_elems()
        for v in self.bbg.vertices():
            vbb = self.bbg.vp.bb[v]
            velem = [molelem[i] for i in self.decomp_map_bb2atoms[vbb]]
            organic = True
            for e in velem:
                if e in elems.metals:
                    organic = False
            if organic:
                orgbbs.append(vbb)
        # iterate over edges in the molgraph ... those which bridge two organic bbs will be cleared
        bb2v = list(self.bbg.vp.bb.a)
        merge = False
        for e in self.moldg.edges():
            if self.moldg.ep.filt[e] == False:
                ibb = self.moldg.vp.bb[e.source()]
                jbb = self.moldg.vp.bb[e.target()]
                if (ibb in orgbbs) and (jbb in orgbbs):
                    # print ("DEBUG: BBs %d and %d are both organic: merging" % (ibb, jbb))
                    self.moldg.ep.filt[e] = True
                    merge = True
        # now invalidate bbg
        if merge: 
            self.decomp_bbg = False
        return merge


    def get_net(self, mode="coc", plot=False):
        """this method takes a sliced moldg and creates a graph of the underlying net

        TBI: currently no pconn is detected .. we rely on being able to construct it from the embedding
             this needs to be tested for very small thigs like a 1x1x1 pcu based MOF

        Args:
            mode (str, otional): Defaults to "coc", mode how to compute the position of the BBs vertex position
                                 coc - center of connectors
        """        
        # now reduce the graph to its basis without the 2c edges (edge bbs go into the ep.bb)
        remove_v = []
        for v in self.bbg.vertices():
            if v.out_degree() == 2:
                remove_v.append(v)
                i, j = v.all_neighbors()
                new_edge = self.bbg.add_edge(self.bbg.vertex(i), self.bbg.vertex(j))
                # we have to assign the satom/tatom for this new edge from the (to be destroyed) edges i-v-j
                # using the order i smaller than j (for the bb index)
                iv = self.bbg.edge(i,v)
                jv = self.bbg.edge(j,v)
                if i<v:
                    ina = self.bbg.ep.satom[iv]
                else:
                    ina = self.bbg.ep.tatom[iv]
                if j<v:
                    jna = self.bbg.ep.satom[jv]
                else:
                    jna = self.bbg.ep.tatom[jv]
                if i<j:
                    self.bbg.ep.satom[new_edge] = ina
                    self.bbg.ep.tatom[new_edge] = jna
                else:
                    self.bbg.ep.tatom[new_edge] = ina
                    self.bbg.ep.satom[new_edge] = jna
                self.bbg.ep.bb[new_edge] = self.bbg.vp.bb[v]
                #print ("REMOVE edge BB %d - connected to BBs %d %d - edge atoms %d %d" % (int(v), i, j, self.bbg.ep.satom[new_edge], self.bbg.ep.tatom[new_edge]))
        # self.bbg.remove_vertex(remove_v)
        for v in reversed(sorted(remove_v)):
            self.bbg.remove_vertex(v)
        self.decomp_nv = self.bbg.num_vertices()
        if plot:
            self.plot_graph("bbg_after", g=self.bbg, label="bb")
        # generate xyz coordinates of the original system where the BBs are properly wrapped into one image
        wrap_xyz = self._mol.get_xyz().copy()
        for i in range(self.decomp_nbb):
            bb_xyz = wrap_xyz[self.decomp_map_bb2atoms[i]]
            bb_xyz = self._mol.apply_pbc(xyz=bb_xyz)
            wrap_xyz[self.decomp_map_bb2atoms[i]] = bb_xyz
        self.wrap_mol = self._mol.clone()
        self.wrap_mol.set_xyz(wrap_xyz)
        #DEBUG
        #self.wrap_mol.write("DEBUG_wrap_mol.mfpx")
        # now compute the vertex positions of the BBs as the scaled embedding depending on the mode
        vpos = np.zeros([self.decomp_nv, 3], dtype="float64")
        for v in self.bbg.vertices():
            if mode == "coc":
                # compute by centroid of connectors (note that this considers pbc correctly but we might have to wrap)
                bbid = self.bbg.vp.bb[v]
                # vpos[int(v)] = self._mol.get_coc(idx = self.decomp_map_bb2cons[bbid])
                vpos[int(v)] = self.wrap_mol.get_coc(idx = self.decomp_map_bb2cons[bbid])
            else:
                print("unknown mode in get_net")
                raise
        # now make a mol object
        self.decomp_net = molsys.mol.from_array(vpos)
        # since this is the unscaled embedding we can directly take the cell from the parent mol
        self.decomp_net.set_cell(self._mol.get_cell())
        # self.decomp_net.apply_pbc()
        # now generate a net_conn from the edges of the bbg graph
        self.decomp_net_conn = []
        self.decomp_net_pconn = []
        for i in range(self.decomp_nv):
            self.decomp_net_conn.append([])
            self.decomp_net_pconn.append([])
        for e in self.bbg.edges():
            i = int(e.source())
            j = int(e.target())
            self.decomp_net_conn[i].append(j)
            self.decomp_net_conn[j].append(i)
            ia = self.bbg.ep.satom[e] # atom indices of the original molg
            ja = self.bbg.ep.tatom[e] # that form this BB connection
            # print ("bond between BBs atom: %3d %3d   BBs: %2d %2d" % (ia, ja, i, j))
            dist_ij = wrap_xyz[ja] - wrap_xyz[ia]
            fracdist_ij = np.dot(dist_ij, self._mol.inv_cell)
            pconn = np.around(fracdist_ij)
            if i<j:
                self.decomp_net_pconn[i].append(-pconn)
                self.decomp_net_pconn[j].append(pconn)
            else:
                self.decomp_net_pconn[i].append(pconn)
                self.decomp_net_pconn[j].append(-pconn)
        vert_elems = []
        for i in range(self.decomp_nv):
            vert_elems.append(elems.topotypes[len(self.decomp_net_conn[i])])
        self.decomp_net.set_elems(vert_elems)
        self.decomp_net.set_conn(self.decomp_net_conn)
        # this is a bit hacky ... we want to use pconn but can not rely on check_need_pconn
        self.decomp_net.set_pconn(self.decomp_net_pconn)
        self.decomp_net.use_pconn = True
        self.decomp_net.make_topo(check_flag=False)
        return self.decomp_net

    def get_bbs(self, get_all=False, write_mfpx=True):
        """this method generates the unique bbs from moldg and bbg
        """
        assert self.bbg is not None
        if get_all == True:
            print "PLEASE IMPLEMENT get_all option!!"
            raise
        self.decomp_vbb = [] # vertex BB mol objects (only unique)
        self.decomp_ebb = [] # edge/linker BB mol objects (only unique)
        self.decomp_vbb_map = [] # mapping of BBs to the vertices
        self.decomp_ebb_map = [] # mapping of the BBS to the edges
        # VERTICES
        # start with the vertices ... check all
        vbb_graphs = []
        vbb_map   = {}
        vbb_remap = {}
        nvbbs = 0
        for v in self.bbg.vertices():
            bb = self.bbg.vp.bb[v]
            cur_bbsg = GraphView(self.moldg, directed=False, vfilt = self.moldg.vp.bb.a == bb)
            # now check if we have this bb already
            known = False
            for i in range(len(vbb_graphs)):
                old_bbsg = vbb_graphs[i]
                if old_bbsg.num_vertices() != cur_bbsg.num_vertices():
                    continue
                # check if isomorphic
                # if isomorphism(cur_bbsg, old_bbsg, vertex_inv1 = cur_bbsg.vp.elem, vertex_inv2 = old_bbsg.vp.elem):
                if isomorphism(cur_bbsg, old_bbsg):
                    known = True
                    break
            if not known:
                # add graph
                vbb_graphs.append(Graph(cur_bbsg, prune=True))
                vbb_map[nvbbs]   = [bb]
                vbb_remap[nvbbs] = [int(v)]
                nvbbs += 1
            else:
                # known already .. i equals nvbb
                vbb_map[i].append(bb)
                vbb_remap[i].append(int(v))
        # we assume that all the BBs have a similar structure and we pick the first from the list to generate the BB object
        # TBI: with flag get_all==True convert them all
        for i in range(nvbbs):
            bb   = vbb_map[i][0] # take the first
            vbbg = vbb_graphs[i]
            # now we need to convert bb into a mol object and store it in self.decomp_vbb
            mol_bb = self._convert_bb2mol(bb, vbbg)
            # now add to the final list
            self.decomp_vbb.append(mol_bb)
            # store the map from the dicitonary into the nested list
            self.decomp_vbb_map.append(vbb_remap[i])
            # DEBUG
            if write_mfpx:
                mol_bb.write("bb_%d.mfpx" % i)
        # EDGES
        # now analyze the edges .. check all edges if their edge property bb is larger than -1 (-1 means no edge bb) 
        ebb_graphs = []
        ebb_map = {}
        ebb_remap = {}
        nebbs = 0
        for e in self.bbg.edges():
            if self.bbg.ep.bb[e] > -1:
                bb = self.bbg.ep.bb[e]
                cur_bbsg = GraphView(self.moldg, directed=False, vfilt = self.moldg.vp.bb.a == bb)
                # now check if we have this bb already
                known = False
                for i in range(len(ebb_graphs)):
                    old_bbsg = ebb_graphs[i]
                    if old_bbsg.num_vertices() != cur_bbsg.num_vertices():
                        continue
                    # check if isomorphic
                    # if isomorphism(cur_bbsg, old_bbsg, vertex_inv1 = cur_bbsg.vp.elem, vertex_inv2 = old_bbsg.vp.elem):
                    if isomorphism(cur_bbsg, old_bbsg):
                        known = True
                        break
                if not known:
                    # add graph
                    ebb_graphs.append(Graph(cur_bbsg, prune=True))
                    ebb_map[nebbs] = [bb]
                    ebb_remap[nebbs] = [e]
                    nebbs += 1
                else:
                    # known already .. i equals nvbb
                    ebb_map[i].append(bb)
                    ebb_remap[i].append(e)
        # we assume that all the BBs have a similar structure and we pick the first from the list to generate the BB object
        # TBI: with flag get_all==True convert them all
        for i in range(nebbs):
            bb   = ebb_map[i][0] # take the first
            ebbg = ebb_graphs[i]
            # now we need to convert bb into a mol object and store it in self.decomp_ebb
            mol_bb = self._convert_bb2mol(bb, ebbg)
            # now add to the final list
            self.decomp_ebb.append(mol_bb)
            # store the map from the dicitonary into the nested list
            self.decomp_ebb_map.append(ebb_remap[i])
            # DEBUG
            if write_mfpx:
                mol_bb.write("ebb_%d.mfpx" % i)
        self.decomp_bb_exist = True
        return (self.decomp_vbb, self.decomp_vbb_map, self.decomp_ebb, self.decomp_ebb_map)

    def _convert_bb2mol(self, bb, bbg):
        """helper function to convert a bb (index) and bbg (graph) to a mol object
        
        Args:
            bb (int): index of building block to convert
            bbg (graph): graph of the bb

        Return:
            mol object
        """
        bb_atoms = self.decomp_map_bb2atoms[bb]
        bb_cons   = self.decomp_map_bb2cons[bb]
        xyz = self.wrap_mol.get_xyz(bb_atoms)
        #xyz = self._mol.get_xyz(bb_atoms)
        # generate a mol object just from positions (inherit cell)
        mol_bb = molsys.mol.from_array(xyz)
        cell = self._mol.get_cell()
        mol_bb.set_cell(cell)
        mol_bb.set_elems([bbg.vp.elem[v] for v in bbg.vertices()])
        mol_bb.set_atypes([bbg.vp.atype[v] for v in bbg.vertices()])
        mol_bb.set_real_mass()
        ctab = [[int(e.source()), int(e.target())] for e in bbg.edges()]
        mol_bb.set_conn_from_tab(ctab)
        # the positions could be split by pbc -> fix it and then remove periodicity
        mol_bb.center_com(check_periodic=False)
        mol_bb.apply_pbc()
        mol_bb.make_nonperiodic()
        mol_bb.addon("bb")
        # compute connector indices of local BB
        connector = [bb_atoms.index(j) for j in bb_cons]
        # compute atypes of atoms bonded to the connectors
        connector_atypes = [[self.moldg.vp.atype[int(self.moldg.vp.con[j])]] for j in bb_cons]
        mol_bb.bb.setup(connector, connector_atypes=connector_atypes)
        return mol_bb
