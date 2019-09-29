# -*- coding: utf-8 -*-

"""

       module to implement an addon feature: graphs using the graph_tool library

       NOTE: this is only imported by __init__.py if graph_tool is present

"""

from graph_tool import Graph
from graph_tool.topology import *
import copy
import numpy as np
import molsys

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

    """decomposition

    Since the molg graph above is really used for fragment finding with a special way to treat 
    hydrogens we use a seperate graph called moldg (mol decomposition graph) which contains all atoms
    and uses only lower case element symbols as vertex symbols. 
    """

    def decompose(self, mode="ringsize"):
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
        else:
            print("unknown decomosition mode")
            return
        # now get the BBs and the net and collect the output
        net = self.get_net()
        return net

    def make_decomp_graph(self):
        """make a mol graph for decomposition
        """
        self.moldg = Graph(directed=False)
        self.moldg.vp.bb   = self.moldg.new_vertex_property("int")
        self.moldg.vp.con  = self.moldg.new_vertex_property("string")
        self.moldg.vp.filt = self.moldg.new_vertex_property("bool")
        # generate vertices and edges
        self.moldg.vp.elem = self.moldg.new_vertex_property("string")
        elems = self._mol.get_elems()
        conn  = self._mol.get_conn()
        for i in range(self._mol.natoms):
            self.moldg.add_vertex()
            self.moldg.vp.elem[i] = elems[i]
            for j in conn[i]:
                if j<i:
                    self.moldg.add_edge(self.moldg.vertex(i), self.moldg.vertex(j))
        # add all further properties
        self.moldg.ep.Nk   = self.moldg.new_edge_property("int")
        self.moldg.ep.filt = self.moldg.new_edge_property("bool")
        # init properties
        for v in self.moldg.vertices():
            self.moldg.vp.con[v] = ""
        self.moldg.ep.Nk.a[:] = 0
        for e in self.moldg.edges():
            self.moldg.ep.filt[e] = True
        # set up flags
        self.decomp_split = False
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
        self.plot_graph("topo", g=self.moldg, label="elem", method="sfdb")
        return

    def get_net(self, mode="coc"):
        """this method takes a sliced moldg and creates a graph of the underlying net

        TBI: currently no pconn is detected .. we rely on being able to construct it from the embedding
             this needs to be tested for very small thigs like a 1x1x1 pcu based MOF

        Args:
            mode (str, otional): Defaults to "coc", mode how to compute the position of the BBs vertex position
                                 coc - center of connectors
        """
        assert self.decomp_split is True
        label_components(self.moldg, vprop=self.moldg.vp.bb)
        self.moldg.set_edge_filter(None)
        # prepare lists with mappings for the vertices to the BBs (which atom is in which vertex BB etc)
        self.decomp_nv = self.moldg.vp.bb.a.max()+1
        self.decomp_net_conn = []
        self.decomp_map_atoms = []
        self.decomp_map_cons = []
        for i in range(self.decomp_nv):
            self.decomp_net_conn.append([])
            self.decomp_map_atoms.append([])
            self.decomp_map_cons.append([])
        # now check all split edges and register con and add to net connection table
        for e in self.moldg.edges():
            if self.moldg.ep.filt[e] == False:
                i = e.source()
                j = e.target()
                # also add the connector info as a string because we could have multiple
                self.moldg.vp.con[i] += "%d " % j
                self.moldg.vp.con[j] += "%d " % i
                ibb = self.moldg.vp.bb[i]
                jbb = self.moldg.vp.bb[j]
                self.decomp_net_conn[ibb].append(jbb)
                self.decomp_net_conn[jbb].append(ibb)
        # now compute the vertex positions of the BBs as the scaled embedding depending on the mode
        # first generate the mappings
        for v in self.moldg.vertices():
            vbb = self.moldg.vp.bb[v]
            self.decomp_map_atoms[vbb].append(int(v))
            if len(self.moldg.vp.con[v]) > 0:
                # this is a connector
                self.decomp_map_cons[vbb].append(int(v))
        # now compute the centers
        vpos = np.zeros([self.decomp_nv, 3], dtype="float64")
        for i in range(self.decomp_nv):
            if mode == "coc":
                # compute by centroid of connectors (note that this considers pbc correctly but we might have to wrap)
                vpos[i] = self._mol.get_coc(idx = self.decomp_map_cons[i])
            else:
                print("unknown mode in get_net")
                raise
        # now make a mol object
        self.decomp_net = molsys.mol.from_array(vpos)
        # since this is the unscaled embedding we can directly take the cell from the parent mol
        self.decomp_net.set_cell(self._mol.get_cell())
        # IMPROVE: this is the element map stolen from lqg.py --> lqg should evetnually go into a topo addon .. make more smart
        elems_map = {2:'x',3:'n',4:'s',5:'p',6:'o',8:'c'}
        vert_elems = []
        for i in range(self.decomp_nv):
            vert_elems.append(elems_map[len(self.decomp_net_conn[i])])
        self.decomp_net.set_elems(vert_elems)
        self.decomp_net.set_conn(self.decomp_net_conn)
        self.decomp_net.make_topo()
        return self.decomp_net



    def get_bbs(self):
        """this method takes a split moldg graph and generates the unique bbs
        """
        return
