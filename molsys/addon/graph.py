# -*- coding: utf-8 -*-

"""

       module to implement an addon feature: graphs using the graph_tool library

       NOTE: this is only imported by __init__.py if graph_tool is present

"""

from graph_tool import Graph
from graph_tool.topology import *
import copy

import logging
logger = logging.getLogger("molsys.graph")

class graph:

    def __init__(self, mol):
        """
        instantiate a graph object which will be attached to the parent mol

        :Parameter:

             - mol : a mol type object (can be a derived type like bb or topo as well)
        """
        self._mol = mol
        logger.debug("generated the graph addon")
        return

    def make_graph(self, idx = None, hashes = True):
        """
        generate a graph for the mol object (atoms should be typed)
        we use the atomtype name with the "_" and everything after it (rule=2) truncated.
        in other words the vertex property is the element plus the coordination number

        """
        if idx == None: idx = range(self._mol.natoms)
        self.molg = Graph(directed=False)
        # now add vertices
        self.molg.vp.type = self.molg.new_vertex_property("string")
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

    def plot_graph(self, fname, g = None, size=1000, fsize=16, vsize=8, ptype = "pdf",method='arf'):
        """
        plot the grap (needs more tuning options

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
        if g:
            draw_g = g
        else:
            draw_g = self.molg
        import graph_tool.draw
        import graph_tool.draw as gt
        g=draw_g
        if method=='arf':
            pos = graph_tool.draw.arf_layout(draw_g, max_iter=0)
        elif method=='frucht':
            pos = graph_tool.draw.fruchterman_reingold_layout(draw_g, n_iter=1000)
        elif method=='radtree':
            pos = gt.radial_tree_layout(g, g.vertex(0))
        elif method=='sfdp':
            pos = gt.sfdp_layout(g)
        elif method=='sfdp':
            pos = gt.random_layout(g)
        else:
            pos=None
        from graph_tool.draw import graph_draw
        graph_draw(draw_g,pos=pos, vertex_text=draw_g.vp.type, vertex_font_size=fsize, vertex_size=vsize, \
            output_size=(size, size), output=fname+"."+ptype, bg_color=[1,1,1,1])
        return

    def find_subgraph(self, graph, subg):
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
        generate a generate a graph with vertices and connectivity in conn
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




