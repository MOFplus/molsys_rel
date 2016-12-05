# -*- coding: utf-8 -*-

"""

       module to implement an addon feature: graphs using the graph_tool library

       NOTE: this is only imported by __init__.py if graph_tool is present

"""

from graph_tool import Graph
from graph_tool.topology import *

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

    def make_graph(self):
        """
        generate a graph for the mol object (atoms should be typed)
        we use the atomtype name with the "_" and everything after it (rule=2) truncated.
        in other words the vertex property is the element plus the coordination number

        """
        self.molg = Graph(directed=False)
        # now add vertices
        self.molg.vp.type = self.molg.new_vertex_property("string")
        self.vert2atom = [] # this list maps vertex indices to the real atoms becasue we omit the hydrogens in the graph
        ig = 0
        for i in xrange(self._mol.natoms):
            if self._mol.elems[i] != "x":
                self.molg.add_vertex()
                self.vert2atom.append(i)
                vtype = self._mol.atypes[i]
                # extract element and coordination number
                if "_" in vtype:
                    vtype = vtype.split("_")[0]
                # if the coordination number is one replace the element by a #
                if vtype[-1] == "1":
                    vtype = "#"
                self.molg.vp.type[ig] = vtype
                ig += 1
        self.nvertices = len(self.vert2atom)
        logger.info("generated a graph for a mol object with %d vertices" % self.nvertices)
        # now add edges ... only bonds between vertices
        for i in xrange(self.nvertices):
            ia = self.vert2atom[i]
            for ja in self._mol.conn[ia]:
                if ja>=ia:   #we need a .le. here for those atoms/vertices connected to itself twice in different boxes
                if ja in self.vert2atom:
                    # print "bond from %d to %d" % (ia, ja)
                    # print self._mol.atypes[ia], self._mol.atypes[ja]
                    self.molg.add_edge(self.molg.vertex(i), self.molg.vertex(self.vert2atom.index(ja)))
        return

    def plot_graph(self, fname, size=800, fsize=10):
        """
        plot the grap (needs more tuning options

        :Parameter:
            - fname : filename (will write filename.pdf)
            - size  : outputsize will be (size, size) in px [default 800]
            - fsize : font size [default 10]

        """
        from graph_tool.draw import graph_draw
        graph_draw(self.molg, vertex_text=self.molg.vp.type, vertex_font_size=10,  \
            output_size=(size, size), output=fname+".pdf")
        return

    def find_sub(self, subg):
        """
        use graph_tools subgraph_isomorphism tool to find substructures

        :Parameter:

            - subg : graph object (from another molsys) to be searched

        :Returns:

            a list of lists with the (sorted) vertex indices of the substructure
        """
        maps = subgraph_isomorphism(subg.molg, self.molg, vertex_label=(subg.molg.vp.type, self.molg.vp.type))
        subs = []
        for m in maps:
            sl = list(m)
            sl.sort()
            if sl not in subs: subs.append(sl)
        return subs

    def find_fragment(self, frag):
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
                for ca in self._mol.conn[a]:
                    if self._mol.elems[ca] == "h":
                        f.append(ca)
            frags.append(f)
        return frags



