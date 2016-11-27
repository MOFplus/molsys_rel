# -*- coding: utf-8 -*-

"""

       module to implement an addon feature: graphs using the graph_tool library

       NOTE: this is only imported by __init__.py if graph_tool is present

"""

from graph_tool import Graph
from graph_tool.draw import graph_draw

class graph:

    def __init__(self, mol):
        """
        instantiate a graph object which will be attached to the parent mol

        :Parameter:

             - mol : a mol type object (can be a derived type like bb or topo as well)
        """
        self._mol = mol
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
            if self._mol.elems[i] != "h":
                self.molg.add_vertex()
                self.vert2atom.append(i)
                vtype = self._mol.atypes[i]
                if "_" in vtype:
                    vtype = vtype.split("_")[0]
                self.molg.vp.type[ig] = vtype
                ig += 1
        self.nvertices = len(self.vert2atom)
        # now add edges ... only bonds between vertices
        for i in xrange(self.nvertices):
            ia = self.vert2atom[i]
            for ja in self._mol.conn[ia]:
                if ja>ia:
                    if ja in self.vert2atom:
                        print "bond from %d to %d" % (ia, ja)
                        print self._mol.atypes[ia], self._mol.atypes[ja]
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
        graph_draw(self.molg, vertex_text=self.molg.vp.type, vertex_font_size=10,  \
            output_size=(size, size), output=fname+".pdf")
        return


