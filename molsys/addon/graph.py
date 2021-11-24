# -*- coding: utf-8 -*-

"""

       module to implement an addon feature: graphs using the graph_tool library

       NOTE: this is only imported by __init__.py if graph_tool is present

"""

import graph_tool
from graph_tool import Graph, GraphView
from graph_tool.topology import *
import copy
import numpy as np
import molsys
from molsys.util import elems

import uuid

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
        self.vert2atom = [] # this list maps vertex indices to the real atoms because we omit the hydrogens in the graph
        ig = 0
        self.molg.clear() # allways start from a fresh graph
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

    def make_comp_graph(self, elem_list = ["c"], idx = None):
        """
        Like make_graph this creates a graph for the mol object, but here we focus on a graph 
        for comparing molecular species in reactions. In the current graph we only look at the C-graph

        """
        if idx is None: idx = range(self._mol.natoms)
        # now add vertices
        self.vert2atom = []  
        ig = 0
        self.molg.clear()  # allways start from a fresh graph
        for i in idx:
            if self._mol.elems[i] in elem_list:
                self.molg.add_vertex()
                self.vert2atom.append(i)
                vtype = self._mol.atypes[i]
                # extract element and coordination number
                if "_" in vtype:
                    vtype = vtype.split("_")[0]
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
                        self.molg.add_edge(self.molg.vertex(i), self.molg.vertex(self.vert2atom.index(ja)))
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
        #GS TODO call plot_molgraph to do actual plotting

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
    def plot_mol_graph(fname, g, label=None, edge_label=None, size=1000, fsize=16, vsize=8, ptype = "pdf",method='arf'):
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
        if label:
            vlabel = label
        else:
            vlabel = "type"

        draw_g = g

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
    def find_subgraph(graph, subg, graph_property = None, subg_property = None, N=0):
        """
        use graph_tools subgraph_isomorphism tool to find substructures

        :Parameter:

            - graph : parent graph to be searched
            - subg  : graph to be found
            - N (int): number of subgraphs to find, N=0 is all (defaults to N=0)

        :Returns:

            a list of lists with the (sorted) vertex indices of the substructure
        """

        if graph_property is None: graph_property = graph.vp.type
        if subg_property is None: subg_property = subg.vp.type
        property_maps = (subg_property, graph_property)
        maps = subgraph_isomorphism(subg, graph, vertex_label=property_maps,max_n=N)
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

    # the following static methods are used in the fragments addon (to avoid graphtool dependence of the rest of fragments)
    @staticmethod
    def get_kcore(graph):
        return kcore_decomposition(graph).get_array()


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

    def check_sub(self, subg):
        """check if subg found in self.graph
        
        Args:
            subg (mol.graph objects): subgraph to be tested
        """
        subs = self.find_subgraph(self.molg, subg.molg, N=1)
        if subs != []:
            return True
        else:
            return False
        

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

    def util_graph(self, vertices, conn, atom_map=None, vtypes2 = None):
        """generates a fragment or atom graph

        Args:
            vertices (list of strings): vertex identifier
            conn (list of lists): connectivity 
            atom_map (list of list of ints, optional): atoms mapped by this vertex. Defaults to None.
            vtypes2 (list of strings, optional): alternative vertex identifiers. Defaults to None.

        Returns:
            graph: graph object    
        
        RS: this is currently just a helper function that produces the graph and returns it.
            we could consider to store this graph with a name in the graph addon

        """
        if vtypes2 is not None:
            assert len(vtypes2)==len(vertices)
        if atom_map is not None:
            assert len(atom_map) == len(vertices)
        g = Graph(directed=False)
        # now add vertices
        g.vp.type = g.new_vertex_property("string")
        if vtypes2 is not None:
            g.vp.types2 = g.new_vertex_property("string")
        if atom_map is not None:
            g.vp.atom_map = g.new_vertex_property("vector<int>")
        for i, v in enumerate(vertices):
            g.add_vertex()
            g.vp.type[i] = v
            if vtypes2 is not None:
                g.vp.types2[i] = vtypes2[i]
            if atom_map is not None:
                g.vp.atom_map[i] = atom_map[i]
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

    def unfilter_graph(self):
        self.molg.clear_filters()
        return

    def get_components(self, fidx= None):
        """Get all the components aka molecules from the atomic graph

        it adds a property map molid to the graph

        Returns:
            number of components found

        """
        if fidx is not None:
            # filter 
            self.filter_graph(fidx)
        label_components(self.molg, vprop=self.molg.vp.molid)
        ncomp = max(self.molg.vp.molid.a)+1
        if fidx is not None:
            self.unfilter_graph()
        comp = list(self.molg.vp.molid.a)
        return ncomp, comp


    @staticmethod
    def is_equal(molg1, molg2, use_fast_check=False):
        """helper function to identify if two molecular graphs are equal or not
        
        Args:
            molg1 (molecular graph): molecular graph to be compared to molg2 
            molg2 (molecular graph): molecular graph to be compared to molg1 
            use_fast_check (bool): will enforce a fast check based on the similarity rather than a graph isomorphisim

        Return:
            bool
        """

        error_code = 0

        if use_fast_check:

            similarity = graph_tool.topology.similarity(molg1,molg2)

            is_equal = (similarity > 0.99)

        else:

            try:

                e1 = molg1.get_edges()
                e2 = molg2.get_edges()
                
                if e1.shape[0] > 0 and e2.shape[0] > 0:

                    # quick exist?
                    vert1 = molg1.get_vertices()
                    vert2 = molg2.get_vertices()

                    if len(vert1) != len(vert2) or len(e1) != len(e2):
                       is_equal = False
                       return is_equal, error_code
                     
                    masterg = Graph(molg2)
                    masterg.add_vertex() 

                    vertex_maps = graph_tool.topology.subgraph_isomorphism(molg1, masterg, max_n=0, vertex_label=(molg1.vp.type,masterg.vp.type), edge_label=None, induced=False, subgraph=True, generator=False)

                    is_equal = len(vertex_maps) > 0

                    if is_equal:
                       for vi,vj in zip(molg1.vertices(),vertex_maps[0]):
                           if molg1.vp.type[vi] != molg2.vp.type[vj]:
                               is_equal = False
                               break
                       
                else:
                    # We don't have any edges... 

                    # compare vertices
                    vert1 = molg1.get_vertices()
                    vert2 = molg2.get_vertices()

                    if len(vert1) != len(vert2):
                       is_equal = False
                       return is_equal, error_code

                    is_equal = True

                    for v1, v2 in zip(vert1,vert2):
                        if molg1.vp.type[v1] != molg2.vp.type[v2]:
                            is_equal = False
                            break

            except:
                print("An error occured during the graph comparison")
                is_equal = False
                error_code = 1 


        return is_equal, error_code

