#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import molsys
from graph_tool import Graph
from graph_tool.topology import *
import numpy

class conngraph:

    def __init__(self,mol):
        self.mol = mol

    def make_graph(self, forbidden = None):
        self.molg = Graph(directed=False)
        ig = 0
        # setup vertices
        self.molg.vp.fix = self.molg.new_vertex_property("short")
        self.molg.vp.midx = self.molg.new_vertex_property("short")
        self.molg.vp.elem = self.molg.new_vertex_property("string")
        self.molg.vp.coord = self.molg.new_vertex_property("vector<double>")
        self.molg.vp.inter = self.molg.new_vertex_property("bool") # what is this good for? i have no idea -marco
        self.molg.vp.filled = self.molg.new_vertex_property("bool") # boolean for flood fill
        for i in xrange(self.mol.natoms):
            ig = self.molg.add_vertex()
            self.molg.vp.coord[ig] = self.mol.xyz[i,:]
            self.molg.vp.elem[ig] = self.mol.elems[i]
            self.molg.vp.inter[ig] = False
            self.molg.vp.midx[ig] = i
            if type(forbidden) != type(None) and int(ig) in forbidden:
                self.molg.vp.fix[ig] = 1
            else:
                self.molg.vp.fix[ig] = 0
        # setup edges
        self.molg.ep.act = self.molg.new_edge_property("bool")
        self.molg.ep.Nk = self.molg.new_edge_property("short")
        for i in xrange(self.mol.natoms):
            for j in self.mol.conn[i]:
                if j > i:
                    e = self.molg.add_edge(self.molg.vertex(i),self.molg.vertex(j))
                    self.molg.ep.act[e] = True
        # create Backup of the original graph for comparison
        self.keep = Graph(self.molg, directed=False)
        return
    
    def cut_to_2core(self):
        """
        Cuts graph to its 2-core
        
        Returns: Graph object
        """
        k = kcore_decomposition(self.molg).get_array()
        idxlist = numpy.argwhere(k==1).tolist()
        idx = []
        for i in idxlist:
            if self.molg.vp.fix[i[0]] == 0:
                idx.append(i[0])
        for v in reversed(sorted(idx)):
            self.molg.remove_vertex(v)
        #self.molg.remove_vertex(idx)
        return idxlist

    def graph2topo(self):
        m = molsys.topo()
        xyz = []
        natoms = 0
        conn = []
        elems = []
        for i,v in enumerate(self.molg.vertices()):
            xyz.append(self.molg.vp.coord[i])
            elems.append(self.molg.vp.elem[i])
            lconn = []
            for c in v.all_neighbours():
                #print v, c
                lconn.append(int(str(c)))
            conn.append(lconn)
            natoms +=1
        m.natoms = natoms
        m.set_xyz(numpy.array(xyz))
        m.set_conn(conn)
        m.set_elems(elems)
        m.set_atypes(natoms*['1'])
        m.set_cell(self.mol.get_cell())
        m.add_pconn()
        return m

    def remove_2conns(self):
        """
        From self.molg, removes all vertices with 2 connecting edges.
        
        Returns: boolean if the method could find any 2-connected vertices.
        """
        found_2conns = False
        for v in self.molg.vertices():
            if self.molg.vp.fix[v] == 0:
                neighbours = []
                for i in v.out_neighbours(): neighbours.append(i)
                if len(neighbours) == 2:
                    found_2conns = True
                    if not self.molg.edge(neighbours[0], neighbours[1]):
                        self.molg.add_edge(neighbours[0], neighbours[1])
                    self.molg.remove_vertex(v)
                    self.remove_2conns()
                    break
        return found_2conns
    
    def dissolve_ngon(self, max_n=100):
        """
        Dissolves a n-gon into a center vertex and n 2-connected vertices.
        
        Returns: boolean if the method could find any n-gons.
        """
        n = 3
        # find ngon
        iso = []
        while iso == []:
            ngon = self.create_ngon(n)
            ngon.vp.fix = ngon.new_vertex_property("short")

            for v in ngon.vertices():
                ngon.vp.fix[v] = 0
            iso = subgraph_isomorphism(ngon, self.molg, 1, (ngon.vp.fix, self.molg.vp.fix))
            n += 1
            if n > max_n:
                break
        # dissolve ngon
        if iso != []:
            new_xyz = numpy.array([0.0,0.0,0.0])
            coords = []
            subg = list(iso[0])
            new_v = self.molg.add_vertex()
            self.molg.vp.fix[new_v] = 0
            for i,v in enumerate(subg):
                self.molg.add_edge(v, new_v)
                self.molg.remove_edge(self.molg.edge(subg[i], subg[i-1]))
                new_xyz += self.molg.vp.coord[int(str(v))]
                coords.append(self.molg.vp.coord[int(str(v))])
            new_xyz /= n
            self.molg.vp.coord[new_v] = new_xyz
            return True
        else:
            return False
    
    def create_ngon(self, n):
        """
        Creates a n-gon as a Graph object.
        """
        ngon = Graph(directed=False)
        v = ngon.add_vertex()
        v_first = v
        for i in range(n-1):
            v_old = v
            v = ngon.add_vertex()
            ngon.add_edge(v_old, v)
        ngon.add_edge(v_first, v)
        return ngon

    def determine_Nk(self):
        """
        Assigns size of minimal ring to every edge
        """
        i = 0
        for e in self.molg.edges():
            i += 1
            self.molg.ep.act[e] = False
            self.molg.set_edge_filter(self.molg.ep.act)
            dist = shortest_distance(self.molg, source=e.source(), target= e.target())
            self.molg.ep.act[e] = True
            if dist < 2147483647:
                self.molg.ep.Nk[e] = dist+1
            else:
                self.molg.ep.Nk[e] = 0
            #print self.molg.ep.Nk[e]
        return

    def find_cluster_threshold(self):
        """
        Finds thresholds.
        Needs Nk values of the edges -> determine_Nk() has to be called before calling this method.
        """
        self.threshes = []
        Nk = [i for i in self.molg.ep.Nk.get_array() if i != 0]
        Nk.sort()
        for i in range(len(Nk)-1):
            if Nk[i+1]-Nk[i]>2:
                self.threshes.append(Nk[i+1])
        return self.threshes

    def get_clusters(self):
        """
        Get the clusters of the MOF.
        """
        try:
            assert self.threshes
        except:
            self.find_cluster_threshold()
        # remove side chains
        def forloop():
            broken = False
            for e in self.molg.edges():
                # since the indices of the edges will change when one edge is removed, the for loop has to be
                # restarted every time an edge is deleted, or else there will be problems when the next edge 
                # should be deleted...
                if self.molg.ep.Nk[e] == 0:
                    # found a sidechain... now use flood_fill to find out which side of the
                    # edge belongs to the sidechain and which to the main framework
                    src = e.source()
                    trg = e.target()
                    self.molg.remove_edge(e)
                    self.molg.vp.filled.set_value(False)
                    struc1 = self.flood_fill(self.molg, src, [])
                    struc2 = self.flood_fill(self.molg, trg, [])
                    # the larger structure is the framework
                    if len(struc1) < len(struc2):
                        for i in reversed(sorted(struc1)):
                            self.molg.remove_vertex(i)
                    elif len(struc2) < len(struc1):
                        for i in reversed(sorted(struc2)):
                            self.molg.remove_vertex(i)
                    else:
                        # both structures have same size
                        pass
                    broken = True
                    break
            return broken
        while forloop():
            pass
        ### set edge filter
        self.molg.vp.filled.set_value(False)
        self.molg.ep.act.set_value(True)
        if self.threshes != []:
            thresh = self.threshes[0]
        else:
            thresh = 0
        for e in self.molg.edges():
            if self.molg.ep.Nk[e] >= thresh:
                self.molg.ep.act[e] = False
        self.molg.set_edge_filter(self.molg.ep.act)
        ### perform flood fill
        clusters = []
        while False in list(self.molg.vp.filled.get_array()):
            vidx = list(self.molg.vp.filled.get_array()).index(0)
            vstart = self.molg.vertex(vidx)
            cluster = self.flood_fill(self.molg, self.molg.vertex(vstart), [])
            cluster = map(int, cluster)
            clusters.append(cluster)
        self.molg.clear_filters()
        self.clusters = clusters
        return clusters
    
    def get_cluster_atoms(self):
        """
        Get all atoms in the clusters of the conngraph.
        Needs clusters in self.clusters -> call get_clusters before calling this method!
        
        Returns: 
        -List of list of vertices of each cluster in the backup "self.keep" graph.
        -List of list of atom ids of each cluster
        """
        try:
            assert self.clusters
        except:
            self.get_clusters()
        midx_list = self.keep.vp.midx.get_array().tolist()
        # set edge filter
        self.keep.vp.filled.set_value(False)
        self.keep.ep.act.set_value(True)
        thresh = self.threshes[0]  # temporary...
        for e in self.molg.edges():
            if self.molg.ep.Nk[e] >= thresh:
                src = e.source()
                trg = e.target()
                midx = (self.molg.vp.midx[src], self.molg.vp.midx[trg])
                atomid = (midx_list.index(midx[0]), midx_list.index(midx[1]))
                keepedge = self.keep.edge(atomid[0], atomid[1])
                self.keep.ep.act[keepedge] = False
        self.keep.set_edge_filter(self.keep.ep.act)
        # then find all connecting atoms
        clusters_vertices = []
        clusters_atoms = []
        for c in self.clusters:
            for vid in c:
                v = self.molg.vertex(vid)
                midx = self.molg.vp.midx[v]
                atomid = midx_list.index(midx)
                atom = self.keep.vertex(atomid)
                if self.keep.vp.filled[atom] == False:
                    this_cluster_verts = self.flood_fill(self.keep, atom, [])
                    clusters_vertices.append(this_cluster_verts)
                    this_cluster_atoms = []
                    for i in this_cluster_verts:
                        this_cluster_atoms.append(self.keep.vp.midx[i])
                    clusters_atoms.append(this_cluster_atoms)
        self.molg.clear_filters()
        return clusters_vertices, clusters_atoms
    
    def make_topo_graph(self, verbose=True):
        try:
            assert self.clusters
        except:
            self.get_clusters()
        tm = molsys.topo()
        tm.natoms = len(self.clusters)
        tm.set_empty_conn()
        xyz = []
        elems = []
        for i, c in enumerate(self.clusters):
            ext_bond = []
            cluster_atoms = self.clusters[i]
            cidx = []
            for ia in cluster_atoms:
                cidx.append(self.molg.vp.midx[ia])
                via = self.molg.vertex(ia)
                for j in via.all_neighbours():
                    if j not in cluster_atoms:
                        # thus bond is an external bond
                        ext_bond.append(int(str(j)))
            xyz.append(self.mol.get_com(cidx))
            #xyz.append(self.center(cxyz))
            if verbose: print "cluster %s consisting of %d atoms is %d times connected" % (str(i), 
                    len(cluster_atoms), len(ext_bond))
            # now check to which clusters these external bonds belong to
            for ea in ext_bond:
                for ji, j in enumerate(self.clusters):
                    if ea in j:
                        if verbose: print " -> bonded to cluster ", ji
                        tm.conn[i].append(ji)
                        break
        ### check for consistence of conn
        for i in xrange(tm.natoms):
            if len(tm.conn[i]) == 4:
                elems.append('c')
            elif len(tm.conn[i]) == 2:
                elems.append('o')
            else:
                elems.append('n')
            for j in tm.conn[i]:
                if j>i:
                    if not i in tm.conn[j]:
                        if verbose: print "Fragment topology is inconsitent"
        tm.set_xyz(numpy.array(xyz))
        tm.set_elems(elems)
        tm.set_atypes(tm.natoms*['0'])
        tm.set_cell(self.mol.get_cell())
        tm.add_pconn()
        topograph = self.__class__(tm)
        topograph.make_graph()
        if verbose: print self.threshes
        return topograph

    def get_all_cs(self, depth, tg=None):
        """
        Calculates (coordination sequence) all cs values of the conngraph.
        depth = maximum level
        """
        if tg == None:
            tg = self.make_topo_graph(False)
        cs_list = []
        for i in range(len(self.clusters)):
            cs = self.get_cs(depth, i, tg=tg)
            if cs not in cs_list:
                cs_list.append(cs)
        return cs_list

    def get_cs(self, depth, start_vertex=0, start_cell=numpy.array([0,0,0]), tg=None):
        """
        Calculates the cs (coordination sequence) values of the vertex specified in start_vertex and start_cell.
        depth = maximum level
        tg = topograph of this class - pre-calculate it if you want a faster execution speed when
             calculating multiple cs values on the same conngraph.
        """
        def contains(l, obj):
            # the normal construction "if x not in y" does not work if arrays are somehow involved inside a list
            # thus, we need this helper function which is probably terribly slow (if anyone has a better solution, please let me hear it)
            found = False
            for i in l:
                if i[0] == obj[0] and (i[1] == obj[1]).all():
                    found = True
            return found
        # create cs with an appropriate length.
        cs = depth*[0]
        if tg == None:
            # topograph should only be generated once, if it's not there yet, so the program runs more quickly.
            tg = self.make_topo_graph(False)
            # the data of the topo object is found at tg.mol
        # now, if we want to have a cs_n value, we need to get all neighbours, then do the same thing for all neighbour's
        # neighbours, and continue as many times as necessary.
        # however, we must not start looking at the neighbour's neighbours before we looked at ALL neighbours, and set them
        # into the ignore list!
        ignore_list = [[start_vertex, start_cell]]
        neighbours = [[start_vertex, start_cell]]
        for level in range(depth):
            neighbours2 = []
            #print "--------- level "+str(level)+" -----------"
            #print "neighbours: "+str(neighbours)
            for n in neighbours:
                # get the neighbours of the neighbours and add them to the list neighbours2
                visited = self.get_cs1(n[0], n[1], tg)
                for v in visited:
                    if not contains(neighbours2, v):
                        if not contains(ignore_list, v):
                            neighbours2.append(v)
            #print "ignore_list: " +str(ignore_list)
            # put the neighbours2 into the ignore list
            for n2 in neighbours2:
                ignore_list.append(n2)
            #print "neighbours2: "+str(neighbours2)
            # the neighbours2 are all the vertices which can be reached with the cs_level operation.
            cs[level] = len(neighbours2)
            # if we want to repeat the procedure for the cs_(level+1) operation we have to make the neighbours2 to the neighbours
            neighbours = neighbours2
            # and remove neighbours2 because of side effects
            del(neighbours2)
        return cs
    
    def get_cs1(self, start_vertex=0, start_cell=numpy.array([0,0,0]), tg=None):
        """
        This function will return all vertices, which are connected to the vertices start_vertex in the cell start_cell.
        """
        assert self.clusters
        visited = []
        if tg==None:
            tg = self.make_topo_graph(False)
        # loop over all neighbouring clusters and add them to the list
        for nj, j in enumerate(tg.mol.conn[start_vertex]):
            current_vertex = j
            current_cell = start_cell + tg.mol.pconn[start_vertex][nj]
            visited.append([current_vertex, current_cell])
        return visited
    
    def get_vertex_symbol(self, start_vertex):
        """
        In the vertex symbol, the number of the shortest ring at each angle of a vertex is given with
        the number of such rings in brackets (originally: subscript, but this can't be realized in
        python).
        Relevant literature: O. Delgado-Friedrichs, M. O'Keeffe, Journal of Solid State Chemistry 178, 2005, p. 2480ff.
        """
        self.molg.vp.filled.set_value(False)
        self.molg.vp.filled[start_vertex] = True
        # ACHTUNG BAUSTELLE !!!!
        return
    
    def center(self,xyz):
        cell_abc = self.mol.cellparams[:3]
        fix = xyz[0,:]
        a = xyz[1:,:] - fix
        xyz[1:,:] -= cell_abc*numpy.around(a/cell_abc)
        center = numpy.sum(xyz, axis = 0)/numpy.shape(xyz)[0]
        return center
    
    def handle_islands(self, thresh_small=0.2, silent=False):
        """
        Handles parts of the structure which are not connected to the rest ("islands"):
        - if the island is smaller than the size of largest island times thresh_small, it will be deleted.
        - else, a warning will be printed out, since the structure might be bugged or an interpenetrated network.
        Returns:
        - n_removed_atoms: Number of removed atoms
        - multiple_large_islands: boolean if there are multiple large islands
        """
        multiple_large_islands = False
        islands = self.get_islands()
        if len(islands) < 2:
            # we only have 1 island, so nothing has to be done
            return 0, False
        else:
            # ok we have multiple islands, time to find out what we have here.
            # first we sort the islands by size
            island_sizes = map(len, islands)
            sort_indices = numpy.array(island_sizes).argsort()
            biggest_size = island_sizes[sort_indices[-1]]
            remove_list = []
            # then we compare the sizes of the smaller islands with the size of the largest island
            for index, i in enumerate(island_sizes):
                # don't compare the biggest island with itself though...
                if index != sort_indices[-1]:
                    ratio = float(i)/float(biggest_size)
                    if ratio < thresh_small:
                        for v in reversed(sorted(islands[index])):
                            remove_list.append(v)
                    else:
                        multiple_large_islands = True
                        if not silent: 
                            print("Warning: Two large parts of the structure are not connected to each other. Check the structure.")
            # DELETE THE SMALL ISLANDS LIKE GLOBAL WARMING !!!!
            self.molg.remove_vertex(remove_list)
            n_removed_atoms = len(remove_list)
            return n_removed_atoms, multiple_large_islands
    
    def get_islands(self):
        """
        Finds all islands.
        """
        self.molg.vp.filled.set_value(False)
        remain_list = range(self.molg.num_vertices())
        islands = []
        while len(remain_list) > 0:
            l = self.flood_fill(self.molg, self.molg.vertex(remain_list[0]), [])
            #print len(l) # debug message - prints sizes of all islands...
            islands.append(l)
            for i in l:
                remain_list.remove(int(i))
        return islands
        
    def flood_fill(self, graph, vertex, return_list=[]):
        """
        Uses flood fill to find all vertices that are connected to the starting vertex.
        Caution: You might want to reset the vertex property "filled" before calling flood_fill!
        
        Parameters:
        - vertex: starting vertex
        - return_list: list of vertices which have already been iterated (when calling the function, you might have to force this to be [])
        - graph: The graph in which flood fill should be performed
        
        Returns:
        - list of all vertices that could be reached by flood fill.
        """
        # perform flood fill
        graph.vp.filled[vertex] = True
        return_list.append(vertex)
        for n in vertex.all_neighbours():
            if graph.vp.filled[n] == False:
                self.flood_fill(graph, n, return_list)
        return return_list
    
    def search_and_destroy(self, printout=False):
        """
        Uses the cut_to_2core, remove_2conns and dissolve_ngon methods to minimize the graph of a building block.
        
        Parameters:
        - printout (bool): if True, this will create a series of PNG images, which show (and slow) the process.
        """
        pc = 0
        def printout_func(printout_counter):
            graph_draw(self.molg, vertex_text=self.molg.vp.fix, vertex_font_size=10, output_size=(500, 500), \
                output="test%d.png"%printout_counter)
            printout_counter += 1
            return printout_counter
        if printout: pc = printout_func(pc)
        self.cut_to_2core()
        if printout: pc = printout_func(pc)
        self.remove_2conns() # initial removal of 2conns
        if printout: pc = printout_func(pc)
        while True:
            # dissolve ngons and remove 2conns until one of the functions returns False.
            if not self.dissolve_ngon(): break
            if printout: pc = printout_func(pc)
            if not self.remove_2conns(): break
            if printout: pc = printout_func(pc)
        return


