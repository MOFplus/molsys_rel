#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import molsys
from graph_tool import Graph
from graph_tool.topology import *
import numpy
import copy
from weaver import user_api
from string import ascii_lowercase

class conngraph:
    """ This is the "conngraph" class
    It is the basis of the molgraph and topograph classes, which use graph theory for the deconstruction of
    MOF structures, or the analysis of topologies, respectively. It provides the important ground functions
    used by both classes, like make_graph, flood_fill and center.
    """

    def __init__(self, mol):
        """
        :Parameters:
        - mol: molsys.mol or molsys.topo object
        """
        self.mol = mol
        self.make_graph()

    def make_graph(self, forbidden = None):
        """ Create a Graph object from a molsys.mol object. """
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
        """
        k = kcore_decomposition(self.molg).get_array()
        idxlist = numpy.argwhere(k==1).tolist()
        idx = []
        for i in idxlist:
            if self.molg.vp.fix[i[0]] == 0:
                idx.append(i[0])
        for v in reversed(sorted(idx)):
            self.molg.remove_vertex(v)
        return idxlist

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
    
    def handle_islands(self, thresh_small=0.2, silent=False):
        """
        Handles parts of the structure which are not connected to the rest ("islands").
        
        :Parameters:
        - thresh_small: If the island is smaller than the size of the largest island multiplied by
          thresh_small, it will be deleted.
        - silent: if True, a warning will be printed out when there are multiple large islands.
        
        :Returns:
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
        """ Finds all islands. """
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
        
        :Parameters:
        - vertex: starting vertex
        - return_list: list of vertices which have already been iterated (when calling the function, you might have to force this to be [])
        - graph: The graph in which flood fill should be performed
        
        :Returns:
        - list of all vertices that could be reached by flood fill.
        """
        # perform flood fill
        graph.vp.filled[vertex] = True
        return_list.append(vertex)
        for n in vertex.all_neighbours():
            if graph.vp.filled[n] == False:
                self.flood_fill(graph, n, return_list)
        return return_list
    
class molgraph(conngraph):
    """ This class handles the deconstruction of a MOF structure into a graph. """
    
    def determine_Nk(self):
        """ Assigns size of minimal ring to every edge """
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

    def find_cluster_treshold(self):
        self.find_cluster_threshold()
        return self.threshes

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
        
        The clusters (saved in self.clusters) are lists of atoms, which belong in one group or form 
        one 'building block'. They are firstly defined by ring sizes: If the Nk value is higher or equal 
        than the threshold, this bond will be defined as an inter-cluster bond.
        The clusters can be later modified to create different representations of the framework, for 
        example by summarizing all organic clusters which are bonded to each other.
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
        
        :Returns: 
        -List of list of vertices of each cluster in the backup "self.keep" graph.
        -List of list of atom ids of each cluster
        """
        try:
            assert self.clusters
        except:
            self.get_clusters()
        midx_list = self.keep.vp.midx.get_array().tolist()
        # set edge filter
        for i, c in enumerate(self.clusters):
            cluster_atoms = self.clusters[i]
            for ia in cluster_atoms:
                via = self.molg.vertex(ia)
                for j in via.all_neighbours():
                    if j not in cluster_atoms:
                        # found external bond, set edge filter here
                        midx = (self.molg.vp.midx[via], self.molg.vp.midx[j])
                        atomid = (midx_list.index(midx[0]), midx_list.index(midx[1]))
                        keepedge = self.keep.edge(atomid[0], atomid[1])
                        self.keep.ep.act[keepedge] = False
        self.keep.set_edge_filter(self.keep.ep.act)
        # then find all connecting atoms
        clusters_vertices = []
        clusters_atoms = []
        self.keep.vp.filled.set_value(False)
        for ic, c in enumerate(self.clusters):
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
                        # set atomtypes in self.mol
                        self.mol.atypes[self.keep.vp.midx[i]] = ic
                    clusters_atoms.append(this_cluster_atoms)
        self.keep.clear_filters()
        return clusters_vertices, clusters_atoms
    
    def cluster_conn(self):
        """
        Helper function that returns a list, which describes, which clusters are connected with each other.
        """
        assert self.clusters
        cluster_conn = []
        for i, cluster_atoms in enumerate(self.clusters):
            this_cluster_conn = []
            ext_bond = []
            for ia in cluster_atoms:
                via = self.molg.vertex(ia)
                for j in via.all_neighbours():
                    if j not in cluster_atoms:
                        ext_bond.append(int(str(j)))
            #print "cluster %s consisting of %d atoms is %d times connected" % (str(i), len(cluster_atoms), len(ext_bond))
            # now check to which clusters these external bonds belong to
            for ea in ext_bond:
                for ji, j in enumerate(self.clusters):
                    if ea in j:
                        this_cluster_conn.append(ji)
            #            print " -> bonded to cluster ", ji
                        break
            cluster_conn.append(this_cluster_conn)
        return cluster_conn
    
    def get_bbs(self):
        """
        Returns the building blocks (BBs) of the MOF.
        Since the definition of "building block" is arbitrary, this function will do the following things:
        - if multiple 2-connected clusters are connected to each other, these clusters will be summarized
          as one BB.
        The function will OVERWRITE self.clusters and replace it with the newly generated building blocks!
        """
        try:
            assert self.clusters
        except:
            self.get_clusters()
        # summarize 2-connected clusters:
        stop = False
        while not stop:
            stop = True
            # find external bonds
            cluster_conn = self.cluster_conn()
            # find out if 2-connected clusters are bonded to other 2-connected clusters
            for i, c in enumerate(cluster_conn):
                if len(c) == 2:
                    for b in c:
                        if len(cluster_conn[b]) == 2:
                            # and if they are, create new clusters which contain everything the first clusters contained
                            self.clusters.append(self.clusters[i] + self.clusters[b])
                            # and then remove those clusters
                            if i > b:
                                del self.clusters[i]
                            del self.clusters[b]
                            if i < b:
                                del self.clusters[i]
                            stop = False
                            break
                if not stop:
                    break
        return self.clusters
    
    def detect_organicity(self, organic_elements = ["h", "b", "c", "n", "o", "f", "si", "p", "s", "cl", "as", "se", "br", "i"]):
        """
        Finds out whether a cluster classifies as "organic" or "inorganic".
        
        :Parameters:
        - organic_elements: This is a list which defines, which elements are allowed inside
          an "organic" building block. All element symbols must be lower case. Default:
          H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I
          
        :Returns:
        - self.cluster_organicity: list, each element is True (for organic) or False (for 
          inorganic), and the index is equivalent to the index of the cluster in self.clusters
        """
        try:
            assert self.clusters
        except:
            self.get_clusters()        
        # go through each atom in every cluster and compare the elements
        self.cluster_organicity = []
        for cluster in self.clusters:
            org = True
            for i in cluster:
                v = self.molg.vertex(i)
                midx = self.molg.vp.midx[v]
                if self.mol.elems[midx].lower() not in organic_elements:
                    org = False
                    break
            self.cluster_organicity.append(org)
        return self.cluster_organicity
    
    def get_bbs_by_organicity(self):
        """
        Groups all organic clusters, which are connected with each other, into a single cluster.
        Like self.get_bbs, this method will OVERWRITE self.clusters.
        """
        try:
            assert self.cluster_organicity
        except:
            self.detect_organicity()
        stop = False
        while not stop:
            stop = True
            # find external bonds
            cluster_conn = self.cluster_conn()
            # find out if organic clusters are bonded to other organic clusters
            for i, c in enumerate(cluster_conn):
                if self.cluster_organicity[i] == True:
                    for b in c:
                        if self.cluster_organicity[b] == True:
                            # and if they are, create new clusters which contain everything the first clusters contained
                            self.clusters.append(self.clusters[i] + self.clusters[b])
                            # and then remove those clusters
                            if i > b:
                                del self.clusters[i]
                            del self.clusters[b]
                            if i < b:
                                del self.clusters[i]
                            stop = False
                            # recompute organicity, because of the changed indices
                            self.detect_organicity()
                            break
                if not stop:
                    break
        return self.clusters
    
    def make_topograph(self, verbose=True, allow_2conns=False):
        """
        Create the topograph of the topology of the MOF.
        
        :Parameters:
        - verbose: if True, there will be some information printed out
        - allow_2conns: if True, 2-connected vertices will be allowed in the topograph
        """
        try:
            assert self.clusters
        except:
            self.get_clusters()
        tm = molsys.topo()
        tm.natoms = len(self.clusters)
        tm.set_empty_conn()
        xyz = []
        elems = []
        for i, cluster_atoms in enumerate(self.clusters):
            ext_bond = []
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
        tg = topograph(tm, allow_2conns)
        tg.make_graph()
        if verbose: print self.threshes
        return tg
    

class topograph(conngraph):
    """ This class handles the analysis of nets using graph theory. """
    
    def __init__(self, mol, allow_2conns=False):
        """
        :Parameters:
        - mol: molsys.topo object
        - allow_2conns: if False, all 2-connected vertices will be deleted.
        """
        self.mol = mol
        if not allow_2conns:
            self.midx_2conns = self.remove_2conns_from_mol()
        self.make_graph()
        return
    
    def remove_2conns_from_mol(self):
        """ Removes all vertices with 2 connecting edges from self.mol """
        # delete atoms
        delete_list = []
        for i in range(self.mol.natoms):
            if len(self.mol.conn[i]) == 2:
                delete_list.append(i)
        for i in reversed(sorted(delete_list)):
            # retain connectivity information
            connected = []
            for j in self.mol.conn[i]:
                connected.append(j)
            self.mol.conn[connected[0]].append(connected[1])
            self.mol.conn[connected[1]].append(connected[0])
            # now delete the atom
            self.mol.delete_atom(i)
        # recompute pconn
        self.mol.add_pconn()
        return delete_list
    
    def graph2topo(self):
        """
        Returns a molsys.topo object (topology) created from the Graph,
        which can be used to create a graphical representation of the Graph viewable
        in molden, VMD or similar programs.
        """
        return copy.deepcopy(self.mol)
    
    def get_all_cs(self, depth=10, use_atypes=False):
        """
        Calculates all cs (coordination sequence) values of the graph.
        This function just loops over all vertices and calls get_cs for each one
        
        :Parameters:
        - depth: maximum level
        - use_atypes: if this is True, then every vertex with the same atomtype will only be calculated once. 
                      (do NOT use this for topographs deconstructed from a molgraph !!!)
        """
        if use_atypes:
            vertexlist = []
            found_atypes = []
            for i in range(self.mol.natoms):
                if self.mol.atypes[i] not in found_atypes:
                    found_atypes.append(self.mol.atypes[i])
                    vertexlist.append(i)
        else:
            vertexlist = range(self.mol.natoms)
        cs_list = []
        for i in vertexlist:
            cs = self.get_cs(depth, i)
            cs_list.append(cs)
        return cs_list

    def get_cs(self, depth, start_vertex=0, start_cell=numpy.array([0,0,0])):
        """
        Calculates the cs (coordination sequence) values of the vertex specified in start_vertex and start_cell.
        
        :Parameters:
        - depth: maximum level
        - start_vertex: ID of vertex, of which the cs value should be calculated
        - start_cell: starting cell of the function
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
                visited = self.get_cs1(n[0], n[1])
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
    
    def get_cs1(self, start_vertex=0, start_cell=numpy.array([0,0,0])):
        """ This function will return all vertices, which are connected to the vertex start_vertex in the cell start_cell. """
        visited = []
        # loop over all neighbouring clusters and add them to the list
        for nj, j in enumerate(self.mol.conn[start_vertex]):
            current_vertex = j
            current_cell = start_cell + self.mol.pconn[start_vertex][nj]
            visited.append([current_vertex, current_cell])
        return visited

    def get_all_vs(self, use_atypes=False, wells = False):
        """
        Calculates all vertex symbols of the graph.
        
        :Parameters:
        - use_atypes: if this is True, then every vertex with the same atomtype will only be calculated once.
        - wells: If True, the Wells and the long symbols will be returned. If False, only the long symbols are used
        """
        if use_atypes:
            vertexlist = []
            found_atypes = []
            for i in range(self.mol.natoms):
                if self.mol.atypes[i] not in found_atypes:
                    found_atypes.append(self.mol.atypes[i])
                    vertexlist.append(i)
        else:
            vertexlist = range(self.mol.natoms)
        vs_list = []
        supercells = [None, copy.deepcopy(self.mol)]
        keep = copy.deepcopy(self.mol)
        for i in vertexlist:
            success = False
            supercell_size = 2
            while not success:
                if supercell_size > len(supercells)-1:
                    self.mol = copy.deepcopy(keep)
                    self.mol.make_supercell([supercell_size]*3)
                    supercells.append(copy.deepcopy(self.mol))
                else:
                    self.mol = supercells[supercell_size]
                success = True
                self.make_graph()
                try:
                    ws, ls = self.get_vertex_symbol(i)
                except ValueError:
                    supercell_size += 1
                    success = False
            self.mol = keep
            self.make_graph()
            if wells: 
                vs = (ws, ls)
            else:
                vs = ls
            vs_list.append(vs)
        return vs_list
    
    def get_vertex_symbol(self, start_vertex):
        """
        In the vertex symbol, the number of the shortest ring at each angle of a vertex is given with
        the number of such rings in brackets (originally: subscript, but this can't be realized in
        python).
        Relevant literature: O. Delgado-Friedrichs, M. O'Keeffe, Journal of Solid State Chemistry 178, 2005, p. 2480ff.
        """
        if self.mol.natoms == 1:
            raise ValueError("Topology only consists of only one vertex")
        self.molg.vp.filled.set_value(False)
        self.molg.vp.filled[start_vertex] = True
        vertex_symbol = []
        paths = []
        for source in self.molg.vertex(start_vertex).all_neighbours():
            for target in self.molg.vertex(start_vertex).all_neighbours():
                if source < target:
                    self.molg.set_vertex_filter(self.molg.vp.filled, inverted=True)
                    asp = all_shortest_paths(self.molg, source, target)
                    self.molg.clear_filters()
                    append_list = []
                    for p1 in asp:
                        path = p1.tolist()
                        path = map(int, path)
                        p2 = [start_vertex]+path
                        vol = self.get_cycle_voltage(p2)
                        if vol.any() != numpy.zeros(3).any():
                            raise ValueError("Cycle with non zero voltage detected")
                        path.append(start_vertex)
                        append_list.append(path)
                    if len(append_list) != 0:
                        vertex_symbol.append((len(append_list[0]), len(append_list)))
        ws = self.compute_wells_symbol(vertex_symbol)
        ls =  self.compute_long_symbol(vertex_symbol)
        return ws, ls

    def compute_wells_symbol(self, clist):
        symbol = ""
        clist = numpy.array(clist)[:,0].tolist()
        sclist = sorted(set(clist))
        for i, s in enumerate(sclist):
            count = clist.count(s)
            if count != 1:
                symbol += "%s^%s." % (s,count)
            else:
                symbol += "%s." % s
        return symbol[:-1]

    def compute_long_symbol(self,clist):
        symbol = ""
        dtype = [("length",int),("number",int)]
        clist = numpy.array(clist,dtype=dtype)
        sclist = numpy.sort(clist, order=["length","number"]).tolist()
        for i, s in enumerate(sclist):
            if s[1]==1:
                symbol += "%s." % s[0]
            else:
                symbol += "%s(%s)." % (s[0],s[1]) 
        return symbol[:-1]

    def get_cycle_voltage(self, cycle):
        cycle.append(cycle[0])
        vol = numpy.zeros(3)
        for i in range(len(cycle)-1):
            cidx = self.mol.conn[cycle[i]].index(cycle[i+1])
            vol += self.mol.pconn[cycle[i]][cidx]
        return vol

    def get_unique_vd(self, cs, vs, atype = True):
        """
        Returns the unique cs values and vertex symbols ("vertex descriptors", vd). 
        :Parameters:
        - cs: list of cs values
        - vs: list of vertex symbols
        - atype: If this is True, the function also changes the atomtypes in the 
          topograph (i.e. the "vertex types"), according to the different vertex 
          descriptors.
        """
        assert type(cs) == list
        assert type(vs) == list
        assert len(vs) == len(cs)
        atypes = []
        uvd = []
        atcount = 0
        for c,v in zip(cs,vs):
            vd = tuple([tuple(c),v])
            if vd not in uvd: 
                uvd.append(vd)
                atypes.append(str(atcount))
                atcount += 1
            else:
                atypes.append(str(uvd.index(vd)))
        ucs = []
        uvs = []
        for i in uvd:
            ucs.append(i[0])
            uvs.append(i[1])
        if atype: self.mol.set_atypes(atypes)
        return ucs, uvs

    def build_coordination_pattern(self,pattern, cn):
        assert type(pattern) == list
        assert len(pattern) == 2
        assert type(cn) == list
        assert len(cn) == pattern[0]-1 + pattern[1]-1
        assert cn.count(pattern[0]) == cn.count(pattern[1]) == 0
        ### build subgraph
        patg = Graph(directed=False)
        patg.vp.cn = self.molg.new_vertex_property("short")
        for i in pattern:
            v = patg.add_vertex()
            patg.vp.cn[v] = i
        patg.add_edge(patg.vertex(0),patg.vertex(1))
        for i, c in enumerate(pattern):
            for j in range(c-1):
                v = patg.add_vertex()
                patg.vp.cn[v] = cn[i+j]
                patg.add_edge(v, patg.vertex(i))
        return patg

    def search_coordination_pattern(self,patg):
        assert type(patg) == Graph
        self.molg.vp.cn = self.molg.new_vertex_property("short")
        for v in self.molg.vertices():
            self.molg.vp.cn[v] = len(list(v.out_neighbours()))
        maps = subgraph_isomorphism(patg, self.molg, vertex_label =
                (patg.vp.cn, self.molg.vp.cn))
        subs = []
        for m in maps:
            sl = list(m)
            sl.sort()
            if sl not in subs: subs.append(sl)
        return subs

    def collapse_subs(self, subs, pattern):
        dl = []
        for s in subs:
            center = []
            v = self.molg.add_vertex()
            for vidx in s: 
                vi = self.molg.vertex(vidx)
                if self.molg.vp.cn[vi] in pattern:
                    center.append(vidx)
                    if vidx not in dl: dl.append(vidx)
            self.mol.set_unit_mass()
            xyz = self.mol.get_com(center)
            self.molg.vp.coord[v] = xyz
            self.mol.insert_atom('c','n',xyz,center[0],center[1])
            ### coordinates
            for vidx in s:
                vi = self.molg.vertex(vidx)
                if self.molg.vp.cn[vi] not in pattern:
                    self.molg.add_edge(vi, v)
                    self.mol.conn[-1].append(vidx)
                    self.mol.conn[vidx].append(self.mol.natoms-1)
        for v in reversed(sorted(dl)):
            self.molg.remove_vertex(v)
            self.mol.delete_atom(v)
        self.mol.add_pconn()


class topotyper(object):
    """ Wrapper class which combines molgraph and topograph for the deconstruction of MOF structures. """
  
    def __init__(self, mol, split_by_org=True):
        """
        :Parameters:
        - mol: molsys.mol object which should be deconstructed
        - split_by_org: if True, organicity is used for defining building blocks.
        """
        self.mg = molgraph(mol)
        self.api = None
        self.deconstruct(split_by_org)
        return
 
    def deconstruct(self, split_by_org=True):
        """ perform deconstruction """
        self.mg.handle_islands()
        self.mg.determine_Nk()
        self.mg.find_cluster_threshold()
        self.mg.get_clusters()
        if split_by_org:
            self.mg.get_bbs_by_organicity()
        else:
            self.mg.get_bbs()
        self.tg = self.mg.make_topograph(False)
        cs = self.tg.get_all_cs()
        vs = self.tg.get_all_vs()
        self.cs, self.vs = self.tg.get_unique_vd(cs, vs)
        return

    def get_net(self):
        """ Connects to mofplus API to compare the cs and vs value to the database, and return the topology. """
        self.api = user_api(experimental=False)
        self.nets = self.api.search_cs(self.cs, self.vs)
        return self.nets

    def tg_atypes_by_isomorphism(self):
        """
        This method will modify the atomtypes in the topograph (i.e. "vertex types"): If
        two vertices with the same atomtype have different structures, the atomtype will be
        changed.
        Possible todo: compare the elements with each other
        """
        cv, ca = self.mg.get_cluster_atoms()
        bb_molgraphs = []
        try:
            assert self.mg.cluster_organicity
        except:
            self.mg.detect_organicity()
        # Remove all 2-connected clusters, since they aren't in the topograph
        list2c = self.tg.midx_2conns
        for i in reversed(sorted(list2c)):
            del ca[i]
        # check for isomorphism
        for i, atoms in enumerate(ca):
            m = self.mg.mol.new_mol_by_index(atoms)
            mg = molgraph(m)
            #print "i "+str(i)
            j = len(bb_molgraphs) - 1
            while j >= 0:
                mg2 = bb_molgraphs[j]
                #print "  j "+str(j)
                if self.tg.mol.atypes[j] == self.tg.mol.atypes[i]:
                    if not isomorphism(mg.molg, mg2.molg):
                        self.tg.mol.atypes[i] += "#"
                        j = len(bb_molgraphs)
                j -= 1
            bb_molgraphs.append(mg)
        for i in range(len(self.tg.mol.atypes)):
            n = self.tg.mol.atypes[i].count("#")
            if n > 0:
                self.tg.mol.atypes[i] = self.tg.mol.atypes[i].replace("#", "")
                self.tg.mol.atypes[i] += ascii_lowercase[n-1]
        return
    
    def write_bbs(self, foldername, org_flag="_ORG", ino_flag="_INO"):
        """
        Write the clusters of the molgraph into the folder specified in the parameters.
        The names of the clusters written out will be those of the atomtypes of the vertices 
        in the topograph, if they are "vertex" building blocks (i.e. they have more than
        2 neighbours). If they are "edge" building blocks (they have exactly 2 neighbours), 
        their name will be the atomtypes of the vertex BBs they are connected to.
        Additionally, at the end of the filename, a string defined by the parameters org_flag
        and ino_flag will be appended to denote whether the BB is organic.
        
        :Parameters:
        - foldername: Name of the folder, in which the mfpx files of the building blocks 
          should be saved.
        - org_flag: String, which will be added to the filename if the BB is organic.
        - ino_flag: String, which will be added to the filename if the BB is inorganic.
        """
        self.tg_atypes_by_isomorphism()
        cv, ca = self.mg.get_cluster_atoms()
        bbs = []
        organicity = []
        for i, atoms in enumerate(ca):
            m = self.mg.mol.new_mol_by_index(atoms)
            bbs.append(m)
            if self.mg.cluster_organicity[i] == True:
                organicity.append(org_flag)
            else:
                organicity.append(ino_flag)
        # prepare vertex_bb_list (to translate indices of a list with 2-connected clusters to those of one without them)
        list2c = self.tg.midx_2conns
        vertex_bb_list = range(len(self.mg.clusters))
        for i in reversed(sorted(list2c)):
            del vertex_bb_list[vertex_bb_list.index(i)]
        # Use the atomtypes to identify "vertex" BBs
        atomtype_dict = {}
        for i, atype in enumerate(self.tg.mol.atypes):
            try:
                atomtype_dict[atype]
            except KeyError:
                atomtype_dict[atype] = []
            atomtype_dict[atype].append(vertex_bb_list[i])
        unique_bbs = []
        cluster_names = []
        for key in atomtype_dict.keys():
            unique_bbs.append(atomtype_dict[key])
            cluster_names.append(key)
        # determine which of the clusters are "edge" BBs and
        # Check to which vertex BBs the edge BBs are connected
        cluster_conn = self.mg.cluster_conn()
        neighbour_atypes = []
        for i in list2c:
            neighbour_atype = []
            for cc in cluster_conn[i]:
                for key in atomtype_dict.keys():
                    if cc in atomtype_dict[key]:
                        neighbour_atype.append(key)
                        break
            neighbour_atype = list(sorted(neighbour_atype))
            if neighbour_atype not in neighbour_atypes:
                neighbour_atypes.append(neighbour_atype)
                unique_bbs.append([i])
                cluster_names.append(neighbour_atype[0]+"-"+neighbour_atype[1])
        # Now write building blocks
        #print "============"
        #print unique_bbs
        #print atomtype_dict
        #print organicity
        #try:
        #    print neighbour_atypes
        #except:
        #    print "no 2-conns"
        #print "============"
        if not os.path.exists(foldername):
            os.mkdir(foldername)
        for n, i in enumerate(unique_bbs):
            m = bbs[i[0]]
            m.write(foldername+"/"+cluster_names[n]+organicity[i[0]]+".mfpx", "mfpx")
        return


