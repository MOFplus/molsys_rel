#!/usr/bin/env python
# -*- coding: utf-8 -*-
import string
from graph_tool import Graph
from graph_tool.topology import *
import numpy
import pdb
import molsys.topo as topo
import copy

class reader(object):

    def __init__(self):
        return

    def load_keys_from_file(self,fname):
        self.keys = {}
        f = open(fname, 'r')
        for line in f.xreadlines():
            sline = string.split(line)
            if len(sline)>0:
                if sline[0] == 'key':
                    dim = int(sline[1])
                    key = string.join(sline[2:])
                if sline[0] == 'id':
                    name = sline[1]
                    self.keys[name] = [dim,key]
        return

    def dump_keys_to_pickle(self,pname):
        return

    def load_keys_from_pickle(self,pname):
        return

    def __call__(self,name):
        return self.keys[name]


class lqg(object):

    def __init__(self, dim = 3):
        self.dim = dim
        return


    def read_systre_key(self, skey):
        skey = string.split(skey)
        self.nedges = len(skey)/5
        self.nvertices = 1
        self.edges = []
        self.labels = []
        for i in range(self.nedges):
            edge = map(int, skey[i*5:i*5+2])
            for j in edge:
                if j > self.nvertices: self.nvertices = j
            edge = list(numpy.array(edge)-1)
            label = map(int, skey[i*5+2:i*5+5])
            self.edges.append(edge)
            self.labels.append(label)
        return

    def get_lqg_from_topo(self,topo):
        # be careful not working for nets where an vertex is connected to itself
        self.nvertices = topo.get_natoms()
        self.nedges = 0
        self.edges = []
        self.labels = []
        for i in xrange(self.nvertices):
            for j,v in enumerate(topo.conn[i]):
                if v > i:
                    self.nedges += 1
                    self.edges.append([i,v])
                    #pdb.set_trace()
                    self.labels.append(list(topo.pconn[i][j]))
        return


    def build_lqg(self):
        self.molg = Graph(directed=True)
        self.molg.ep.label  = self.molg.new_edge_property("vector<double>")
        self.molg.ep.number = self.molg.new_edge_property("int")
        for i in xrange(self.nvertices):
            iv = self.molg.add_vertex()
        for i,e in enumerate(self.edges):
            ie = self.molg.add_edge(self.molg.vertex(e[0]),self.molg.vertex(e[1]))
            self.molg.ep.label[ie] = self.labels[i]
            self.molg.ep.number[ie] = i
        return

    def get_cyclic_basis(self):
        nbasevec = self.nedges - self.nvertices + 1
        self.nbasevec = nbasevec
        basis = numpy.zeros([nbasevec,self.nedges], dtype="int")
        self.molg.set_directed(False)
        tree = min_spanning_tree(self.molg)
        i = 0
        for e in self.molg.edges():
            if tree[e] == 0:
                self.molg.set_edge_filter(tree)
                vl, el = shortest_path(self.molg, self.molg.vertex(int(e.source())), self.molg.vertex(int(e.target())))
                self.molg.set_edge_filter(None)
                basis[i, self.molg.ep.number[e]] = 1
                for eb in el:
                    idx = self.molg.ep.number[eb]
                    ebt = self.get_edge_with_idx(idx)
                    if ebt.source() == e.source():
                        basis[i, self.molg.ep.number[eb]] = -1
                    else:
                        basis[i, self.molg.ep.number[eb]] = 1
                    e = ebt
                i += 1
        self.basis = basis
        self.molg.set_directed(True)
        return self.basis

    def get_ncocycles(self,n):
        self.molg.set_directed(False)
        cocycles = numpy.zeros([n, self.nedges])
        i = 0
        for v in self.molg.vertices():
            el = v.out_edges()
            for eb in el:
                idx = self.molg.ep.number[eb]
                ebt = self.get_edge_with_idx(idx)
                if ebt.source() == v:
                    cocycles[i, idx] = 1
                else:
                    cocycles[i, idx] = -1
            i+=1
            if i == n: break
        return cocycles

    def get_B_matrix(self):
        n = self.nedges - (self.nedges - self.nvertices +1)
        if n > 0: 
             cocycles = self.get_ncocycles(n)
             self.B = numpy.append(self.basis, cocycles, axis = 0)
        else:
            self.B = self.basis
        return self.B

    def get_alpha(self):
        vimg = []
        labels = numpy.array(self.labels)
        for i in range(numpy.shape(self.basis)[0]):
            img = numpy.sum(self.basis[i]* labels.T,axis = 1)
            vimg.append(img)
        for i in range(self.nedges-self.nbasevec):
            vimg.append([0,0,0])
        self.alpha = numpy.array(vimg)
        return self.alpha

    def get_fracs(self):
        self.fracs = numpy.dot(numpy.linalg.inv(self.B),self.alpha)
        return self.fracs


    def get_cell(self):
        if self.nbasevec == -5:
            self.cell = numpy.dot(self.basis, self.basis.T)
        else:
            k = numpy.zeros([self.nbasevec-3+self.nvertices-1,self.nedges])
            idx = self.find_li_vectors(self.alpha)
            latbase = self.alpha[idx]
            Lr = self.basis[idx]
            ### we need to orthonormalize the latbase ###
            L = numpy.zeros([3,self.nedges])
            olatbase = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
            for i in range(3):
                b = numpy.linalg.solve(latbase.T, olatbase[i,:])
                for j in range(3):
                    L[i,:]+= b[j]*Lr[j,:]
            counter = 0
            ### TODO: switsch to other basis to make it more beautiful
            for i in range(self.nbasevec):
                if i not in idx:
                    b = numpy.linalg.solve(latbase.T,self.alpha[i])
                    bb = numpy.zeros(self.nedges)
                    for j in range(3):
                        bb += b[j]*self.basis[idx[j]]
                    k[counter] = self.basis[i]-bb
                    counter += 1
            if self.nvertices > 1:
                k[self.nbasevec-3:,:] = self.get_ncocycles(self.nvertices-1)
            ### do projection of L ###
            S = numpy.dot(k,k.T)
            P = numpy.eye(self.nedges,self.nedges) - numpy.dot(k.T, 
                    numpy.dot(numpy.linalg.inv(S), k))
            self.cell = numpy.dot(L, numpy.dot(P,L.T))
        print self.cell
        import molsys.util.unit_cell as uc
        print uc.abc_from_vectors(self.cell)
        return self.cell

    def place_vertices(self, first = numpy.array([0,0,0])):
        frac_xyz = numpy.zeros([self.nvertices,3])
        frac_xyz[0,:] = first
        done = [0]
        counter = 0
        while len(done) != self.nvertices:
            for i,e in enumerate(self.edges):
                if self.labels[i] == [0,0,0]:
                    if ((e[0] in done) and (e[1] not in done)):
                        frac_xyz[e[1],:] = frac_xyz[e[0],:] + self.fracs[i,:]
                        done.append(e[1])
                    elif ((e[1] in done) and (e[0] not in done)):
                        frac_xyz[e[0],:] = frac_xyz[e[1],:] - self.fracs[i,:]
                        done.append(e[0])
            counter += 1
            if counter > 10: break
        if len(done) != self.nvertices:
            for i,e in enumerate(self.edges):
                if ((e[0] in done) and (e[1] not in done)):
                    frac_xyz[e[1],:] = frac_xyz[e[0],:] + self.fracs[i,:]
                    done.append(e[1])
                elif ((e[1] in done) and (e[0] not in done)):
                    frac_xyz[e[0],:] = frac_xyz[e[1],:] - self.fracs[i,:]
                    done.append(e[0])
        ### perhaps a flooring has to be performe
        self.frac_xyz = frac_xyz
        return self.frac_xyz

    def make_mol(self):
        t = topo()
        t.natoms = self.nvertices
        t.set_cell(self.cell)
        t.set_xyz_from_frac(self.frac_xyz)
        t.set_atypes(self.nvertices*['1'])
        t.set_empty_conn()
        for i,e in enumerate(self.edges):
            t.conn[e[0]].append(e[1])
            t.conn[e[1]].append(e[0])
        #t.wrap_in_box()
        t.set_elems_by_coord_number()
        t.write('test.xyz', 'txyz')
        return t


    def get_edge_with_idx(self, idx):
        for i in self.molg.edges():
            if self.molg.edge_index[i] == idx: return i
            #if self.molg.ep.number[i] == idx: return i

    def find_li_vectors(self,R):
        rank = numpy.linalg.matrix_rank(R)
        idx = []
        idx.append(0)
        for i in range(1,R.shape[0]):
            indep = True
            for j in idx:
                if i != j:
                    inner_product = numpy.dot( R[i,:], R[j,:] ) #compute the scalar product
                    norm_i = numpy.linalg.norm(R[i,:]) #compute norms
                    norm_j = numpy.linalg.norm(R[j,:])
                    if abs(inner_product - norm_j * norm_i) < 1e-4:
                        # vector i is linear dependent, iterate i
                        indep = False
                        break
            if indep == True:
                idx.append(i)
                if numpy.linalg.matrix_rank(R[idx]) != len(idx):
                    idx.pop()
                if len(idx)==rank: break
        return idx


    def find_li_vectors_old(self,dim, R):
        print R 
        r = numpy.linalg.matrix_rank(R) 
        index = numpy.zeros(r) #this will save the positions of the li columns in the matrix
        counter = 0
        index[0] = 0 #without loss of generality we pick the first column as linearly independent
        j = 0 #therefore the second index is simply 0
        for i in range(R.shape[0]): #loop over the columns
            print i, j 
            if i != j: #if the two columns are not the same
                inner_product = numpy.dot( R[i,:], R[j,:] ) #compute the scalar product
                norm_i = numpy.linalg.norm(R[i,:]) #compute norms
                norm_j = numpy.linalg.norm(R[j,:])
                #inner product and the product of the norms are equal only if the two vectors are parallel
                #therefore we are looking for the ones which exhibit a difference which is bigger than a threshold
                if abs(inner_product - norm_j * norm_i) > 1e-4:
                    print 'lala'
                    counter += 1 #counter is incremented
                    print counter, index
                    index[counter] = i #index is saved
                    j = i #j is refreshed
                    #do not forget to refresh j: otherwise you would compute only the vectors li with the first column!!
        R_independent = numpy.zeros((r, dim))
        i = 0
        print index
        #now save everything in a new matrix
        while( i < r ):
            R_independent[i,:] = R[index[i],:] 
            i += 1
        return R_independent


