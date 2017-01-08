#!/usr/bin/env python
# -*- coding: utf-8 -*-
import string
from graph_tool import Graph
from graph_tool.topology import *
import numpy

class lqg(object):

    def __init__(self):
        self.dim = 3
        return


    def read_systre_key(self, skey):
        assert self.dim == 3
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
                vl, el = shortest_path(self.molg, e.source(), e.target())
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

    def get_cyclic_lattice_basis(self):
        vlat = [
                numpy.array([0,0,1]),
                numpy.array([0,0,-1]),
                numpy.array([0,1,0]),
                numpy.array([0,-1,0]),
                numpy.array([1,0,0]),
                numpy.array([-1,0,0])
                ]
        ### first vector
        for i in vlat:
            idx = numpy.argwhere(numpy.all((self.alpha-i==0), axis = 1))
            print idx


    def get_cell(self):
        if self.nbasevec == 3:
            self.cell = numpy.dot(self.basis, self.basis.T)
        return self.cell


    def get_edge_with_idx(self, idx):
        for i in self.molg.edges():
            if self.molg.ep.number[i] == idx: return i





            


