#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy

class rcsr(object):

    def __init__(self):
        self._nets = {}
        return
    
    def read_arcs(self, fname):
        pass
   
    def read_3dall(self, fname):
        txt = open(fname,'r').read().split('start')[1:-1]
        for i,t in enumerate(txt):
            self.parse_3dall(t)
        return

    def parse_3dall(self, txt):
        ndic = {}
        lines = txt.split('\n')[1:]
        # jump over line cotaining the id
        lines.pop(0)
        # get the netname
        name = lines.pop(0).split()[0]
        ndic['name'] = name
        ndic['embed_type'] = lines.pop(0)

        nsymbols = int(lines.pop(0).split('!')[0])
        ndic['symbols'] = [lines.pop(0) for i in range(nsymbols)]

        nnames = int(lines.pop(0).split('!')[0])
        ndic['knownas'] = [lines.pop(0) for i in range(nnames)]
        nnames = int(lines.pop(0).split('!')[0])
        ndic['knownas'] += [lines.pop(0) for i in range(nnames)]
        nkeys = int(lines.pop(0).split('!')[0])
        ndic['keywords'] = [lines.pop(0) for i in range(nkeys)]
        nrefs = int(lines.pop(0).split('!')[0])
        ndic['refs'] = [lines.pop(0) for i in range(nrefs)]
        t = lines.pop(0).split()
        ndic['sg_name'] = t[0]
        ndic['sg_number'] = t[1]
        ndic['cell'] = numpy.array(map(float,lines.pop(0).split()))
        #self.make_cellvec()
        
        nverts = int(lines.pop(0))
        
        ndic['symbolic']   = []
        ndic['wyckoff']     = []
        ndic['symmetry']    = []
        ndic['order']       = []
        ndic['node']        = []
        ndic['node_coordination'] = []
        ndic['cs']          = []
        ndic['vs']          = []

        for i in range(nverts):
            ndic['node_coordination'].append(int(lines.pop(0).split()[-1]))
            ndic['node'].append(map(float,lines.pop(0).split()))
            ndic['symbolic'].append(lines.pop(0).split()[0])
            ndic['wyckoff'].append(lines.pop(0).split()[0])
            ndic['symmetry'].append(lines.pop(0).split()[0])
            ndic['order'].append(int(lines.pop(0)))
    
        nedges = int(lines.pop(0))
        
        ndic['center_symbolic']   = []
        ndic['center_wyckoff']    = []
        ndic['center_symmetry']   = []
        ndic['edge_center']       = []
        for i in range(nedges):
            temp = lines.pop(0)
            ndic['edge_center'].append(map(float,lines.pop(0).split()))
            ndic['center_symbolic'].append(lines.pop(0))
            ndic['center_wyckoff'].append(lines.pop(0))
            ndic['center_symmetry'].append(lines.pop(0))
        # jump over the next 5 lines
        # read coord seqences and vertex symbols
        for i in range(5): lines.pop(0)
        for i in range(nverts):
            ndic['cs'].append(map(int, lines.pop(0).split())[:-1])
        for i in range(nverts):
            ndic['vs'].append(lines.pop(0).split()[0])
        # put ndic into the overall _nets dictionaray
        if self._nets.keys().count(name) == 0:
            self._nets[name] = ndic
        else:
            self._nets[name].update(ndic)
        return 
