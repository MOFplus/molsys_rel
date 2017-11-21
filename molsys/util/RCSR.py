#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
import string

class rcsr(object):

    def __init__(self):
        self._nets = {}
        return
    
    def read_arc(self,fname):
        f = open(fname, 'r')
        for line in f.xreadlines():
            sline = line.split()
            if len(sline)>0:
                if sline[0] == 'key':
                    dim = int(sline[1])
                    key = string.join(sline[2:])
                if sline[0] == 'id':
                    name = sline[1]
                    if name.count('*') > 0:
                        name = name.replace('*','s')
        if name in self._nets: 
            self._nets[name]['dim'] = dim
            self._nets[name]['key'] = key
        else:
            self._nets[name] = {'dim':dim, 'key':key}
        return

    def write_arc(self,fname):
        pass

    def read_cgd(self,fname):
        entries = open(fname, 'r').read().split('CRYSTAL')[1:]
        #entries = open(fname, 'r').read().split('CRYSTAL')[1:-1]
        #entries = open(fname, 'r').read()
        for i, e in enumerate(entries):
            self.parse_cgd(e)
        pass

    def parse_cgd(self, entry):
        dic = {}
        nodes = []
        edges = []
        lines = entry.split("\n")[1:]
        for i,l in enumerate(lines):
            sline = l.split()
            if sline[0] == 'NAME':
                name = sline[1]
                if name.count('*') > 0:
                    name = name.replace('*','s')
                dic['NAME'] = name
            elif sline[0] == 'GROUP':
                dic['GROUP'] = sline[1]
            elif sline[0] == 'CELL':
                dic['CELL'] = sline[1:7]
            elif sline[0] == 'NODE':
                nodes.append(map(float,sline[3:]))
            elif sline[0] == 'EDGE':
                edges.append(map(float,sline[1:]))
            elif sline[0] == 'END':
                break
        dic['nodes'] = nodes
        dic['edges'] = edges
        if name in self._nets.keys():
            self._nets[name]['cgd'] = dic
        else:
            self._nets[name] = {'cgd':dic}
        return
   
    def write_cgd(self, names, fname):
        def write_entry(f, name):
            try:
                cgd = self._nets[name]['cgd']
            except KeyError:
                return
            f.write('CRYSTAL\n')
            f.write('  NAME %s\n'  % name)
            f.write('  GROUP %s\n' % cgd['GROUP'])
            f.write('  CELL %s\n' % string.join(cgd['CELL']))
            for i in cgd['nodes']:
                f.write('  NODE %s\n' % string.join(i))
            for i in cgd['edges']:
                f.write('  EDGE %s\n' % string.join(i))
            f.write('END\n\n')
            return

        f = open(fname, 'w')
        for i,n in enumerate(names):
            write_entry(f, n)
        f.close()
        return

    def read_3dall(self, fname):
        txt = open(fname,'r').read().split('start')[1:-1]
        for i,t in enumerate(txt):
            self.parse_3dall(t)
        return

    def parse_3dall(self, txt):
        ndic = {}
#        print(txt)
        lines = txt.split('\n')[1:]
        # jump over line cotaining the id
        lines.pop(0)
        # get the netname
        name = lines.pop(0).split()[0]
        if name.count('*') > 0:
            name = name.replace('*','s')
        
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
