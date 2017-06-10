#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 18:17:13 2017

@author: johannes
"""

import numpy
import string
from molsys.util.units import angstrom

def read(mol, f):
    coord = False
    xyz = []
    elems = []
    for line in f:
        sline = string.split(line)
        if sline[0][0] == "$":
            if sline[0] == "$coord": 
                coord = True
                continue
            else:
                coord = False
                continue
        if coord:
            xyz.append(map(float,sline[:3]))
            elems.append(sline[3])
    f.close()
    mol.natoms = len(elems)
    mol.xyz = numpy.array(xyz)/angstrom
    mol.elems = elems
    mol.atypes = elems
    mol.set_empty_conn()
    mol.set_nofrags()
    return


def write(mol, fname):
    f = open(fname, "w")
    f.write("$coord\n")
    c = mol.xyz*angstrom
    for i in xrange(mol.natoms):
        f.write("  %19.14f %19.14f %19.14f   %-2s\n" % 
                (c[i,0],c[i,1], c[i,2], mol.elems[i]))
    f.write("$end\n")
    f.close()
    return
