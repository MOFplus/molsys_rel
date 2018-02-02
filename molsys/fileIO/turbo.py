#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 18:17:13 2017

@author: johannes
"""

import numpy
import string
from molsys.util.units import angstrom, kcalmol

def read(mol, f, gradient = False, cycle = -1):
    if gradient:
        return read_gradfile(mol,f,cycle)
    coord = False
    xyz = []
    elems = []
    for line in f:
        sline = line.split()
        if sline[0][0] == "$":
            if sline[0] == "$coord": 
                coord = True
                continue
            else:
                coord = False
                continue
        if coord:
            xyz.append([float(i) for i in sline[:3]])
            elems.append(sline[3])
    f.close()
    mol.natoms = len(elems)
    mol.xyz = numpy.array(xyz)/angstrom
    mol.elems = elems
    mol.atypes = elems
    mol.set_empty_conn()
    mol.set_nofrags()
    return

def read_gradfile(mol, f, cycle):
    elems    = []
    xyz      = []
    grad     = []
    ncycle   = 0
    found    = False
    ### get how many cylces are in the file
    for line in f:
        sline = line.split()
        if sline[0] == "cycle": ncycle += 1
    f.seek(0)
    scycle = range(ncycle)[cycle]
    for line in f:
        sline = line.split()
        if sline[0] == "cycle" and int(sline[2])-1 == scycle:
            ### found read in 
            energy = float(sline[6])
            found  = True
        elif sline[0] == "cycle" and int(sline[2])-1 != scycle:
            found = False
        if found:
            if len(sline) == 4:
                ### coord info
                xyz.append([float(i) for i in sline[:3]])
                elems.append(sline[3])
            elif len(sline) == 3:
                ### grad info
                #grad.append([lambda a: float(a.replace("D","E"))(i) for i in sline[:3]])
                grad.append([float(i.replace("D","E")) for i in sline[:3]])
    f.close()
    mol.natoms = len(elems)
    mol.xyz = numpy.array(xyz)/angstrom
    mol.elems = elems
    mol.atypes = elems
    mol.set_empty_conn()
    mol.set_nofrags()
    mol.gradient = numpy.array(grad)*(angstrom/kcalmol)
    mol.energy = energy/kcalmol
    return


def write(mol, fname):
    f = open(fname, "w")
    f.write("$coord\n")
    c = mol.xyz*angstrom
    for i in range(mol.natoms):
        f.write("  %19.14f %19.14f %19.14f   %-2s\n" % 
                (c[i,0],c[i,1], c[i,2], mol.elems[i]))
    f.write("$end\n")
    f.close()
    return
