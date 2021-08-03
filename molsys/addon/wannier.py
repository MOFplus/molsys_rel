#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import numpy
import copy
import string

from molsys.addon import base
from molsys.util.constants import debye, angstrom

eA2debye = 1./(debye/angstrom)

def read_until(f,s):
    while True:
        line = f.readline()
        if len(line)==0:
            return ""
        if string.find(line,s) >= 0:
            return line.strip()


class wannier(base):


    def __init__(self, mol, specifier = "x", dcharges = {"na":9., "o":6., "n":5.}):
        super(wannier,self).__init__(mol)
        self._specifier = specifier
        self._centers = []
        self._shells = []
        self._atoms = []
        self._charges = numpy.zeros([self._mol.natoms])
        self._dcharges = dcharges
        self._dcharges[self._specifier]=-2.
        return

    def portion_centers(self):
        atoms = self.atoms
        atom2centers = [[] for i in range(self.natoms)]
        for i,idx in enumerate(self.centers):
            a = self._mol.xyz[atoms]-self._mol.xyz[idx]
            if self._mol.periodic:
                #if self._mol.bcond <= 2:
                #    cell_abc = self._mol.cellparams[:3]
                #    a -= cell_abc * numpy.around(a/cell_abc)
                #if self._mol.bcond == 3:
                frac = numpy.dot(a, self._mol.inv_cell)
                frac -= numpy.around(frac)
                a = numpy.dot(frac, self._mol.cell)
            dist = numpy.sqrt((a*a).sum(axis=1)) # distance to all atoms
            # sort distances
            atom2centers[numpy.argsort(dist)[0]].append(idx)
        self._atom2centers = atom2centers
        return

    def find_shells(self):
        self._shells  = []
        for i,e in enumerate(self._mol.elems):
            if e[0].lower()=="x":
                self._shells.append(i)

    def prep_molecules(self):
        self._mol.addon("molecules")
        self._mol.molecules()
        #self.dipoles = {}
        #diptypes = {}
        #for m in self._mol.molecules.mols:
        #    if len(m) not in self.dipoles.keys():
        #        diptypes[len(m)]=1
        #    else:
        #        diptypes[len(m)]+=1
        #for i,k in diptypes.items():
        #    self.dipoles[i]=numpy.zeros((k))


    def calc_mol_dipoles(self):
        #self._mol.addon("molecules")
        #self._mol.molecules()
        # distinct species
        results = []
        self.portion_centers()
        for m in self._mol.molecules.mols:
            if len(m)>1:
            #if len(m)==1 and self._mol.elems[m[0]].lower() == "na":
                if len(self._atom2centers)!=0:
                    belongs = []
                    for i in m:
                        belongs += self._atom2centers[i]
                    results.append(self.calc_mol_dipole(m+belongs))
                else:
                    results.append(self.calc_mol_dipole(m))
            #if len(m) > 1: print(m,self.calc_mol_dipole(m))
        #print(results)
        return results
        

    @property
    def ncenters(self):
        return len(self.centers)


    @property
    def centers(self):
        centers = []
        for i,e in enumerate(self._mol.elems):
            if e == self._specifier: centers.append(i)
        self._centers = centers
        return self._centers

    @property
    def natoms(self):
        return len(self._atoms)

    @property
    def atoms(self):
        atoms = []
        for i,e in enumerate(self._mol.elems):
            if e != self._specifier: atoms.append(i)
        self._atoms = atoms
        return self._atoms
    
    @property
    def charges(self):
        for i, e in enumerate(self._mol.elems):
            self._charges[i] = self._dcharges[e]
        return self._charges

    def calc_mol_dipole(self, idx):
        molidx = []
        for i in idx:
            if i not in self._centers and i not in self._shells:
                molidx.append(i)
        com = self._mol.get_com(molidx)
        #print (self._charges[idx])
        #print (com)
        #print (self._mol.xyz[idx])
        # for periodic structures it has to be checked first if
        # all are centers are in the same image relative to the 
        # com
        a = self._mol.xyz[idx]-com
        if self._mol.periodic:
            #if self._mol.bcond <=2:
            #    cell_abc = self._mol.cellparams[:3]
            #    a[:,:] -= cell_abc*numpy.around(a/cell_abc)
            #if self._mol.bcond <=3:
            frac = numpy.dot(a, self._mol.inv_cell)
            a[:,:] -= numpy.dot(numpy.around(frac),self._mol.cell)
        d = a*self._charges[idx][:,numpy.newaxis]
        d = numpy.sum(d, axis = 0)*eA2debye
        #print (d)
        #m = numpy.linalg.norm(d)
        return d

    def calc_total_monopole(self):
        return numpy.sum(self._charges)


    def read_charges_from_cp2k(self,f, ctype = "MPA"):
        if ctype == "MPA":
            if len(read_until(f, "Mulliken Population Analysis")) == 0:
                raise IOError("No charges found!")
        elif ctype == "HPA":
            if len(read_until(f, "Hirshfeld Charges")) == 0:
                raise IOError("No charges found!")
        charge = []
        for i in range(2): f.readline()
        for i in range(self._mol.natoms):
            sline = f.readline().split()
            if ctype=="MPA":
                charge.append(float(sline[4]))
            else:
                charge.append(float(sline[5]))
        self._charges=numpy.array(charge)
        return


    def read_charges_from_tm(self,f,ctype="ESP"):
        if ctype == "ESP":
            if len(read_until(f, "charges resulting from fit:")) == 0:
                raise IOError("No charges found!")
        elif ctype == "MPA":
            if len(read_until(f, "atomic populations from total density:")) == 0:
                raise IOError("No charges found!")
        charge = []
        if ctype == "MPA":
            for i in range(2): f.readline()
        else:
            for i in range(3): f.readline()
        for i in range(self._mol.natoms):
            sline = f.readline().split()
            if ctype == "MPA":
                charge.append(float(sline[1]))
            else:
                charge.append(float(sline[3]))
        self._charges=numpy.array(charge)
