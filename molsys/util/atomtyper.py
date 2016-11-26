import numpy
import string
import copy
import math
import os

class atom:

    def __init__(self, element,connatoms):
        self.element = element
        connatoms.sort()
        self.connatoms = connatoms
        self.nconns = len(self.connatoms)
        self.type = ''
        return

    def __cmp__(self, other):
        if other.element != self.element:
            return cmp(self.element, other.element)
        else:
            return cmp(self.connatoms, other.connatoms)

    def get_diff(self,other):
        diff = [0,0,0]
        if __cmp__(self,other) == 0:
            return diff 
        else:
            if other.element != self.element:
                diff[0] = 1
                diff[1] = None
                diff[2] = None
                return diff
            diff[1] = other.nconns - self.nconns
            if diff[1] != 0:
                diff[2] = None
                return diff
            #diff[2] = len(set(other.connatoms)^set(self.connatoms))/float(self.nconns)
            diff[2] = len(set(other.connatoms)^set(self.connatoms))
            return diff


    def __repr__(self):
        rep =  "atom: element = %s, type = %s, bonded atoms = " % (self.element, self.type)
        rep += (self.nconns*"%s ") % tuple(self.connatoms)
        return rep

    def set_type(self,type):
        self.type = type
        return

    def get_type(self):
        return self.type

class atomtyper:

    def __init__(self, mol):
        self.mol = mol
        self.elements = mol.elems
        self.avail_e  = mol.get_elemlist()
        self.xyz = mol.xyz
        self.cnct = mol.conn
        self.natoms = mol.natoms
        self.atoms = []
        self.atypes = []
        self.setup_atoms()
        return

    def setup_atoms(self):
        for i in range(self.natoms):
            self.atoms.append(atom(self.elements[i], map(self.elements.__getitem__, self.cnct[i])))
        return

    def __call__(self,rules = 2): # depending on element
        """
        0 : only elements
        1 : element+coordnumber
        2 : element+coordnumber+bonded atoms
        """
        # create a rule dictionary
        self.atypes = []
        if isinstance(rules, int):
            rules = dict(zip(self.avail_e, len(self.avail_e) * [rules]))
        self.rules = rules
        # loop over all atoms
        for i,a  in enumerate(self.atoms):
            type = self.apply_rules(a)
            a.set_type(type)
            self.atypes.append(a.get_type())
        self.mol.atypes = self.atypes

    def metacall(self, short = True):
        """
        assigns atomtypes on the basis of previous assigned atomtypes,
        uses rules = 2 as basis. Can be called iteratively. 
        if short = False, it will list atomtypes of the next atoms.
        """

        types = list(set(self.atypes))
        for i in range(len(types)):
            type = types[i]
            tatoms = numpy.where(numpy.array(self.atypes)==type)[0].tolist()
            mtypes = []
            for j, t in enumerate(tatoms):
                l = []
                for k in self.cnct[t]:
                    l.append(self.atypes[k])
                mtypes.append(tuple(numpy.sort(l)))
            used = list(set(mtypes))
            nmtypes = len(used)
            if nmtypes > 1:
                ### assign new atomtypes ###
                for l in range(len(tatoms)):
                    oldtype = self.atypes[tatoms[l]]
                    if short:
                        newtype = oldtype + str(used.index(mtypes[l])+1)
                    else:
                        add = ""
                        for m in mtypes[l]:
                            add += "(" + m + ")"
                        newtype = oldtype + add
                    ### distribute new atomtypes ###
                    self.atypes[tatoms[l]] = newtype
        return


    def apply_rules(self, atom):
        rules_iml = [0,1,2]
        try:
            rule = self.rules[atom.element]
        except KeyError:
            print 'No rule found for element %s!' % atom.element
            exit()
        if rule not in rules_iml:
            print 'Rule %s not known' % rule
        type = self.apply_rule(atom, rule)
        return type

    def apply_rule(self,atom,rule):
        if rule == 0:
            type = atom.element
        elif rule == 1:
            type = atom.element+'%s' % str(atom.nconns)
        elif rule == 2:
            #type = (atom.element+'%s'+'_'+atom.nconns*"%s")\
            #        % tuple([str(atom.nconns)]+atom.connatoms)
            type = (atom.element+'%s'++atom.nconns*("%s"))\
                    % tuple([str(atom.nconns)]+atom.connatoms)
        return type
