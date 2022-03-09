""" Group Addon

The group addon is a simple tool to maintain groups of (identical) subgroups.

We can generate such subgroups from the molecules or the fragemnts addon or by giving indices.
Since there are different ways to handle and generate such subgroups on the molsys level, the groups
serve as the single interface to control groups in pylmps/lammps

Note: in lammps, groups refer to a "set of atoms", but these groups contain multiple subsets of atoms
and are closer to chunks in lammps.

Note: atom indices in molsys start with 0 whereas in lammps atom indexing starts with 1

A group does not have to consist of equal subsets (like guest molecules or MOF SBUs etc.) but they are
meant to be used that way. The tools to make such groups will enforce that but if you know what you are doing 
feel free to generate anything by indices
"""

import numpy as np

class groups(object):

    def __init__(self, mol):
        """init groups addon

        for the moment simply pass the mol instance the addon was loaded

        Args:
            mol (molsys.mol): mol object
        """
        self._mol = mol
        self._groups = {}
        self._modes = {}
        return

    def add_group(self, name, mode, idx):
        self._groups[name] = idx
        self._modes[name] = mode

    def get_N(self, name):
        return len(self._groups[name])

    def get_idx(self, name):
        return self._groups[name]
    
    def get_mode(self, name):
        return self._modes[name]

    @property
    def groupnames(self):
        return list(self._groups.keys())

    def get_flat_string(self, name, incr = 1):
        idx = self._groups[name]
        s = ""
        for m in idx:
            ml = np.array(m) + incr
            s += (len(m)*"%d ") % tuple(ml)
        return s



