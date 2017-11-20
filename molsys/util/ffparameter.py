#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np

import logging
logger = logging.getLogger('molsys.ff')

class potentials(dict):
    """
    Class to store the parameter values, multiple ff objects can use the same par instance
    in order to perform multistruc fits
    """
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
        return

    def attach_variables(self):
        if hasattr(self, 'variables') == False:
            self.variables = varpars()

class varpar(object):

    """
    Class to hold information of parameters marked as variable in order to be fitted.
    """

    def __init__(self, par, name, val = 1.0, range = [0.0,2.0], bounds = ["h","i"]):
        assert len(bounds) == 2
        assert bounds[0] in ["h", "i", "z"] and bounds[1] in ["h", "i"]
        self._par     = par
        self.name    = name
        self._val     = val
        self.range   = range
        self.pos     = []
        self.bounds   = bounds

    def __repr__(self):
        """
        Returns objects name
        """
        return self.name

    def __call__(self, val = None):
        """
        Method to set a new value to the varpar oject and write this 
        into the ff.par dictionary.
        :Parameters:
            - val(float): new value of the varpar object
        """
        if val != None: self.val = val
        for i,p in enumerate(self.pos):
            ic, pottype, parindex  = p
            self._par[ic][pottype][1][parindex] = self.val
        return

    @property
    def val(self):
        return self._val

    @val.setter
    def val(self,val):
        assert (type(val) == float) or (type(val) == np.float64) or (val[0] == "$") 
        self._val = val
        return

class varpars(dict):

    """
    Class inherited from dict, that holds all varpar objects
    """

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
        return

    def __setitem__(self,k,v):
        assert type(v) == varpar
        # loop over all items and check if the variable is already in there 
        if k in self.keys():
            k += 'd'
            v.name = k
            return self.__setitem__(k,v)
        else:
            super(varpars,self).__setitem__(k,v)
            return k

    @property
    def ranges(self):
        ranges = []
        for k,v in self.items():
            ranges.append(v.range)
        return np.array(ranges)

    @ranges.setter
    def ranges(self, ranges):
        assert len(ranges) == len(self.keys())
        for i,v in enumerate(self.values()):
            v.range = ranges[i]

    @property
    def vals(self):
        vals   = []
        for k,v in self.items():
            vals.append(v.val)
        return vals

    def cleanup(self):
        """
        Method to delete all unsused varpar objects
        """
        rem = []
        for k,v in self.items():
            if len(v.pos) == 0:
                logger.warning("varpar %s is not used --> will be deleted!" % k)
                rem.append(k)
        for k in rem: del[self[k]]
        return

    @property
    def varpots(self):
        varpots = []
        for k,v in self.items():
            for i in range(len(v.pos)):
                varpot = (v.pos[i][0], v.pos[i][1])
                if varpot not in varpots: varpots.append(varpot)
        return varpots

    def varpotnums(self,ff):
        """
        Property which gives a dictionary telling in how much terms a
        varpot is involved
        """
#        ff = self.values()[0]._par
        varpots = self.varpots
        varpotnums = {}
        for i in varpots: varpotnums[i] = 0
        ics = ["bnd", "ang", "dih", "oop"]
        for ic in ics:
            for pi in ff.parind[ic]:
                for p in pi:
                    if (ic,p) in varpots: varpotnums[(ic,p)]+=1
        return varpotnums


    def __call__(self, vals = None):
        """
        Method to write new values to the varpar objects and in the 
        ff.par dictionary
        :Parameters:
            -vals(list of floats): list holding the new parameters 
        """
        if type(vals) == type(None): vals = len(self)*[None]
        assert len(vals) == len(self)
        for i,v in enumerate(self.values()): v(vals[i])
        return


