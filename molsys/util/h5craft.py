#   -*- coding: utf-8 -*-
###########################################
#
#   auxiliary module to forge hdf data v1.0
#   originally implemented for the
#       pydlpoly testing framework
#
#   by Roberto Amabile (2016)
#
#   PRE-ALPHA VERSION: USE AT YOUR OWN RISK
#       NO WARRANTY, NO LIABILITY
#   The author will use best judgement and
#       reasonable effort to keep the code
#       clear and tidy
#   NOT YET LICENSED
#
###########################################

import h5py, numpy as np

class DataReference:
    """This class builds the binary hdf reference file"""
    def __init__(self, fileref = None, modeaccess = None):
        self.fileref = fileref
        self.modeaccess = modeaccess
        self.h5file = h5py.File(self.fileref, self.modeaccess)
    """Recursively builds up the binary file in a tree of groups and datasets"""
    def build_rec_dataset(self, ref_dict, h5file = None):
        if h5file is None:
            h5file = self.h5file
        new_dict={}
        for key in ref_dict:
            value = ref_dict[key]
            if type(value) is dict:
                group = h5file.create_group(key)
                self.build_rec_dataset(value, group)
            else:
                h5file.create_dataset(key, data=value)
