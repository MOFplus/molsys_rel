#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xyz, txyz, mfpx, cif

read = {
        'xyz':xyz.read,
        'txyz':txyz.read,
        'mfpx':mfpx.read}

write = {
        'xyz':xyz.write,
        'txyz':txyz.write,
        'mfpx':mfpx.write,
        'cif':cif.write} 

def read(mol, filename, format, **kwargs):
    read[format](mol,filename, **kwargs)

def write(mol, filename, format, **kwargs):
    write[format](mol,filename, **kwargs)
