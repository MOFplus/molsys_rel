#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xyz, txyz, mfpx, cif, plain, array, turbo, lammpstrj, mol2

read = {
        'xyz':xyz.read,
        'txyz':txyz.read,
        'mfpx':mfpx.read,
        'cif':cif.read,
        'plain':plain.read,
        'array':array.read,
        'turbo':turbo.read,
        'mol2':mol2.read}

write = {
        'xyz':xyz.write,
        'txyz':txyz.write,
        'mfpx':mfpx.write,
        'cif':cif.write,
        'plain':plain.write,
        'lammpstrj':lammpstrj.write,
        'turbo':turbo.write}

#def read(mol, filename, fmt, **kwargs):
    #reads[fmt](mol,filename, **kwargs)

#def write(mol, filename, fmt, **kwargs):
    #writes[fmt](mol,filename, **kwargs)
