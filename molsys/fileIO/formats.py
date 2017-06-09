#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xyz, txyz, mfpx, cif, turbo

read = {
        'xyz':xyz.read,
        'txyz':txyz.read,
        'mfpx':mfpx.read,
        'cif':cif.read,
        'turbo':turbo.read}

write = {
        'xyz':xyz.write,
        'txyz':txyz.write,
        'mfpx':mfpx.write,
        'cif':cif.write,
        'turbo': turbo.write} 

#def read(mol, filename, fmt, **kwargs):
    #reads[fmt](mol,filename, **kwargs)

#def write(mol, filename, fmt, **kwargs):
    #writes[fmt](mol,filename, **kwargs)
