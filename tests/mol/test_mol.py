import pytest

import molsys
import sys
import os

m = molsys.mol.from_file(fname)

def test_manipulate_atoms():
    pass

def test_manipulate_coordinates():
    pass

def test_manipulate_cell():
    #supercell included
    m.make_supercell(scell)
    pass

def test_read():
    pass

def test_write():
    pass

def test_load_addon(addon):
    m.addon(addon)
