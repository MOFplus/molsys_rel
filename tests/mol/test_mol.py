import pytest

import molsys
import molsys.addon
import molsys.fileIO.formats

from molsys.util.sysmisc import _makedirs

import sys
import os
from os.path import splitext

rundir = _makedirs("run")

fname = "1x1x1"
name = os.path.splitext(fname)[0]

m = molsys.mol.from_file(fname)

def test_manipulate_atoms():
    pass

def test_manipulate_coordinates():
    pass

scells = [
    [1,1,1],
    [2,1,1],
    [1,2,1],
    [1,1,2],
    [2,2,1],
    [1,2,2],
    [2,1,2],
]
@pytest.mark.slow
@pytest.mark.xpass(reason="too slow")
@pytest.mark.parametrize("scell", scells)
def test_manipulate_cell(scell):
    m.make_supercell(scell)

writefmts = molsys.fileIO.formats.write.keys()
@pytest.mark.parametrize("fmt", writefmts)
def test_write(fmt):
    m.write("%s/%s.%s" % (rundir, name, fmt))

readfmts = set(molsys.fileIO.formats.read.keys())
readfmts -= set(["mol2", "freq", "array", "cif", "plain", "cell", "aginp"])
@pytest.mark.parametrize("fmt", readfmts)
def test_read(fmt):
    m.read("%s/%s.%s" % (rundir, name, fmt))

@pytest.mark.parametrize("addon", molsys.addon.__all__)
def test_load_addon(addon):
    m.addon(addon)
