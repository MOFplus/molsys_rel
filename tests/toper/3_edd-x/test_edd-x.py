import pytest

import molsys
import molsys.util.toper as toper
from molsys.util.color import make_mol

import os

fnames = [  "1",  "2",  "3",  "4",  "5"]
nets   = ["rtl","eea","eea","apo","eea"]

@pytest.mark.parametrize("fname,net", zip(fnames,nets))
def test_compute_colors(fname, net):
    name = fname.rstrip('.cif')
    m = molsys.mol()
    mofpath = "%s%s%s.cif" % ("mofs", os.sep, fname)
    m.read(mofpath)

    tt = toper.topotyper(m)
    assert tt.get_net() == [net]

    folder = tt.write_bbs("%s%s%s" % ("colors", os.sep, name))
    tt.compute_colors()

    # TBI: it deserves new methods in topograph or
    # a specific colorgraph class
    ecolors = tt.tg.molg.ep.color.a
    n = make_mol(tt.tg.mol, alpha=3, ecolors=ecolors)
    n.write("%s%s%s.mfpx" % (folder, os.sep, name))
    n.write("%s%s%s.txyz" % (folder, os.sep, name), pbc=False)
