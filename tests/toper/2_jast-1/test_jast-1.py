import pytest

import molsys
import molsys.util.toper as toper
from molsys.util.color import make_mol

import glob
import os
from molsys.util.sysmisc import _makedirs

moffolder = "mofs"

frameworks = glob.glob("weave/mofs/333/*.mfpx")
@pytest.mark.parametrize("mofpath", frameworks)
def test_compute_colors(mofpath):
    mofname = mofpath.rstrip(".mfpx").split(os.sep)[-1]
    m = molsys.mol.from_file(mofpath)

    tt = toper.topotyper(m)
    assert tt.get_net() == ["pcu"], "JAST-1 has pcu net: something went wrong!"

    tt.compute_colors()
    folder = tt.write_bbs("%s%s%s" % ("run", os.sep, mofname), index_run=True)

    ecolors = tt.tg.molg.ep.color.a
    n = make_mol(tt.tg.mol, alpha=3, ecolors=ecolors)
    n.write("%s%s%s.mfpx" % (folder, os.sep, mofname.lstrip("mof_")))
    n.write("%s%s%s.txyz" % (folder, os.sep, mofname.lstrip("mof_")), pbc=False)
