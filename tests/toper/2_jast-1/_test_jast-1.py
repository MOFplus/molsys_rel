import pytest

import molsys
import molsys.util.toper as toper

import os

moffolder = "mofs"

def test_compute_parallel_colors():
    m = molsys.mol()
    mofpath = "%s%s%s" % (moffolder, os.sep, "JAST-1_parallel.txyz")
    m.read(mofpath)

    tt = toper.topotyper(m)
    assert tt.get_net() == ["pcu"]

    tt.compute_colors()
    folder = tt.write_bbs("parallel")

    from molsys.util.color import make_mol
    ecolors = tt.tg.molg.ep.color.a
    n = make_mol(tt.tg.mol, alpha=3, ecolors=ecolors)
    n.write("%s%s%s" % (folder, os.sep, "all.mfpx"))
    n.write("%s%s%s" % (folder, os.sep, "all.txyz"), pbc=False)

def test_compute_skew_colors():
    m = molsys.mol()
    mofpath = "%s%s%s" % (moffolder, os.sep, "JAST-1_skew.txyz")
    m.read(mofpath)

    tt = toper.topotyper(m)
    assert tt.get_net() == ["pcu"]

    tt.compute_colors()
    folder = tt.write_bbs("skew")

    from molsys.util.color import make_mol
    ecolors = tt.tg.molg.ep.color.a
    n = make_mol(tt.tg.mol, alpha=3, ecolors=ecolors)
    n.write("%s%s%s" % (folder, os.sep, "skew.mfpx"))
    n.write("%s%s%s" % (folder, os.sep, "skew.txyz"), pbc=False)

### TBI ###
#def test_compute_skews_colors():
#    m = molsys.mol()
#    m.read("test.mfpx")
#
#    tt = toper.topotyper(m)
#    print(tt.get_net())
#
#    tt.compute_colors()
#    tt.write_bbs("bbs_skew")
#
#    from molsys.util.color import make_mol
#    ecolors = tt.tg.molg.ep.color.a
#    n = make_mol(tt.tg.mol, alpha=3, ecolors=ecolors)
#    n.write("_test.mfpx")
#    n.write("_test.txyz", pbc=False)
