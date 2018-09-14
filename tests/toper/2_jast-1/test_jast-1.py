import pytest

import molsys
import molsys.util.toper as toper

@pytest.mark.slow
def test_compute_parallel_colors():
    m = molsys.mol()
    m.read("JAST-1_parallel.txyz")

    tt = toper.topotyper(m)
    print(tt.get_net())

    tt.compute_colors()
    tt.write_bbs("bbs_parallel")

    from molsys.util.color import make_mol
    ecolors = tt.tg.molg.ep.color.a
    n = make_mol(tt.tg.mol, alpha=3, ecolors=ecolors)
    n.write("all.mfpx")
    n.write("all.txyz", pbc=False)

def test_compute_skew_colors():
    m = molsys.mol()
    m.read("JAST-1_skew.txyz")

    tt = toper.topotyper(m)
    print(tt.get_net())

    tt.compute_colors()
    tt.write_bbs("bbs_skew")

    from molsys.util.color import make_mol
    ecolors = tt.tg.molg.ep.color.a
    n = make_mol(tt.tg.mol, alpha=3, ecolors=ecolors)
    n.write("skew.mfpx")
    n.write("skew.txyz", pbc=False)
