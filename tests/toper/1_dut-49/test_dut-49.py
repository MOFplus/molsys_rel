import pytest

import molsys
import molsys.util.toper as toper

m = molsys.mol()
m.read("DUT-49.mfpx")

@pytest.mark.slow
def test_standard_bbs():
    tt = toper.topotyper(m)
    print(tt.get_net())
    tt.write_bbs("standard")

@pytest.mark.slow
def test_no_organicity_bbs():
    tt = toper.topotyper(m, split_by_org=False)
    print(tt.get_net())
    tt.write_bbs("not_split_by_org")

