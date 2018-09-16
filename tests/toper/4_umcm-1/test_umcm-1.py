import pytest

import molsys
import molsys.util.toper as toper

# Original UMCM-1.cif file has broken disorder and unclear fractional occupancies
# conn_thresh = 1.2 is a rule of thumb to get rid of half-occ. C atoms... but you remove H atoms too!
#m = molsys.mol.from_file("UMCM-1.cif", conn_thresh=1.2) # RULE OF THUMB: actually we want H atoms

# Workaround: use another structure w/ same connectivity to get original UMCM-1 connectivity
m = molsys.mol.from_file("UMCM-1_dioxole.cif")
m.elems = ["h" if e=="o" and m.conn[i] == ["c"] else e for i,e in enumerate(m.elems)]

#scell = [2,2,2]
#strscell = "".join([str(i) for i in scell])
#m.make_supercell(scell)
#m.write("UMCM-1_"+strscell+".mfpx")
#m.write("UMCM-1_"+strscell+".txyz", pbc=False)

@pytest.mark.slow
def test_standard_bbs():
    tt = toper.topotyper(m)
    print(tt.get_net())
    tt.write_bbs("standard")

@pytest.mark.slow
def test_write_no_organicity_bbs():
    tt = toper.topotyper(m, split_by_org=False)
    print(tt.get_net())
    tt.write_bbs("not_split_by_org")
