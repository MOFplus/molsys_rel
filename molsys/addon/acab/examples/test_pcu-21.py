import pytest

from math import pi
import molsys

def base(supercell):
    m = molsys.mol.from_file("pcu")
    m.make_supercell(supercell)
    m.addon("acab")
    m.acab.setup_model()
    m.acab.setup_ecratio_per_vertex([2,1])
    m.acab.setup_vcratio_per_edge([1])
    m.acab.setup_angle_btw_edges(color=1, theta=pi)
    m.acab.cycle_loop(alpha=3, write=False, constr_vertex=False)

def test_111(): # it fails
    base([1,1,1])

def test_222():
    base([2,2,2])

def test_333():
    base([3,3,3])

@pytest.mark.slow
def test_444():
    base([4,4,4])
