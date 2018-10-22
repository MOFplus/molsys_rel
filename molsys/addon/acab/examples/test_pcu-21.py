import pytest

from math import pi
import molsys

def axial(supercell):
    m = molsys.mol.from_file("pcu")
    m.make_supercell(supercell)
    m.addon("acab")
    m.acab.setup_model()
    m.acab.setup_ecratio_per_vertex([2,1])
    m.acab.setup_vcratio_per_edge([1])
    m.acab.setup_angle_btw_edges(color=1, theta=pi)
    scell = ''.join([str(s) for s in supercell])
    N = m.acab.cycle_loop(alpha=3, constr_vertex=False, rundir='run/'+scell)
    return N

def non_axial(supercell):
    m = molsys.mol.from_file("pcu")
    m.make_supercell(supercell)
    m.addon("acab")
    m.acab.setup_model()
    m.acab.setup_ecratio_per_vertex([2,1])
    m.acab.setup_vcratio_per_edge([1])
    scell = ''.join([str(s) for s in supercell])
    N = m.acab.cycle_loop(alpha=3, constr_vertex=False, rundir='run/'+scell)
    return N

@pytest.mark.xfail(raises=KeyError, reason="supercell too small (to be fixed)")
def test_axial_111():
    assert axial([1,1,1]) == 1, "Number of colorings is different from expected"

def test_axial_222():
    assert axial([2,2,2]) == 2, "Number of colorings is different from expected"

@pytest.mark.slow
def test_axial_333():
    assert axial([3,3,3]) == 2, "Number of colorings is different from expected"

@pytest.mark.xpass(reason="too slow")
def test_axial_444():
    assert axial([4,4,4]) == 4, "Number of colorings is different from expected"

@pytest.mark.xfail(raises=TypeError, reason="problems with periodic connectivity writing")
def test_non_axial_111():
    assert non_axial([1,1,1]) == 1, "Number of colorings is different from expected"

@pytest.mark.slow
def test_non_axial_222():
    assert non_axial([2,2,2]) == 41, "Number of colorings is different from expected"

