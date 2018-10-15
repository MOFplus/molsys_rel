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
    m.acab.cycle_loop(alpha=3, constr_vertex=False, rundir='run/'+scell)

def non_axial(supercell):
    m = molsys.mol.from_file("pcu")
    m.make_supercell(supercell)
    m.addon("acab")
    m.acab.setup_model()
    m.acab.setup_ecratio_per_vertex([2,1])
    m.acab.setup_vcratio_per_edge([1])
    m.acab.cycle_loop(alpha=3, write=False, constr_vertex=False)

@pytest.mark.xfail(raises=KeyError, reason="supercell too small (to be fixed)")
def test_axial_111():
    axial([1,1,1])

def test_axial_222():
    axial([2,2,2])

@pytest.mark.slow
def test_axial_333():
    axial([3,3,3])

@pytest.mark.xpass(reason="too slow")
def test_axial_444():
    axial([4,4,4])

def test_non_axial_111():
    non_axial([1,1,1])

@pytest.mark.slow
def test_non_axial_222():
    non_axial([2,2,2])

