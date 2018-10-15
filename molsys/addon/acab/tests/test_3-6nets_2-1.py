import pytest

import molsys
import os

folder = "%s%snets" % (os.path.dirname(__file__), os.sep)

nets      = ["apo", "brk", "eea", "pyr", "qom", "rtl", "zzz"]#, "spn-z"]
ncolors   = [   63,   -70,  -70,     43,   -70,     7,   -70]#,     -10]
ncolors_  = [    3,     3,    0,      2,     0,     2,    10]#,       4]
ncolors__ = [    1,     1,    0,      2,     0,     2,    10]#,       4]
maxcycle  = [   70,    70,   70,     70,    70,    70,    70]#,      10]
#N.B.: spn-z is commented due to numerous symmetry operations (it's slow wrt. others!)
#remove the sharp and the closed bracket before the sharp ( ']#' ) to test

@pytest.mark.parametrize("net,ncolors,maxcycle", zip(nets,ncolors,maxcycle))
def test_nets(net, ncolors, maxcycle):
    m = molsys.mol.from_file("%s%s%s" % (folder, os.sep, net))
    m.addon("acab")
    m.acab.setup_model()
    m.acab.setup_ecratio_per_vertex([2,1])
    m.acab.setup_vcratio_per_edge([1])
    N = m.acab.cycle_loop(Nmax=maxcycle, alpha=3, constr_vertex=False,
        rundir="run/"+net)
    assert N == ncolors, "number of found colors %d must be equal to %d" \
        % (N, ncolors)

@pytest.mark.parametrize("net,ncolors,maxcycle", zip(nets,ncolors_,maxcycle))
def test_nets_loose_axis(net, ncolors, maxcycle):
    m = molsys.mol.from_file("%s%s%s" % (folder, os.sep, net))
    m.addon("acab")
    m.acab.setup_model()
    m.acab.setup_ecratio_per_vertex([2,1])
    m.acab.setup_vcratio_per_edge([1])
    vsele = [i for i,e in enumerate(m.conn) if len(e) == 6]
    m.acab.setup_angle_btw_edges(color=1, theta=2.6, vsele=vsele)
    N = m.acab.cycle_loop(Nmax=maxcycle, alpha=3, constr_vertex=False,
        rundir="run/"+net)
    assert N == ncolors, "number of found colors %d must be equal to %d" \
        % (N, ncolors)

@pytest.mark.parametrize("net,ncolors,maxcycle", zip(nets,ncolors__,maxcycle))
def test_nets_strict_axis(net, ncolors, maxcycle):
    m = molsys.mol.from_file("%s%s%s" % (folder, os.sep, net))
    m.addon("acab")
    m.acab.setup_model()
    m.acab.setup_ecratio_per_vertex([2,1])
    m.acab.setup_vcratio_per_edge([1])
    vsele = [i for i,e in enumerate(m.conn) if len(e) == 6]
    m.acab.setup_angle_btw_edges(color=1, theta=3, vsele=vsele)
    N = m.acab.cycle_loop(Nmax=maxcycle, alpha=3, constr_vertex=False,
        rundir="run/"+net)
    assert N == ncolors, "number of found colors %d must be equal to %d" \
        % (N, ncolors)
