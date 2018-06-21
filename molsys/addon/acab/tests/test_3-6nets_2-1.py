import pytest

import molsys
import os

folder = "%s%snets" % (os.path.dirname(__file__), os.sep)

nets =     ["apo", "brk", "eea", "pyr", "qom", "rtl", "zzz"]#, "spn-z"]
ncolors_ = [   63,   -70,  -70,     43,   -70,     7,   -70]#,     -10]
maxcycle = [   70,    70,   70,     70,    70,    70,    70]#,      10]

@pytest.mark.parametrize("net,ncolors,maxcycle", zip(nets,ncolors_,maxcycle))
def test_nets(net, ncolors, maxcycle):
    m = molsys.mol.from_file("%s%s%s" % (folder, os.sep, net))
    m.addon("acab")
    m.acab.setup_model()
    m.acab.setup_ecratio_per_vertex([2,1])
    #m.acab.setup_vcratio_per_edge([1])
    N = m.acab.cycle_loop(Nmax=maxcycle, alpha=3, write=False)
    assert N == ncolors, "number of found colors %d must be equal to %d" \
        % (N, ncolors)
    assert N >= ncolors, "number of found colors %d must be equal to %d" \
        % (N, ncolors)
