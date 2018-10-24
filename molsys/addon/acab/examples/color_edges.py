import molsys
from math import pi

### OPTIONS ###
net = "tbo" # as in mfpx file or in MOF+
supercell = [1,1,1]
scell = ''.join([str(s) for s in supercell])
eratio = [2,1]

### MAIN ###
# this is old
try:
    m = molsys.mol.from_file(net)
except IOError:
    try:
        import mofplus
    except ImportError as e:
        raise ImportError("Download net via https://www.mofplus.org/nets/net/%s" % net)
    api = mofplus.user_api()
    api.get_net(net)
    m = molsys.mol.from_file(net)
m.make_supercell(supercell)

# this is new
m.addon("acab")
m.acab.setup_model()
m.acab.setup_ecratio_per_vertex(eratio_per_vertex) # stricter
#m.acab.setup_ecratio(eratio_per_vertex) # looser
m.acab.setup_vcratio_per_edge([1]) # w/o this: edges are not connected!
#m.acab.setup_angle_btw_edges(color=1, theta=pi) #only straight angles
N = m.acab.cycle_loop(alpha=3, constr_vertex=False, rundir='run/'+scell)
print("Number of unequivalent colorings: %s" % N)

