import molsys
import os
from math import pi

### OPTIONS ####################################################################
#.. new options ................................................................
# angle between second-colored edges
#theta = pi # only straight angles == axes
#theta = 3.1
#theta = 3.0
theta = 2.6
#theta = 0

# sense of the angle constraint
sense="min" # minimum theta allowed
#sense="max" # maximum theta allowed
#sense="close" # close theta allowed

# tolerance
eps=1e-3

#.. old options ................................................................
net = "brk" # must be in the running directory
supercell = [2,2,2]
scell = ''.join([str(s) for s in supercell]) # just a string
ecratio = [2,1] # edge color ratio
rundir = "%s%s%s_%s_%s" % ("run", os.sep, "axes", net, scell)

### MAIN #######################################################################
#.. old part (settings and constraints) ........................................
m = molsys.mol.from_file(net)
m.addon("acab")
m.acab.setup_model()
m.acab.setup_ecratio_per_vertex(ecratio) # applies to EACH vertices (3-c & 6-c)
m.acab.setup_vcratio_per_edge([1]) # keeps edge connected w/ a 0-colored vertex!

#.. new part (new constraint) ..................................................
### select 6-connected vertices which angle constraint holds for
### N.B.: 3-connected vertices are left free (still edge color ratio holds)
vsele = [i for i,e in enumerate(m.conn) if len(e) == 6] # list of vertex indices
# theta in radiants
# color is the order, e.g. 1 refers to the SECOND index (Python starts from 0)
#   so that if ecratio = [2,1] then it applies to the 1-ratioed edges
m.acab.setup_angle_btw_edges(color=1, theta=theta, sense=sense, eps=eps,
    vsele=vsele)

#.. old part (loop) ............................................................
N = m.acab.cycle_loop(constr_vertex=False, rundir=rundir, newrundir=False)
print("Number of unequivalent colorings: %s" % N)

