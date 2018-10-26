import molsys
from math import pi

import sys
import os

### OPTIONS ####################################################################
# net name as in filename or in MOF+ database
net = "pcu"
# number unit cell replica in the 3 directions of space
supercell = [2,2,2]
# string of supercell
scell = ''.join([str(s) for s in supercell])
# vertex color ratio
vcratio = [1,1]
# run directory (with leading running index)
rundir = "%s%s%s_%s" % ('run', os.sep, "vertices", scell)

### UTILITY ####################################################################
# script path
scrpath = os.path.realpath(os.path.dirname(sys.argv[0]))
# download error message
download_message = """Download net by yourself
    via https://www.mofplus.org/nets/net/%s
    at the bottom of the web page
    and move the downloaded net to:
        %s%s""" % (net, scrpath, os.sep)


### MAIN #######################################################################
#.. old part (molsys + mofplus) ................................................
try:
    m = molsys.mol.from_file(net)   # read net file
except IOError:
    # if there is no net file...
    try:
        import mofplus
    except ImportError as e:
        raise ImportError(download_message)
    # ...download net
    api = mofplus.user_api()
    api.get_net(net)
    m = molsys.mol.from_file(net)
# make net supercell
m.make_supercell(supercell)

#.. new part (acab addon) ......................................................
m.addon("acab")     # load addon as attribute of mol instance

m.acab.setup_model()    # setup model options

### setup vertex color ratio per each edge (local, for bipartite)
m.acab.setup_vcratio_per_edge(vcratio)  # stricter
### setup vertex color ratio altogether (global, for disorder)
#m.acab.setup_vcratio(vcratio) # looser
### N.B.: setup_ecratio is not needed to connect vertices here

### cycle loop of solutions
### N is the number of solutions
N = m.acab.cycle_loop(rundir=rundir, newrundir=False)
print("Number of unequivalent colorings: %s" % N)

