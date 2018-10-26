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
# edge color ratio per vertex
ecratio = [2,1]
# run directory (with leading running index)
rundir = "%s%s%s_%s" % ('run', os.sep, "edges", scell)

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

### setup edge color ratio per each vertex (local)
m.acab.setup_ecratio_per_vertex(ecratio) # stricter
### setup edge color ratio altogether (global)
#m.acab.setup_ecratio(ecratio) # looser
### setup vertex color ratio as [1] => vertices are colored all the same!
### N.B.: w/o this: there is no vertex and edges are not connected!
m.acab.setup_vcratio_per_edge([1])  

### cycle loop of solutions
### N is the number of solutions
### N.B. constr_vertex=False to unconstraint colored vertices, otherwise the
###     loop finds just one solution! (there is only one way to color all the
###     vertices by just one color...)
N = m.acab.cycle_loop(constr_vertex=False, rundir=rundir, newrundir=False)
print("Number of unequivalent colorings: %s" % N)

