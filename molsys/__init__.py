import numpy as np
from util import unit_cell 
from util import elems as elements
from util import rotations
from io import formats


images = []
for x in xrange(-1,2):
    for y in xrange(-1,2):
        for z in xrange(-1,2):
            images.append([x,y,z])
images = np.array(images,"d")
