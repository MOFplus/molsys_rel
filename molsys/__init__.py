### 
import molsys
import numpy
molsys = molsys.molsys


images = []
for x in xrange(-1,2):
    for y in xrange(-1,2):
        for z in xrange(-1,2):
            images.append([x,y,z])
images = numpy.array(images,"d")
