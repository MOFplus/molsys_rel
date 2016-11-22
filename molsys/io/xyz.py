import numpy
import string

def read(mol, fname):
    """
    Routine, which reads an xyz file
    :Parameters:
        -fname  (str): name of the xyz file
        -mol    (obj): instance of a molclass
    """
    f = open(fname, 'r')
    natoms = string.atoi(string.split(f.readline())[0])
    f.readline()
    xyz = numpy.zeros((natoms, 3))
    elements = []
    for i in range(natoms):
        line = string.split(f.readline())
        elements.append(string.lower(line[0]))
        xyz[i,:] = map(float,line[1:4])
    self.mol.natoms = natoms
    self.mol.xyz = xyz
    self.mol.elems = elements
    return

def write_xyz(self, mol, fname):
    """
    Routine, which writes an xyz file
    :Parameters:
        -fname  (str): name of the xyz file
        -mol    (obj): instance of a molclass
    """
    natoms = numpy.shape(xyz)[0]
    f = open(fname,"w")
    f.write("%d\n\n" % mol.natoms)
    for i in xrange(natoms):
        f.write("%s %12.6f %12.6f %12.6f\n" % (mol.elems[i], mol.xyz[i,0], mol.xyz[i,1], mol.xyz[i,2]))
    f.close()
    return
