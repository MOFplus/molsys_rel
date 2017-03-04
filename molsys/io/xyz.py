import numpy
import string

def read(mol, f):
    """
    Routine, which reads an xyz file
    :Parameters:
        -fname  (str): name of the xyz file
        -mol    (obj): instance of a molclass
    """
    fline = string.split(f.readline())
    natoms = string.atoi(fline[0])
    if len(fline)>1:
        cellparams = map(string.atof,fline[1:7])
        mol.set_cellparams(cellparams)
    f.readline()
    xyz = numpy.zeros((natoms, 3))
    elements = []
    atypes = []
    for i in range(natoms):
        line = string.split(f.readline())
        elements.append(string.lower(line[0]))
        atypes.append(string.lower(line[0]))
        xyz[i,:] = map(float,line[1:4])
    mol.natoms = natoms
    mol.xyz = numpy.array(xyz)
    mol.elems = elements
    mol.atypes = atypes
    mol.set_empty_conn()
    mol.set_nofrags()
    return

def write(mol, fname):
    """
    Routine, which writes an xyz file
    :Parameters:
        -fname  (str): name of the xyz file
        -mol    (obj): instance of a molclass
    """
    natoms = mol.natoms 
    f = open(fname,"w")
    if mol.periodic:
        f.write("%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n\n" % tuple([mol.natoms]+mol.cellparams))
    else:
        f.write("%d\n\n" % mol.natoms)
    for i in xrange(natoms):
        f.write("%2s %12.6f %12.6f %12.6f\n" % (mol.elems[i], mol.xyz[i,0], mol.xyz[i,1], mol.xyz[i,2]))
    f.close()
    return
