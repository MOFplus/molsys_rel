import numpy
import string

def read(mol, f, delimiter=','):
    """
    Routine, which reads a plain file
    :Parameters:
        -f   (obj): plain file object
        -mol (obj): instance of a molclass
        -delimiter=',' (str): coordinate delimiter
    """
    splits = f.read().splitlines()
    xyz = [s.split(delimiter) for s in splits]
    mol.natoms = len(xyz)
    mol.xyz = numpy.array(xyz, dtype='float')
    mol.elems = ["x"]*mol.natoms
    mol.atypes = ["0"]*mol.natoms
    mol.set_empty_conn()
    mol.set_nofrags()
    return

def write(mol, fname):
    """
    Routine, which writes an xyz file
    :Parameters:
        -fname  (str): name of the plain file
        -mol    (obj): instance of a molclass
    """
    with open(fname,'w') as f:
        for i in range(natoms):
            f.write("%12.6f %12.6f %12.6f\n" % (mol.xyz[i,0], mol.xyz[i,1], mol.xyz[i,2]))
    return
