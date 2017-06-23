import numpy
import string

def read(mol, arr, **kwargs):
    """
    Routine, which reads a coordinate array
    :Parameters:
        -arr (ndarray): name of the coordinate array
        -mol     (obj): instance of a molclass
    """
    mol.natoms = len(arr)
    mol.xyz = numpy.array(arr, dtype='float')
    mol.elems = ["x"]*mol.natoms
    mol.atypes = ["0"]*mol.natoms
    mol.set_empty_conn()
    mol.set_nofrags()
    return
