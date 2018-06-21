import numpy
import string

def read(mol, arr, use_pconn=False, **kwargs):
    """
    Routine, which reads a coordinate array
    :Parameters:
        -arr    (ndarray): name of the coordinate array
        -mol        (obj): instance of a molclass
        -use_pconn (bool): True to set empty pconn
    """
    mol.natoms = len(arr)
    mol.xyz = numpy.array(arr, dtype='float')
    mol.elems = ["x"]*mol.natoms
    mol.atypes = ["0"]*mol.natoms
    mol.set_empty_conn()
    if use_pconn: mol.set_empty_pconn()
    mol.set_nofrags()
    return
