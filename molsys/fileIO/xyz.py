import numpy
import string

def read(mol, f, cycle = 0):
    """
    Routine, which reads an xyz file
    :Parameters:
        -f   (obj): xyz file object
        -mol (obj): instance of a molclass
    """
    ncycle=0
    fline = f.readline().split()
    natoms = int(fline[0])
    line = f.readline()
    periodic = False
    # check for keywords used in extended xyz format
    # in the moment only the lattice keyword is implemented
    if "Lattice=" in line:
        periodic = True
        lattice = [float(i) for i in line.rsplit('"',1)[0].rsplit('"',1)[-1].split() if i != '']
        cell = numpy.array(lattice).reshape((3,3))
    xyz = numpy.zeros((natoms, 3))
    elements = []
    atypes = []
    for i in range(natoms):
        line = f.readline().split()
        elements.append(line[0].lower())
        atypes.append(line[0].lower())
        xyz[i,:] = list(map(float,line[1:4]))
    mol.natoms = natoms
    if periodic: mol.set_cell(cell)
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
        f.write("%d\n" % mol.natoms)
        f.write('Lattice="%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n"' % 
            tuple(mol.cell.ravel()))
    else:
        f.write("%d\n\n" % mol.natoms)
    for i in range(natoms):
        f.write("%2s %12.6f %12.6f %12.6f\n" % (mol.elems[i], mol.xyz[i,0], mol.xyz[i,1], mol.xyz[i,2]))
    f.close()
    return
