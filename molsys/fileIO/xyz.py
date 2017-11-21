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
    ### look how many cycles are in 
#    for line in f.readlines():
#        sline=line.split()
#        if len(sline)==1:
#            if int(sline)==natoms: ncycle+=1
#    ### now seek to cycle
#    cycle=range(ncycle)[cycle]
#    f.seek(natoms*cycle+1)
#    if len(fline)>1:
#        cellparams = map(float,fline[1:7])
#        mol.set_cellparams(cellparams)
    f.readline()
    xyz = numpy.zeros((natoms, 3))
    elements = []
    atypes = []
    for i in range(natoms):
        line = f.readline().split()
        elements.append(line[0].lower())
        atypes.append(line[0].lower())
        xyz[i,:] = list(map(float,line[1:4]))
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
#    if mol.periodic:
#        f.write("%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n\n" % tuple([mol.natoms]+mol.cellparams))
#    else:
    f.write("%d\n\n" % mol.natoms)
    for i in range(natoms):
        f.write("%2s %12.6f %12.6f %12.6f\n" % (mol.elems[i], mol.xyz[i,0], mol.xyz[i,1], mol.xyz[i,2]))
    f.close()
    return
