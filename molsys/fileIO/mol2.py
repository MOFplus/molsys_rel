#### 
import string
import numpy
def read(mol, f, delimiter=','):
    """
    Routine, which reads a mol2 formatted file
    :Parameters:
        -f   (obj): mol2 file object
        -mol (obj): instance of a molclass
        -delimiter=',' (str): coordinate delimiter
    """
    s = f.read().splitlines()
    bonds,atoms=0,0
    mol.xyz,mol.elems,mol.ctab = [],[],[]
    for i,line in enumerate(s): 
        #print(i,len(line),line)
        if len(line) != 0: stsp = line.split()
        if i==3:
            natoms = int(stsp[0])
            nbonds = int(stsp[1])
            atoms,bonds = natoms, nbonds
            continue
        if (bonds != 0) and (atoms == 0): 
            mol.ctab.append([int(stsp[0])-1,int(stsp[1])-1])
            bonds -= 1
            continue
        if atoms != 0: 
            mol.xyz.append(map(float, stsp[0:3]))
            mol.elems.append(stsp[3].lower())
            atoms -= 1
            continue
    import pdb; pdb.set_trace()
    mol.xyz = numpy.array(mol.xyz)
    mol.natoms = len(mol.xyz)
    mol.atypes = ["0"]*mol.natoms
    mol.set_conn_from_tab(mol.ctab)
    mol.set_nofrags()
    return

def write(mol, fname):
    """
    Write mol sytem in mol2 format 
    ### NOT YET IMPLEMENTED
    :Parameters:
        -fname  (str): name of the plain file
        -mol    (obj): instance of a molclass
    """
    print('NOT YET IMPLEMENTED!')
    return
