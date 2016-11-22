import numpy
import string
import unit_cell

def read(self, mol, fname):
    """
    Routine, which reads an mfpx file
    :Parameters:
        -fname  (str): name of the txyz file
        -mol    (obj): instance of a molclass
    """
    elems = []
    xyz = []
    atypes = []
    conn = []
    type = 'xyz'
    f = open(fname, "r")
    ### read header ###
    lbuffer = string.split(f.readline())
    stop = False
    while not stop:
        if lbuffer[0] != '#':
            natoms = int(lbuffer[0])
            stop = True
        else:
            keyword = lbuffer[1]
            if keyword == 'type':
                type = lbuffer2
            if keyword == 'cell':
                pass

    lbuffer = string.split(f.readline())
    natoms = string.atoi(lbuffer[0])
    if len(lbuffer) > 1 and lbuffer[1] != 'molden':
        boundarycond = 3
        if lbuffer[1] == "#":
            # read full cellvectors
            celllist = map(string.atof,lbuffer[2:11])
            cell = numpy.array(celllist)
            cell.shape = (3,3)
            cellparams = unit_cell.abc_from_vectors(cell)
        else:
            cellparams = map(string.atof, lbuffer[1:7])
            cell = unit_cell.vectors_from_abc(cellparams)
        if ((cellparams[3]==90.0) and (cellparams[4]==90.0) and (cellparams[5]==90.0)):
            boundarycond=2
            if ((cellparams[0]==cellparams[1])and(cellparams[1]==cellparams[2])and\
                (cellparams[0]==cellparams[2])):
                    boundarycond=1
    for i in xrange(natoms):
        lbuffer = string.split(f.readline())
        xyz.append(map(string.atof, lbuffer[2:5]))
        elems.append(string.lower(lbuffer[1]))
        t = lbuffer[5]
        atypes.append(t)
        conn.append((numpy.array(map(string.atoi, lbuffer[6:]))-1).tolist())
    # done: wrap up
    xyz = numpy.array(xyz)
    mol.elems = elems
    mol.xyz = xyz
    mol.atypes = atypes
    mol.conn = conn
    if 'cell' in locals():
        mol.cell = cell
        mol.cellparams = cellparams
        mol.bcond = boundarycond
    return 

def write_tinker_xyz(self, mol, fname):
    """
    Routine, which writes an txyz file
    :Parameters:
        -fname  (str): name of the txyz file
        -mol    (obj): instance of a molclass
    """
    elems  = mol.elements
    atypes = mol.atypes
    xyz    = mol.xyz
    cnct   = mol.conn
    natoms = mol.natoms
    cellparams = mol.cellparams
    f = open(fname, 'w')
    if type(cellparams) != type(None):
        f.write("%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" % tuple([natoms]+cellparams))
    else:
        f.write("%5d \n" % natoms)
    for i in xrange(natoms):
        line = ("%3d %-3s" + 3*"%12.6f" + " %5s") % \
            tuple([i+1]+[elems[i]]+ xyz[i].tolist() + [atypes[i]])
        conn = (numpy.array(cnct[i])+1).tolist()
        if len(conn) != 0:
            line += (len(conn)*"%7d") % tuple(conn)
        f.write("%s \n" % line)
    f.close()
    return
