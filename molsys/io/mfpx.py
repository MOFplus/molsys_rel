import numpy
import string
import unit_cell
import txyz
import logging

def read(self, mol, fname):
    """
    Routine, which reads an mfpx file
    :Parameters:
        -fname  (str): name of the txyz file
        -mol    (obj): instance of a molclass
    """
    ftype = 'xyz'
    f = open(fname, "r")
    ### read header ###
    lbuffer = string.split(f.readline())
    stop = False
    while not stop:
        if lbuffer[0] != '#':
            mol.natoms = int(lbuffer[0])
            stop = True
        else:
            keyword = lbuffer[1]
            if keyword == 'type':
                ftype = lbuffer2
            elif keyword == 'cell':
                mol.cellparams = map(string.atof,lbuffer[2:8])
                mol.cell = unit_cell.vectors_from_abc(mol.cellparams)
            elif keyword == 'cellvect'
                celllist = map(string.atof,lbuffer[2:11])
                cell = numpy.array(celllist)
                cell.shape = (3,3)
                mol.cell = cell
                mol.cellparams = unit_cell.abc_from_vectors(mol.cell)
            elif keyword == 'bbcenter':
                mol.centerpoint = lbuffer[2]
                if mol.centerpoint == 'special':
                    mol.special_center_point = np.array(map(float,lbuffer[3:6]))
            elif keyword == 'bbconn':
                con_info = lbuffer[2:]
            lbuffer = string.split(f.readline())
    ### read body
    if ftype == 'xyz':
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                read_body(f,mol.natoms,frags=True)
    elif ftype == 'topo':
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers,\
                mol.pconn = read_body(f,mol.natoms,frags=True, topo = True)
        pass
    else:
        ftype = 'xyz':
        logging.warning('Unknown mfpx file type specified. Using xyz as default')
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                read_body(f,mol.natoms,frags=False)
    ### pass bb info
    try:
        line = f.readline().split()
        if line != [] and line[0][:5] == 'angle':
            self.angleterm = line
    if 'con_info' in locals():
        mol.dummies = []
        mol.dummy_neighbors=[]
        mol.connectors=[]
        mol.connectors_type=[]
        contype = 0
        for c in con_info:
            if c == "/":
                contype_count += 1
            else:
                ss = c.split('*') # ss[0] is the dummy neighbors, ss[1] is the connector atom
                if len(ss) != 2: raise IOError('This is not a proper BB file, convert with script before!')
                stt = ss[0].split(',')
                mol.connectors.append(int(ss[1])-1)
                mol.connectors_type.append(contype_count)
                if string.lower(self.elems[int(ss[1])-1]) == 'x':
                    mol.dummies.append(int(ss[1])-1) # simplest case only with two atoms being the connecting atoms
                    #self.natoms += 1
                mol.dummy_neighbors.append((np.array(map(int,stt)) -1).tolist())
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
