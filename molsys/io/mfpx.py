import numpy
import string
import unit_cell
import txyz
import logging

def read(mol, fname):
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
                ftype = lbuffer[2]
            elif keyword == 'cell':
                mol.cellparams = map(string.atof,lbuffer[2:8])
                mol.cell = unit_cell.vectors_from_abc(mol.cellparams)
            elif keyword == 'cellvect':
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
                txyz.read_body(f,mol.natoms,frags=True)
    elif ftype == 'topo':
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers,\
                mol.pconn = txyz.read_body(f,mol.natoms,frags=True, topo = True)
        pass
    else:
        ftype = 'xyz'
        logging.warning('Unknown mfpx file type specified. Using xyz as default')
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                txyz.read_body(f,mol.natoms,frags=False)
    ### pass bb info
    try:
        line = f.readline().split()
        if line != [] and line[0][:5] == 'angle':
            self.angleterm = line
    except:
        pass
    if 'con_info' in locals():
        mol.dummies = []
        mol.dummy_neighbors=[]
        mol.connectors=[]
        mol.connectors_type=[]
        contype_count = 0
        for c in con_info:
            if c == "/":
                contype_count += 1
            else:
                ss = c.split('*') # ss[0] is the dummy neighbors, ss[1] is the connector atom
                if len(ss) != 2: raise IOError('This is not a proper BB file, convert with script before!')
                stt = ss[0].split(',')
                mol.connectors.append(int(ss[1])-1)
                mol.connectors_type.append(contype_count)
                if string.lower(mol.elems[int(ss[1])-1]) == 'x':
                    mol.dummies.append(int(ss[1])-1) # simplest case only with two atoms being the connecting atoms
                    #self.natoms += 1
                mol.dummy_neighbors.append((numpy.array(map(int,stt)) -1).tolist())
    return

def write(mol, fname, topo = False, bb = False):
    """
    Routine, which writes an mfpx file
    :Parameters:
        -mol   (obj) : instance of a molsys class
        -fname (str) : name of the mfpx file
        -topo  (bool): flag to specify if pconn should be in mfpx file or not
        -bb    (bool): flag to specify if bb info should be in mfpx file or not
    """
    f = open(fname, 'w')
    if topo:
        ftype = 'topo'
    else:
        ftype = 'xyz'
    f.write('# type %s\n' % ftype)
    if type(mol.cellparams) != type(None):
        f.write('# cell %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n' %\
                tuple(mol.cellparams))
    if bb:
        if mol.centerpoint != 'special':
            f.write('# bbcenter %s\n' % mol.centerpoint)
        else:
            f.write('# bbcenter %s %12.6f %12.6f %12.6f\n' % 
                    tuple([mol.centerpoint]+ mol.special_center_point.tolist()))
        connstrings = ''
        for i,d in enumerate(mol.dummy_neighbors):
            for j in d:
                connstrings = connstrings + str(j+1) +','
            connstrings = connstrings[0:-1] + '*' + str(mol.connectors[i]+1)+' '
        f.write('# bbconn %s\n' % connstrings)
    f.write('%i\n' % mol.natoms)
    if ftype == 'xyz':
        txyz.write_body(f,mol)
    else:
        txyz.write_body(f,mol,topo=True)
    f.close()
    return

