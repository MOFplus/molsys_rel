import numpy
import string
import txyz
import logging


logger = logging.getLogger("molsys.io")

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
                cellparams = map(string.atof,lbuffer[2:8])
                mol.set_cellparams(cellparams)
            elif keyword == 'cellvect':
                mol.periodic = True
                celllist = map(string.atof,lbuffer[2:11])
                cell = numpy.array(celllist)
                cell.shape = (3,3)
                mol.set_cell(cell)
            elif keyword == 'bbcenter':
                mol.is_bb = True
                mol.center_point = lbuffer[2]
                if mol.center_point == 'special':
                    mol.special_center_point = np.array(map(float,lbuffer[3:6]))
            elif keyword == 'bbconn':
                mol.is_bb = True
                con_info = lbuffer[2:]
            lbuffer = string.split(f.readline())
    ### read body
    if ftype == 'xyz':
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                txyz.read_body(f,mol.natoms,frags=True)
    elif ftype == 'topo':
        if mol.__class__.__name__ != 'topo':
            logger.warning('Topology information is read to a regular mol object')
#        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers,\
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.pconn = txyz.read_body(f,mol.natoms,frags=True, topo = True)
        pass
    else:
        ftype = 'xyz'
        logger.warning('Unknown mfpx file type specified. Using xyz as default')
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                txyz.read_body(f,mol.natoms,frags=False)
    ### pass bb info
    try:
        line = f.readline().split()
        if line != [] and line[0][:5] == 'angle':
            mol.angleterm = line
    except:
        pass
    if 'con_info' in locals():
        txyz.parse_connstring(mol,con_info)
    return

def write(mol, fname):
    """
    Routine, which writes an mfpx file
    :Parameters:
        -mol   (obj) : instance of a molsys class
        -fname (str) : name of the mfpx file
        -topo  (bool): flag to specify if pconn should be in mfpx file or not
    """
    f = open(fname, 'w')
    if mol.__class__.__name__ == 'topo':
        ftype = 'topo'
    else:
        ftype = 'xyz'
    f.write('# type %s\n' % ftype)
    if type(mol.cellparams) != type(None):
        f.write('# cell %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n' %\
                tuple(mol.cellparams))
    if mol.is_bb:
        if mol.center_point != 'special':
            f.write('# bbcenter %s\n' % mol.center_point)
        else:
            f.write('# bbcenter %s %12.6f %12.6f %12.6f\n' %
                    tuple([mol.center_point]+ mol.special_center_point.tolist()))
        connstrings = ''
        ctype = 0
        for i,d in enumerate(mol.connector_atoms):
            if mol.connectors_type[i] != ctype:
                ctype +=1
                connstrings += '/ '
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

