import numpy
import string
from . import txyz
import logging


logger = logging.getLogger("molsys.io")

def read(mol, f):
    """
    Routine, which reads an mfpx file
    :Parameters:
        -f   (obj): mfpx file object
        -mol (obj): instance of a molclass
    """
    ### read header ###
    ftype = 'xyz'
    lbuffer = f.readline().split()
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
                cellparams = [float(i) for i in lbuffer[2:8]]
                mol.set_cellparams(cellparams)
            elif keyword == 'cellvect':
                mol.periodic = True
                celllist = [float(i) for i in lbuffer[2:11]]
                cell = numpy.array(celllist)
                cell.shape = (3,3)
                mol.set_cell(cell)
            elif keyword == 'bbcenter':
                mol.is_bb = True
                mol.center_point = lbuffer[2]
                if mol.center_point == 'special':
                    mol.special_center_point = numpy.array([float(i) for i in lbuffer[3:6]])
            elif keyword == 'bbconn':
                mol.is_bb = True
                con_info = lbuffer[2:]
            elif keyword == 'orient':
                orient = [int(i) for i in lbuffer[2:]]
                mol.orientation = orient
            lbuffer = f.readline().split()
    ### read body
    if ftype == 'xyz':
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
            txyz.read_body(f,mol.natoms,frags=True)
    elif ftype == 'topo':
        if mol.__class__.__name__ != 'topo':
            logger.warning('Topology information is read to a regular mol object')
#        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers,\
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.pconn, mol.pimages =\
            txyz.read_body(f,mol.natoms,frags=True, topo = True)
    else:
        ftype = 'xyz'
        logger.warning('Unknown mfpx file type specified. Using xyz as default')
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                txyz.read_body(f,mol.natoms,frags=False)
    mol.set_ctab_from_conn()
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

def write(mol, fname, fullcell = True):
    """
    Routine, which writes an mfpx file
    :Parameters:
        -mol   (obj) : instance of a molsys class
        -fname (str) : name of the mfpx file
        -topo  (bool): flag to specify if pconn should be in mfpx file or not
    """
    if mol.fragtypes == []: mol.set_nofrags()
    f = open(fname, 'w')
    if mol.__class__.__name__ == 'topo':
        ftype = 'topo'
    else:
        ftype = 'xyz'
    f.write('# type %s\n' % ftype)
    if type(mol.cellparams) != type(None):
        if fullcell:
#            elif keyword == 'cellvect':
#                mol.periodic = True
#                celllist = [float(i) for i in lbuffer[2:11]]
#                cell = numpy.array(celllist)
#                cell.shape = (3,3)
#                mol.set_cell(cell)
            f.write('# cellvect %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n' %\
                    tuple(mol.cell.ravel()))
        else:
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
    if hasattr(mol, "orientation"):
        o = len(mol.orientation) * "%3d" % tuple(mol.orientation)
        f.write('# orient '+o+"\n")
    f.write('%i\n' % mol.natoms)
    if ftype == 'xyz':
        txyz.write_body(f,mol)
    else:
        txyz.write_body(f,mol,topo=True)
    f.close()
    return

