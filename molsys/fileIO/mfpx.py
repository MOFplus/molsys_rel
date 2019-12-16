import numpy
from . import txyz
import logging
from molsys.util.misc import argsorted


logger = logging.getLogger("molsys.io")

def read(mol, f):
    """
    Read mfpx file
    :Parameters:
        -f   (obj): mfpx file object or mfpx-like readable object
        -mol (obj): instance of a molclass
    """
    ### read check ###
    try:
        f.readline ### do nothing
    except AttributeError:
        raise IOError("%s is not readable" % f)
    ### read func ###
    ftype = 'xyz'
    lbuffer = f.readline().split()
    stop = False
    bbinfo = {}
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
                bbinfo['center_type'] = lbuffer[2]
                if lbuffer[2] == 'special':
                    bbinfo['center_xyz'] = numpy.array([float(i) for i in lbuffer[3:6]])
            elif keyword == 'bbconn':
                mol.is_bb = True
                con_info = lbuffer[2:]
            elif keyword == 'bbatypes':
                mol.is_bb = True
                bbinfo['connector_atypes'] = lbuffer[2:]
            elif keyword == 'orient':
                orient = [int(i) for i in lbuffer[2:]]
                mol.orientation = orient
            lbuffer = f.readline().split()
    if mol.is_bb == True: ftype = 'bb'
    ### read body
    if ftype == 'xyz':
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
            txyz.read_body(f,mol.natoms,frags=True)
    elif ftype == 'topo':
        # topo file so set the corresponding flags
        mol.is_topo =   True
        mol.use_pconn = True
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.pconn, mol.pimages =\
            txyz.read_body(f,mol.natoms,frags=True, topo = True)
    elif ftype == 'cromo':
        mol.is_topo =   True
        mol.is_cromo =   True
        mol.use_pconn = True
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.pconn, mol.pimages, mol.oconn =\
            txyz.read_body(f,mol.natoms,frags=True, topo = True, cromo = True)
    elif ftype == 'bb': # 
        connector, connector_atoms, connector_types = parse_connstring(mol,con_info)
        #bbinfo['connector'] = connector  ## that actually has to be an argument as the setup currently is coded.
        bbinfo['connector_atoms'] = connector_atoms
        bbinfo['connector_types'] = connector_types
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
            txyz.read_body(f,mol.natoms,frags=True, topo = False, cromo = False)
        mol.addon('bb')
        mol.bb.setup(connector,**bbinfo)
    else:
        ftype = 'xyz'
        logger.warning('Unknown mfpx file type specified. Using xyz as default')
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                txyz.read_body(f,mol.natoms,frags=False)
    mol.set_ctab_from_conn(pconn_flag=mol.use_pconn)
    mol.set_etab_from_tabs()
    if ftype == 'cromo':
        mol.set_otab_from_oconn()
    ### pass bb info
    try:
        line = f.readline().split()
        if line != [] and line[0][:5] == 'angle':
            mol.angleterm = line
    except:
        pass
    return

def write(mol, f, fullcell = True):
    """
    Routine, which writes an mfpx file
    :Parameters:
        -mol   (obj) : instance of a molsys class
        -f (obj) : file object or writable object
        -fullcell  (bool): flag to specify if complete cellvectors arre written
    """
    ### write check ###
    try:
        f.write ### do nothing
    except AttributeError:
        raise IOError("%s is not writable" % f)
    ### write func ###
    if len(mol.fragtypes) == 0:
        mol.set_nofrags()
    if mol.is_topo:
        ftype = 'topo'
        if mol.use_pconn == False:
            # this is a topo object but without valid pconn. for writing we need to generate it
            mol.add_pconn()
    elif mol.is_bb:
        ftype = "bb"
    else:
        ftype = 'xyz'
    f.write('# type %s\n' % ftype)
    if mol.bcond>0:
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
        ### BLOCK DEFINITION ###################################################
        ### bbcenter write ###
        if mol.center_point != 'special':
            f.write('# bbcenter %s\n' % mol.center_point)
        else:
            f.write('# bbcenter %s %12.6f %12.6f %12.6f\n' %
                    tuple([mol.center_point]+ mol.special_center_point.tolist()))
        ### bbconn write ###
        # sort by connectors' types
        _argsorted_type = argsorted(mol.connectors_type)
        connectors_type = [mol.connectors_type[_] for _ in _argsorted_type]
        connector_atoms = [sorted(mol.connector_atoms[_]) for _ in _argsorted_type]
        connectors      = [mol.connectors[_] for _ in _argsorted_type]
        # sort by connector atoms after connectors' types
        _sort_priority = [_*mol.natoms for _ in connectors_type]
        connector_atoms_priority = [
            [__+_sort_priority[_k] for __ in _] for _k,_ in enumerate(connector_atoms)
        ] # preparatory for nested sorting
        _argsorted_catoms = argsorted(connector_atoms_priority)
        connector_atoms = [connector_atoms[_] for _ in _argsorted_catoms]
        connectors      = [connectors[_]      for _ in _argsorted_catoms]
        # make bbcon line
        connstrings = ''
        ctype = 0
        for i,d in enumerate(mol.bb.connector_atoms):
            if mol.bb.connector_types[i] != ctype:
                ctype +=1
                connstrings += '/ ' # ctypes sep
            for j in d:
                connstrings = connstrings + str(j+1) +',' # connector_atoms sep
            connstrings = connstrings[0:-1] + '*' + str(connectors[i]+1)+' '
        # write bbconn line
        f.write('# bbconn %s\n' % connstrings)
        if mol.bb.connector_atypes is not None:
            # also write connector atypes
            temp_atypes = []
            for at in mol.bb.connector_atypes:
                if type(at) == type([]):
                    temp_atypes.append(",".join(at))
                else:
                    temp_atypes.append(at)
            f.write('# bbatypes %s\n' % " ".join(temp_atypes))
    if hasattr(mol, "orientation"):
        ### ORIENTATION DEFINITION #############################################
        o = len(mol.orientation) * "%3d" % tuple(mol.orientation)
        f.write('# orient '+o+"\n")
    f.write('%i\n' % mol.natoms)
    if ftype == 'xyz' or ftype == "bb":
        txyz.write_body(f,mol)
    else:
        txyz.write_body(f,mol,topo=True)
    return

def parse_connstring(mol, con_info): 
    """ 
    Routines which parses the con_info string of a txyz or an mfpx file 
    :Parameters: 
        - mol      (obj) : instance of a molclass 
        - con_info (str) : string holding the connectors info 
    """ 
    connector             = [] 
    connector_atoms       = [] 
    connector_types       = [] 
    contype_count = 0 
    for c in con_info: 
        if c == "/": 
            contype_count += 1 
        else: 
            ss = c.split('*') # ss[0] is the dummy neighbors, ss[1] is the connector atom 
            if len(ss) > 2: 
                raise IOError('This is not a proper BB file, convert with script before!') 
            elif len(ss) < 2: # neighbor == connector 
                ss *= 2 
            stt = ss[0].split(',') 
            connector.append(int(ss[1])-1) 
            connector_types.append(contype_count) 
            connector_atoms.append((numpy.array([int(i) for i in stt]) -1).tolist()) 
    return connector, connector_atoms, connector_types

