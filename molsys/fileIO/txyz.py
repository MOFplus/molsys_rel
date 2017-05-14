import numpy
import string
#from molsys import *
import molsys.util.images as images
import logging

logger = logging.getLogger("molsys.io")

def read(mol, f, topo = False):
    """
    Routine, which reads an txyz file
    :Parameters:
        -fname  (str) : name of the txyz file
        -mol    (obj) : instance of a molclass
        -topo   (bool): flag for reading topo information
    """
    lbuffer = string.split(f.readline())
    mol.natoms = string.atoi(lbuffer[0])
    if len(lbuffer) > 1 and lbuffer[1] not in ['special','coc','com']:
        boundarycond = 3
        if lbuffer[1] == "#":
            celllist = map(string.atof,lbuffer[2:11])
            cell = numpy.array(celllist)
            cell.shape = (3,3)
            mol.set_cell(cell)
        else:
            cellparams = map(string.atof, lbuffer[1:7])
            mol.set_cellparams(cellparams)
        if ((cellparams[3]==90.0) and (cellparams[4]==90.0) and (cellparams[5]==90.0)):
            boundarycond=2
            if ((cellparams[0]==cellparams[1])and(cellparams[1]==cellparams[2])and\
                (cellparams[0]==cellparams[2])):
                    boundarycond=1
    elif len(lbuffer) > 1:
        mol.is_bb=True
        mol.center_point = lbuffer[1]
        con_info = lbuffer[2:]
        parse_connstring(mol,con_info, new = False)
    if topo == False:
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                read_body(f,mol.natoms,frags=False)
    else:
        mol.elems, mol.xyz, mol.atypes, mol.conn,\
                mol.pconn = read_body(f,mol.natoms,frags=False, topo = True)
    ### this has to go at some point
    if 'con_info' in locals():
        if mol.center_point == "special":
            line = string.split(f.readline())
            mol.special_center_point = numpy.array(map(string.atof,line[0:3]),"d")
        try:
            line = f.readline().split()
            if line != [] and line[0][:5] == 'angle':
                mol.angleterm = line
        except:
            pass
    
    mol.set_nofrags()
    return

def parse_connstring(mol, con_info, new = True):
    """
    Routines which parses the con_info string of a txyz or an mfpx file
    :Parameters:
        - mol      (obj) : instance of a molclass
        - con_info (str) : string holding the connectors info
        - new      (bool): bool to switch between old and new type of con_info string
    """
    mol.connector_dummies = []
    mol.connector_atoms   = []
    mol.connectors        = []
    mol.connectors_type   = []
    contype_count = 0
    for c in con_info:
        if c == "/":
            contype_count += 1
        else:
            if new:
                ss = c.split('*') # ss[0] is the dummy neighbors, ss[1] is the connector atom
                if len(ss) != 2: raise IOError('This is not a proper BB file, convert with script before!')
                stt = ss[0].split(',')
                mol.connectors.append(int(ss[1])-1)
                mol.connectors_type.append(contype_count)
                if string.lower(mol.elems[int(ss[1])-1]) == 'x':
                    mol.connector_dummies.append(int(ss[1])-1) # simplest case only with two atoms being the connecting atoms
                    #self.natoms += 1
                mol.connector_atoms.append((numpy.array(map(int,stt)) -1).tolist())
            else:
                # in the old format only 1:1 connections exists: dummy_neighbors  are equal to connector atoms
                mol.connectors.append(string.atoi(c)-1)
                mol.connector_atoms = [[c] for c in mol.connectors]
                mol.connectors_type.append(contype_count)
    return

def read_body(f, natoms, frags = True, topo = False):
    """
    Routine, which reads the body of a txyz or a mfpx file
    :Parameters:
        -f      (obj)  : fileobject
        -natoms (int)  : number of atoms in body
        -frags  (bool) : flag to specify if fragment info is in body or not
        -topo   (bool) : flag to specigy if pconn info is in body or not
    """
    elems       = []
    xyz         = []
    atypes      = []
    conn        = []
    fragtypes   = []
    fragnumbers = []
    pconn       = []
    if topo: frags=False
    for i in xrange(natoms):
        lbuffer = string.split(f.readline())
        xyz.append(map(string.atof, lbuffer[2:5]))
        elems.append(string.lower(lbuffer[1]))
        t = lbuffer[5]
        atypes.append(t)
        if frags == True:
            fragtypes.append(lbuffer[6])
            fragnumbers.append(int(lbuffer[7]))
            offset = 2
        else:
            fragtypes.append('0')
            fragnumbers.append(0)
            offset = 0
        if topo == False:
            conn.append((numpy.array(map(string.atoi, lbuffer[6+offset:]))-1).tolist())
        else:
            txt = lbuffer[6+offset:]
            a = [map(int,i.split('/')) for i in txt]
            c,pc = [i[0]-1 for i in a], [images[i[1]] for i in a]
            conn.append(c)
            pconn.append(pc)
    if topo:
#        return elems, numpy.array(xyz), atypes, conn, fragtypes, fragnumbers, pconn
        return elems, numpy.array(xyz), atypes, conn, pconn
    else:
        return elems, numpy.array(xyz), atypes, conn, fragtypes, fragnumbers


def write_body(f, mol, frags=True, topo=False, moldenr=False):
    """
    Routine, which writes the body of a txyz or a mfpx file
    :Parameters:
        -f      (obj)  : fileobject
        -mol    (obj)  : instance of molsys object
        -frags  (bool) : flag to specify if fragment info should be in body or not
        -topo   (bool) : flag to specigy if pconn info should be in body or not
    """
    if topo: frags = False   #from now on this is convention!
    if topo: pconn = mol.pconn
    if frags == True:
        fragtypes   = mol.fragtypes
        fragnumbers = mol.fragnumbers
    elems       = mol.elems
    xyz         = mol.xyz
    cnct        = mol.conn
    natoms      = mol.natoms
    if moldenr: ###TO BE DEBUGGED
        moltypes = {}
        moldentypes = []
        for i in xrange(mol.natoms):
            #if frags: item = '__'.join([mol.atypes[i], mol.fragtypes[i], str(mol.fragnumbers[i])]) 
            if frags: item = '__'.join([mol.atypes[i], mol.fragtypes[i]]) 
            else: item = mol.atypes[i]
            if not item in moltypes:
                moltypes[item]=(len(moltypes)+1, mol.elems[i])
            moldentypes.append( moltypes[item][0] )
        atypes = moldentypes
        frags=False ### => only one column all together
    else:
        atypes      = mol.atypes
    for i in xrange(mol.natoms):
        line = ("%3d %-3s" + 3*"%12.6f" + "   %-24s") % \
            tuple([i+1]+[elems[i]]+ xyz[i].tolist() + [atypes[i]])
        if frags == True: line += ("%-16s %5d") % tuple([fragtypes[i]]+[fragnumbers[i]])
        conn = (numpy.array(cnct[i])+1).tolist()
        if len(conn) != 0:
            if topo:
                pimg = []
                for pc in pconn[i]:
                    for ii,img in enumerate(images):
                        if all(img==pc):
                            pimg.append(ii)
                            break
                for cc,pp in zip(conn,pimg):
                    if pp < 10:
                        line +="%8d/%1d" % (cc,pp)
                    else:
                        line += "%7d/%2d" % (cc,pp)
            else:
                line += (len(conn)*"%7d") % tuple(conn)
        f.write("%s \n" % line)
    if moldenr:
        f.write("#atomtypes \n")
        for i in moltypes:
            f.write("%s %s %s \n" % (moltypes[i][0], i, moltypes[i][1]))
    return


def write(mol, fname, topo = False, frags = False, moldenr=False):
    """
    Routine, which writes an txyz file
    :Parameters:
        -fname  (str) : name of the txyz file
        -mol    (obj) : instance of a molclass
        -topo   (bool): flag top specify if pconn should be in txyz file or not
    """
    cellparams = mol.cellparams
    f = open(fname, 'w')
    if type(cellparams) != type(None):
        f.write("%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" % tuple([mol.natoms]+list(cellparams)))
    else:
        f.write("%5d \n" % mol.natoms)
    write_body(f, mol, topo=topo, frags = frags, moldenr=moldenr)
    f.close()
    return