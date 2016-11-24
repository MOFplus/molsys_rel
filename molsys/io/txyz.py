import numpy
import string
from molsys import *

def read(mol, fname, topo = False):
    """
    Routine, which reads an txyz file
    :Parameters:
        -fname  (str) : name of the txyz file
        -mol    (obj) : instance of a molclass
        -topo   (bool): flag for reading topo information
    """
    f = open(fname, "r")
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
        mol.center_point = lbuffer[1]
        con_info = lbuffer[2:]
        pass_connstring(mol,con_info, new = False)
    if topo == False:
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers =\
                read_body(f,mol.natoms,frags=False)
    else:
        mol.elems, mol.xyz, mol.atypes, mol.conn, mol.fragtypes, mol.fragnumbers,\
                mol.pconn = read_body(f,mol.natoms,frags=False, topo = True)
    ### this has to go at some point
    if 'con_info' in locals():
        if self.center_point == "special":
            line = string.split(f.readline)
            self.special_center_point = numpy.array(map(string.atof,line[0:3]),"d")
        try:
            line = f.readline().split()
            if line != [] and line[0][:5] == 'angle':
                mol.angleterm = line
        except:
            pass
    return 

def pass_connstring(mol, con_info, new = True):
    """
    Routines which passes the con_info string of a txyz or an mfpx file
    :Parameters:
        - mol      (obj) : instance of a molclass
        - con_info (str) : string holding the connectors info
        - new      (bool): bool to switch between old and new type of con_info string
    """
    mol.dummies = []
    mol.dummy_neighbors = []
    mol.connectors = []
    mol.connectors_type = []
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
                    mol.dummies.append(int(ss[1])-1) # simplest case only with two atoms being the connecting atoms
                    #self.natoms += 1
                mol.dummy_neighbors.append((numpy.array(map(int,stt)) -1).tolist())
            else:
                mol.connectors.append(string.atoi(c)-1)
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
        return elems, numpy.array(xyz), atypes, conn, fragtypes, fragnumbers, pconn
    else:
        return elems, numpy.array(xyz), atypes, conn, fragtypes, fragnumbers


def write_body(f, mol, frags=True, topo=False):
    """
    Routine, which writes the body of a txyz or a mfpx file
    :Parameters:
        -f      (obj)  : fileobject
        -mol    (obj)  : instance of molsys object
        -frags  (bool) : flag to specify if fragment info should be in body or not
        -topo   (bool) : flag to specigy if pconn info should be in body or not
    """
    elems       = mol.elems
    atypes      = mol.atypes
    xyz         = mol.xyz
    cnct        = mol.conn
    natoms      = mol.natoms
    if frags == True:
        fragtypes   = mol.fragtypes
        fragnumbers = mol.fragnumbers
    if topo: pconn = mol.pconn
    for i in xrange(mol.natoms):
        line = ("%3d %-3s" + 3*"%12.6f" + " %10s") % \
            tuple([i+1]+[elems[i]]+ xyz[i].tolist() + [atypes[i]])
        if frags == True: line += ("%10s %3d") % tuple([fragtypes[i]]+[fragnumbers[i]])
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
    return
    

def write(mol, fname, topo = False, frags = False):
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
        f.write("%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" % tuple([mol.natoms]+cellparams))
    else:
        f.write("%5d \n" % mol.natoms)
    write_body(f, mol, topo=topo, frags = frags)
    f.close()
    return
