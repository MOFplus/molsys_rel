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
        -f    (obj) : txyz file object
        -mol  (obj) : instance of a molclass
        -topo (bool): flag for reading topo information
    """
    lbuffer = f.readline().split()
    mol.natoms = int(lbuffer[0])
    if len(lbuffer) >1 and lbuffer[1] == "molden": lbuffer = [lbuffer[0]]
    if len(lbuffer) > 1 and lbuffer[1] not in ['special','coc','com']:
        boundarycond = 3
        if lbuffer[1] == "#":
            celllist = [float(i) for i in lbuffer[2:11]]
            cell = numpy.array(celllist)
            cell.shape = (3,3)
            mol.set_cell(cell)
        else:
            cellparams = [float(i) for i in lbuffer[1:7]]
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
            line = f.readline().split()
            mol.special_center_point = numpy.array([float(i) for i in line[0:3]],"d")
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
    mol.connector_dummies     = []
    mol.connector_atoms       = []
    mol.connectors            = []
    mol.connectors_type       = []
    mol.connectors_group      = []
    mol.connectors_complexity = []
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
                if mol.elems[int(ss[1])-1].lower() == 'x':
                    mol.connector_dummies.append(int(ss[1])-1) # simplest case only with two atoms being the connecting atoms
                    #self.natoms += 1
                mol.connector_atoms.append((numpy.array([int(i) for i in stt]) -1).tolist())
            else:
                # in the old format only 1:1 connections exists: dummy_neighbors  are equal to connector atoms
                mol.connectors.append(int(c)-1)
                mol.connector_atoms = [[c] for c in mol.connectors]
                mol.connectors_type.append(contype_count)
    return

def parse_connstring_new(mol, con_info, **kwargs):
    """
    Routines which parses the con_info string of a txyz or an mfpx file
    :Parameters:
        - mol      (obj) : instance of a molclass
        - con_info (str) : string holding the connectors info
    """
    mol.connector_atoms       = []
    mol.connector_dummies     = []
    mol.connectors            = []
    mol.connectors_complexity = []
    mol.connectors_group      = []
    mol.connectors_type       = []
    contype_count = 0
    for icon, con in enumerate(con_info):
        if con == "/":
            contype_count += 1
            continue
        ### PARSE ####
        logger.debug(con)
        line = con[:]
        markcount = line.count("?")
        starcount = line.count("*")
        if markcount == 0:
            complexity = ""
        elif markcount == 1:
            line, complexity = line.split('?')
        else:
            logger.error("More than one question mark in con_info group")
            raise ValueError
        if starcount == 0:
            connectors = ""
        elif starcount == 1:
            line, connectors = line.split('*')
        else:
            logger.error("More than one asterisk in con_info group")
            raise ValueError
        atoms = line.split(",")
        line = con #reset, debugging purpose
        logger.debug("a:%s c:%s t:%s *:%s ?:%s" % 
            (atoms, connectors, complexity, starcount, markcount))
        complexity = complexity.split()
        if complexity != []:
            complexity = complexity[0].split(",")
        connectors = connectors.split()
        if connectors != []:
            connectors = connectors[0].split(",")
        logger.debug("a:%s c:%s t:%s *:%s ?:%s" % 
            (atoms, connectors, complexity, starcount, markcount))
        logger.debug("la:%s lc:%s lt:%s *:%s ?:%s" % 
            (len(atoms), len(connectors), len(complexity), starcount, markcount))
        ### HANDLE ###
        if len(atoms) == 0:
            logger.error("No atoms in con_info group")
            raise ValueError
        if len(connectors) == 0:
            connectors = atoms[:]
        if len(complexity) == 0:
            complexity = [1]*len(atoms)
        elif len(complexity) == 1:
            complexity = complexity*len(atoms)
        if len(complexity) != len(atoms):
            logger.error("Complexity can only be: implicit OR one for all OR assigned per each")
            raise ValueError
        atoms = [int(a) - 1 for a in atoms] #python indexing
        connectors = [int(c) - 1 for c in connectors] #python indexing
        complexity = map(int,complexity)
        logger.debug("a:%s c:%s t:%s *:%s ?:%s" % 
            (atoms, connectors, complexity, starcount, markcount))
        logger.debug("la:%s lc:%s lt:%s *:%s ?:%s" % 
            (len(atoms), len(connectors), len(complexity), starcount, markcount))
        ### SET ######
        mol.connector_atoms.append(atoms) #((numpy.array(map(int,stt)) -1).tolist())
        for a in atoms:
            if mol.elems[a].lower() == "x":
                mol.connector_dummies.append(a) # simplest case only with two atoms being the connecting atoms
        mol.connectors.append(connectors) #(int(ss[1])-1)
        mol.connectors_complexity.append(complexity)
        mol.connectors_group.append(icon - contype_count)
        mol.connectors_type.append(contype_count)
    mol.connectors = numpy.array(mol.connectors)
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
    pimages     = []
    if topo: frags=False
    for i in range(natoms):
        lbuffer = f.readline().split()
        xyz.append([float(i) for i in lbuffer[2:5]])
        elems.append(lbuffer[1].lower())
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
            conn.append((numpy.array([int(i) for i in lbuffer[6+offset:]])-1).tolist())
        else:
            txt = lbuffer[6+offset:]
            a = [[int(j) for j in i.split('/')] for i in txt]
            c,pc,pim = [i[0]-1 for i in a], [images[i[1]] for i in a], [i[1] for i in a]
            conn.append(c)
            pconn.append(pc)
            pimages.append(pim)
    if topo:
#        return elems, numpy.array(xyz), atypes, conn, fragtypes, fragnumbers, pconn
        return elems, numpy.array(xyz), atypes, conn, pconn, pimages
    else:
        return elems, numpy.array(xyz), atypes, conn, fragtypes, fragnumbers


def write_body(f, mol, frags=True, topo=False, pbc=True, moldenr=False):
    """
    Routine, which writes the body of a txyz or a mfpx file
    :Parameters:
        -f      (obj)  : fileobject
        -mol    (obj)  : instance of molsys object
        -frags  (bool) : flag to specify if fragment info should be in body or not
        -topo   (bool) : flag to specigy if pconn info should be in body or not
        -pbc    (bool) : if False, removes connectivity out of the box (meant for visualization)
        -moldenrc (bool) : atypes compatible with molden
    """
    if topo: frags = False   #from now on this is convention!
    if topo: pconn = mol.pconn
    if frags == True:
        fragtypes   = mol.fragtypes
        fragnumbers = mol.fragnumbers
    else:
        fragtypes = None
        fragnumbers = None
    elems       = mol.elems
    xyz         = mol.xyz
    cnct        = mol.conn
    natoms      = mol.natoms
    atypes      = mol.atypes
    if pbc == False:
        cellcond = .5*numpy.array(mol.cellparams[:3])
        cnct = [ [i for i in c if (abs(xyz[i]-xyz[e]) < cellcond).all()] for e,c in enumerate(cnct)]
    if moldenr:
        if fragtypes is None: fragtypes = [None for i in mol.atypes]
        if fragnumbers is None: fragnumbers = [None for i in mol.atypes]
        from collections import Counter
        oldatypes = list(zip(fragnumbers, fragtypes, atypes))
        unique_oldatypes = list(Counter(oldatypes).keys())
        unique_oldatypes.sort()
        old2newatypes = {e:i for i,e in enumerate(unique_oldatypes)}
        new2oldatypes = {i:e for i,e in enumerate(unique_oldatypes)}
        newatypes = [old2newatypes[i] for i in oldatypes]
        atypes = newatypes
        frags = False ### encoded in one column only
    for i in range(mol.natoms):
        try:
            line = ("%3d %-3s" + 3*"%12.6f" + "   %-24s") % \
            tuple([i+1]+[elems[i]]+ xyz[i].tolist() + [atypes[i]])
        except IndexError:
            import pdb; pdb.set_trace()
        if frags == True: line += ("%-16s %5d ") % tuple([fragtypes[i]]+[fragnumbers[i]])
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
                        line +="%8d/%1d " % (cc,pp)
                    else:
                        line += "%7d/%2d " % (cc,pp)
            else:
                line += (len(conn)*"%7d ") % tuple(conn)
        f.write("%s \n" % line)
    if moldenr:
        satoms = len(str(mol.natoms))
        f.write("### type fragment_type fragment_number atom_type\n{\n")
        for key,val in new2oldatypes.items():
            f.write("    %%%ds: %%s,\n" % (satoms,) % (key,val))
        f.write("}\n")
    return


def write(mol, fname, topo = False, frags = False, pbc=True, moldenr=False):
    """
    Routine, which writes an txyz file
    :Parameters:
        -fname  (str) : name of the txyz file
        -mol    (obj) : instance of a molclass
        -topo   (bool): flag top specify if pconn should be in txyz file or not
    """
    cellparams = mol.cellparams
    f = open(fname, 'w')
    if cellparams is not None:
        f.write("%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" % tuple([mol.natoms]+list(cellparams)))
    else:
        f.write("%5d \n" % mol.natoms)
    write_body(f, mol, topo = topo, frags = frags, pbc = pbc, moldenr = moldenr)
    f.close()
    return
