"""systre

This helper module contains functions to convert a systrekey to an embedding by calling systre

REMARKS:
    - systre gives both vertices and edges in terms of fractional coordintes but these can be in the range [0.0, 1.0]
      which means vertices can have either 0.0 or 1.0 as a coordinate which is the same.
      In order to make the edge search working we use convention coords to be in the range [0.0, 1.0[
      (essentially all 1.0 are converted to 0.0)


"""
import numpy as np
import subprocess
import tempfile
import molsys
from molsys.util import unit_cell
from molsys.util import elems


def run_systre(key, debug=False):
    """run systreCmd via jython

    probably this could all be run using jython directly, but i will try to keep these things seperate
    so jypthon calls will be done as a subprocess and we analyze the output.
    
    Args:       
        key (string): systrekey to be converted
    """
    lsk = key.split()
    assert lsk[0] == "3"
    key = " ".join(lsk) # cononicalize key string for later comparison
    nedges = int((len(lsk)-1)/5)
    # now generate a cgd input for systre
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cgd", delete=True) as fcgd:
        fcgd.write("PERIODIC_GRAPH\nEDGES\n")
        for i in range(nedges):
            edge = lsk[i*5+1:i*5+6]
            fcgd.write("%s\n" % " ".join(edge))
        fcgd.write("END\n")
        fcgd.flush()
        try:
            systre_result = subprocess.check_output(args=["jython", "-u", "/home/rochus/code/systreCmd.py", "-u", fcgd.name], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as err:
            raw_err = err.stdout.decode().split("\n")
            print (raw_err)
    # now parse the systre output - we use the regular output but for a P1 embedding (-u)
    if debug:
        print (systre_result.decode())
    systre_result = systre_result.decode().split("\n")
    l = 0 # line pointer
    # first get equivalences
    equiv = []
    line = systre_result[l].split()
    stop = False
    while not stop:
        if len(line)>0:
            if line[1] == "nodes":
                stop = True
                break
        l += 1
        line = systre_result[l].split()
    nvert = int(line[0])
    stop = False
    has_equiv = True
    while not stop:
        if len(line)>0:
            if line[0] == "Equivalences":
                stop = True
            if line[0] == "Coordination":
                stop = True
                has_equiv = False
        l += 1
        line = systre_result[l].split()
    equiv = list(range(nvert))
    if has_equiv:
        stop = False
        while not stop:
            equiv[(int(line[0])-1)] = (int(line[2])-1)
            l += 1
            line = systre_result[l].split()
            if len(line) == 0:
                stop = True
    # equiv is not necessarily continuous from 0 to nunique
    symlabels = list(set(equiv))
    rev_equiv = [symlabels.index(e) for e in equiv]
    # next get space group
    stop = False
    while not stop:
        if len(line) == 5:
            if line[:3] == ["Ideal", "space", "group"]:
                spgroup = line[4][:-1]
                stop = True
        l += 1
        line = systre_result[l].split()
    # get cell params
    stop = False
    while not stop:
        if len(line) == 3:
            if line == ["Relaxed", "cell", "parameters:"]:
                stop = True
        l += 1
        line = systre_result[l].split()
    a = float(line[2][:-1])
    b = float(line[5][:-1])
    c = float(line[8][:-1])
    l += 1
    line = systre_result[l].split()
    alpha = float(line[2][:-1])
    beta  = float(line[5][:-1])
    gamma = float(line[8][:-1])
    # get vertices
    vertices = []
    stop = False
    while not stop:
        if len(line) == 2:
            if line == ["Relaxed", "positions:"]:
                stop = True
        l += 1
        line = systre_result[l].split()
    stop = False
    while not stop:
        if line [0] == "Node":
            vertices.append([float(x) for x in line[2:5]])
            l += 1
            line = systre_result[l].split()
        else:
            stop = True
    nvertices = len(vertices)
    assert line[0] == "Edges:"
    # get edges from to
    edge_from = []
    edge_to   = []
    l += 1
    line = systre_result[l].split()
    stop = False
    while not stop:
        if line [0] != "Edge":
            edge_from.append([float(x) for x in line[0:3]])
            edge_to.append([float(x) for x in line[4:7]])
            l += 1
            line = systre_result[l].split()
        else:
            stop = True
    nedges = len(edge_from)
    # now convert edges to a conn/pconn for our molsys topo object
    cell = unit_cell.vectors_from_abc((a, b, c, alpha, beta, gamma))
    fxyz = np.array(vertices)
    # make sure that all fxyz are in [0.0, 1.0[
    fxyz = np.where(fxyz==1.0, 0.0, fxyz)
    xyz = np.dot(fxyz, cell)
    m = molsys.mol.from_array(xyz)
    m.set_cell(cell)
    m.set_empty_conn()
    m.set_empty_pconn()
    # determine conn and pconn from the edges .. first we need to identify the vertices in edge_from
    v1 = np.array(edge_from)
    v2 = np.array(edge_to)
    # get image and base vectors
    v1_img = v1//1.0
    v1= v1%1.0
    v2_img = v2//1.0
    v2 = v2%1.0
    for e in range(v1.shape[0]):
        # get i vertex (from)
        d = fxyz-v1[e]
        d2 = (d*d).sum(axis=1)
        verti = np.argmin(d2)
        assert d2[verti] < 1.0e-4, "could not identify vertex with a short dist of %12.6f" % d2[verti]
        # get j vertex (to)
        d = fxyz-v2[e]
        d2 = (d*d).sum(axis=1)
        vertj = np.argmin(d2)
        assert d2[vertj] < 1.0e-4, "could not identify vertex with a short dist of %12.6f" % d2[vertj]
        # compute image offset
        img = v2_img[e]-v1_img[e]
        # now set conn and pconn (edges from systre are allready bidirectional)
        img = img.astype("int32")
        m.conn[verti].append(int(vertj))
        m.pconn[verti].append(img)
        #print ("bond %3d %3d img %s" % (verti, vertj, str(img)))
        # m.conn[vertj].append(verti)
        # m.pconn[vertj].append(img*-1)
    # DEBUG DEBUG check bonds 
    #for i in range(m.get_natoms()):
    #    for ij, j in enumerate(m.conn[i]):
    #        if i < j:
    #            ii = m.conn[j].index(i)
    #            print ("bond %3d %3d img %s .. reverse img %s" % (i, j, m.pconn[i][ij], m.pconn[j][ii]))
    # conn/pconn done
    el = []
    for i in range(m.get_natoms()):
        e = elems.topotypes[len(m.conn[i])]
        el.append(e)
    m.set_elems(el)
    m.is_topo = True
    m.use_pconn = True
    # now make sure that the conversion went ok by running systrekey with the built topo .. the key must match the input
    m.addon("topo")
    new_key = m.topo.systrekey
    if key != new_key:
        print ("someting went wrong here")
        print (key)
        print (new_key)
    # at this point we can be sure that the embedding went well and we have all the mappings
    # now the symmetry unique vertices have to be written as strings to atype
    key_mapping = m.topo.skey_mapping
    sym_mapping = []
    for v in key_mapping:
        sym_mapping.append(str(rev_equiv[v]))
    m.set_atypes(sym_mapping)
    m.topo._spgr = spgroup
    return m


def convert_all():
    """this function converts all systrekeys in the systreky db imported from systrekey into embeddings as long as they are not present, yet.
    """
    from molsys.util import systrekey
    import os
    for n in systrekey.db_name2key.keys():
        if not os.path.isfile(n+".mfpx"):
            print ("generating embedding for net %s " % n)
            m = run_systre(systrekey.db_name2key[n])
            m.write(n+".mfpx")
            print ("done!")
    return



# for debugging call this with a RCSR name
if __name__=="__main__":
    import sys
    from molsys.util import systrekey
    name = sys.argv[1]
    # thsi si the key of tbo ... just for testing
    # key = """3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 2 6 0 0 0 2 7 0 0 0 3 8 0 0 0 3 9 0 0 0 4 8 0 1 0 
    #            4 9 1 0 0 5 6 -1 1 1 5 7 0 0 1 6 10 0 0 0 6 11 0 0 0 7 12 0 0 0 7 13 0 0 0 8 10 -1 0 
    #            1 8 11 0 0 0 9 12 0 0 0 9 13 -1 0 1 10 14 0 0 0 11 14 -1 0 0 12 14 -1 0 0 13 14 -1 1 0"""
    assert name in systrekey.db_name2key
    key = systrekey.db_name2key[name]
    m = run_systre(key, debug=True)
    m.write(name+".mfpx")

