"""systre

This helper module contains functions to convert a systrekey to an embedding by calling systre


"""
import subprocess
import tempfile
import molsys


def run_systre(key):
    """run systreCmd via jython

    probably this could all be run using jython directly, but i will try to keep these things seperate
    so jypthon calls will be done as a subprocess and we analyze the output.
    
    Args:       
        key (string): systrekey to be converted
    """
    lsk = key.split()
    assert lsk[0] == "3"
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
    systre_result = systre_result.split("\n")
    l = 0 # line pointer
    # first get equivalences
    equiv = []
    line = systre_result[l].split()
    stop = False
    while not stop:
        if len(line)>0:
            if line[0] == "Equivalences":
                stop = True
        l += 1
        line = systre_result[l].split()
    nunique = int(line[0])-1
    for i in range(nunique):
        equiv.append(i)
    stop = False
    while not stop:
        equiv.append(int(line[2])-1)
        l += 1
        line = systre_result[l].split()
        if len(line) == 0:
            stop = True
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
    

    return











if __name__=="__main__":
    key = """3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 2 6 0 0 0 2 7 0 0 0 3 8 0 0 0 3 9 0 0 0 4 8 0 1 0 
               4 9 1 0 0 5 6 -1 1 1 5 7 0 0 1 6 10 0 0 0 6 11 0 0 0 7 12 0 0 0 7 13 0 0 0 8 10 -1 0 
               1 8 11 0 0 0 9 12 0 0 0 9 13 -1 0 1 10 14 0 0 0 11 14 -1 0 0 12 14 -1 0 0 13 14 -1 1 0"""
    run_systre(key)

