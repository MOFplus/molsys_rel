"""a simple wrapper to run the javascript systrekey from python
"""

import os
import subprocess
import string
import molsys


# get my path 
molsys_path = os.path.dirname(molsys.__file__)
systre_path = os.path.dirname(os.path.dirname(molsys_path)) + "/systreKey"

def run_systrekey(edges, labels):
    """run javascript systrekey
    
    Arguments:
        edges {list of lists of 2 ints} - list of edges
        labels {lost of lists of 3 ints} - list of edge labels
    """
    assert len(edges) == len(labels)
    lqg_string = ""
    for e,l in zip(edges, labels):
        lqg_string += "%d %d" % tuple(e)
        lqg_string += " %d %d %d|" % tuple(l)
    lqg_string = lqg_string[:-1]
    print lqg_string # DEBUG

    key = subprocess.check_output(args=["node", molsys_path+"/util/run_systrekey.js", lqg_string, systre_path])
    key = key[:-1] # remove newline
    return key


if __name__=="__main__":
    edges = [[1,1], [1,1], [1,1]]
    labels = [[1,0,0], [0,1,0], [0,0,1]]
    print run_systrekey(edges, labels)








