"""a simple wrapper to run the javascript systrekey from python

    TBI: i guess it would be much better to use json to communicate with javascript
         -> need to figure out how to do this.

    added a class for systrekey handling

    author: R. Schmid
"""

import os
import subprocess
import string
import molsys


# get my path 
molsys_path = os.path.dirname(molsys.__file__)
systre_path = os.path.dirname(os.path.dirname(molsys_path)) + "/systreKey"

# path to the curent arc file -- should be called RCSR.arc by default (softlink to current version)
arc_file = os.path.dirname(__file__)+"/RCSR.arc"

# we read the arc file once when the first instance is made 
arc_read = False

db_key2name = {}
db_name2key = {}

# shift_map  = {"-1" : "-", "0": "0", "1": "+"}
# shift_rmap = {"-" : "-1", "0": "0", "+": "1"}

class systre_db:

    def __init__(self):
        if not arc_read:
            self.read_arc()
        return

    def read_arc(self):
        global db_key2name, db_name2key, arc_read
        if os.path.isfile(arc_file):
            f = open(arc_file, 'r')
            for line in f:
                sline = line.split()
                if len(sline)>0:
                    if sline[0] == 'key':
                        key = sline[1:]
                    elif sline[0] == 'id':
                        name = sline[1]
                    elif sline[0] == "end":
                        # end of record .. store in directories only 3dim nets
                        if key[0] == "3":
                            key = string.join(key)
                            db_key2name[key] = name
                            db_name2key[name] = key
                    else:
                        pass
            f.close()
        else:
            print("""
            WARNING: the file RCSR.arc is not available.
            Please download the current arc file from RCSR.net/systre into
            the molsys/util directory and softlink it to RCSR.arc
            """)
        arc_read=True
        return

    def get_key(self, name):
        return db_name2key[name]

    def get_name(self, key):
        return db_key2name[key]

    # def compress_key(self, skey):
    #     """compress a systrekey
        
    #     Args:           
    #         self (string or list): original systre key
    #     """
    #     if type(skey) == type(""):
    #         lsk = skey.split()
    #     else:
    #         lsk = skey
    #     lsk.reverse()
    #     dim = lsk.pop() # dimension
    #     assert dim == "3"
    #     csk = ""
    #     while len(lsk)>0:
    #         i = lsk.pop()
    #         j = lsk.pop()
    #         csk += "%s:%s:" % (i,j)
    #         for i in range(3):
    #             csk += shift_map[lsk.pop()]
    #         if len(lsk)>0:
    #             csk+= ":"
    #     return csk


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
    # print lqg_string # DEBUG

    key = subprocess.check_output(args=["node", molsys_path+"/util/run_systrekey.js", lqg_string, systre_path])
    key = key[:-1] # remove newline
    return key


if __name__=="__main__":
    edges = [[1,1], [1,1], [1,1]]
    labels = [[1,0,0], [0,1,0], [0,0,1]]
    print run_systrekey(edges, labels)








