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

# check presence of node in order to run javascript code
node_avail = (subprocess.call(args=["which", "node"], stdout=subprocess.DEVNULL) == 0)


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
                            key = " ".join(key)
                            if key in db_key2name:
                                print ("WARNING the following systrekey is already registered for %s" % db_key2name[key])
                                print (key)
                            db_key2name[key] = name
                            db_name2key[name] = key
                    else:
                        pass
            f.close()
        else:
            return
        arc_read=True
        return

    def report_missing_arc(self):
        print("""
            WARNING: the file RCSR.arc is not available.
            Please download the current arc file from RCSR.net/systre into
            the molsys/util directory and softlink it to RCSR.arc
            """)
        return

    def get_key(self, name):
        if arc_read:
            return db_name2key[name]
        else:
            self.report_missing_arc()
            return "None"

    def get_name(self, key):
        if arc_read:
            return db_key2name[key]
        else:
            self.report_missing_arc()
            return "None"

def run_systrekey(edges, labels):
    import json
    import copy
    """run javascript systrekey
    
    revised version using json strings to pass data to javascript

    Arguments:
        edges {list of lists of 2 ints} - list of edges
        labels {list of lists of 3 ints} - list of edge labels
    """
    if not node_avail:
        # return nothing to avoid an error
        print("""
        WARNING: systrekey was called but node is not installed to run javascript
        """
        )
        return "None", []
    assert len(edges) == len(labels)
    lqg = []
    for e,l in zip(edges, labels):
        el = []
        el.append(int(e[0])+1)
        el.append(int(e[1])+1)
        el.append(l)
        lqg.append(el)
    json_lqg = json.dumps(lqg)

    try:
        json_result = subprocess.check_output(args=["node", molsys_path+"/util/run_systrekey.js", json_lqg, systre_path], stderr=subprocess.STDOUT).decode()
    except subprocess.CalledProcessError as err:
        raw_err = err.stdout.decode().split("\n")
        for el in raw_err:
            if el[:6] == "Error:":
                err = el
                break
        print ("systrekey says -> %s" % err)
        return err, []

    result = json.loads(json_result)
    key = result["key"]
    mapping = result["mapping"]
    return key, mapping


if __name__=="__main__":
    edges = [[1,1], [1,1], [1,1]]
    labels = [[1,0,0], [0,1,0], [0,0,1]]
    print (run_systrekey(edges, labels))








