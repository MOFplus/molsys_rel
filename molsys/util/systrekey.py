"""a simple wrapper to run the javascript systrekey from python

    added a class for systrekey handling

    author: R. Schmid
"""

import os
import subprocess
import string
import molsys

"""
INSTALLATION Instructions:

For this to work you need to have the javascript systrekey code from Olaf Delgado-Friedrichs installed.

You need to have a working node.js environment working (do not use the distro version .. they are often too old) from https://nodejs.org

then clone the systreKey repo from here https://github.com/odf/systreKey to the same directory where your molsys repo is located (e.g. ~/sandbox)
because systre_path is determined from the root path of molsys.
For the installation follow the instructions from the systreKey repos README 

    git clone https://github.com/odf/systreKey.git
    cd systreKey
    npm install
    npm run build

Do not forget to run "npm run build" every time you pulled updates from Olaf's repo.

"""


# get my path 
molsys_path = os.path.dirname(molsys.__file__)
systre_path = os.path.dirname(os.path.dirname(molsys_path)) + "/systreKey"

# check presence of node in order to run javascript code
node_avail = (subprocess.call(args=["which", "node"], stdout=subprocess.DEVNULL) == 0)


class systre_db:

    def __init__(self, arc_file):
        self.arc_file = arc_file
        self.key2name = {}
        self.name2key = {}
        self.arc_read = False
        self.read_arc()
        return

    def read_arc(self):
        global db_key2name, db_name2key, arc_read
        if os.path.isfile(self.arc_file):
            f = open(self.arc_file, 'r')
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
                            if key in self.key2name:
                                print ("WARNING the following systrekey is already registered for %s" % self.key2name[key])
                                print (key)
                            self.key2name[key] = name
                            self.name2key[name] = key
                    else:
                        pass
            f.close()
            self.arc_read = True
        else:
            print("""
            WARNING: the file %s
            is not available. Please link it into the molsys/util directory.
            """ % self.arc_file)
        return

    def get_key(self, name):
        if self.arc_read:
            return self.name2key[name]
        else:
            return "None"

    def get_name(self, key):
        if self.arc_read:
            return self.key2name[key]
        else:
            return "None"

# path to the curent arc file -- should be called RCSR.arc by default (softlink to current version)
RCSR_arc_file = os.path.dirname(__file__)+"/RCSR.arc"
RCSR_db = systre_db(RCSR_arc_file)

db_key2name = RCSR_db.key2name
db_name2key = RCSR_db.name2key


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








