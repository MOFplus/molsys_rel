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


class lqg:

    def __init__(self, edges, labels):
        self.edges = edges
        self.labels = labels
        self.is_systrekey = None
        self.systrekey = None
        self.ne = len(edges)
        nv = 0
        for e in self.edges:
            for v in e:
                if v > nv:
                    nv = v
        self.nv = nv+1
        return

    @classmethod
    def from_string(cls, lqgs):
        lqg = lqgs.split()
        assert len(lqg)%5 == 0
        ne = len(lqg)//5
        edges = []
        labels = []
        for i in range(ne):
            edges.append([int(lqg[i*5])-1, int(lqg[i*5+1])-1])
            labels.append([int(lqg[i*5+2]), int(lqg[i*5+3]), int(lqg[i*5+4]),])
        return cls(edges, labels)

    def __repr__(self):
        out = ""
        for e,l in zip(self.edges, self.labels):
            out += "%d %d %d %d %d " % (e[0]+1, e[1]+1, l[0], l[1], l[2])
        return out

    def get_edge(self, i, j):
        if i < j:
            e = [i, j]
        else:
            e = [j, i]
        if not (self.edges.count(e) == 1):
            print (e)
            print (self.edges)

        ind = self.edges.index(e)
        return ind, self.edges[ind], self.labels[ind]

    def get_systrekey(self):
        """run javascript systrekey
        
        revised version using json strings to pass data to javascript
        """
        import json
        if self.is_systrekey:
            return self
        if self.systrekey is not None:
            return self.systrekey
        # the systrekey is not yet available -> compute it
        if not node_avail:
            # return nothing to avoid an error
            print("""
            WARNING: systrekey was called but node is not installed to run javascript
            """
            )
            return
        lqg_in = []
        for e,l in zip(self.edges, self.labels):
            el = []
            el.append(int(e[0])+1)
            el.append(int(e[1])+1)
            el.append(l)
            lqg_in.append(el)
        json_lqg = str(lqg_in)
        print (json_lqg)
        try:
            json_result = subprocess.check_output(args=["node", molsys_path+"/util/run_systrekey.js", json_lqg, systre_path], stderr=subprocess.STDOUT).decode()
        except subprocess.CalledProcessError as err:
            raw_err = err.stdout.decode().split("\n")
            for el in raw_err:
                if el[:6] == "Error:":
                    err = el
                    break
            print ("ERROR: systrekey says -> %s" % err)
            return
        result = json.loads(json_result)
        print (result)
        skey = result["key"][2:] # cut off the "3 " for the 3D
        self.systrekey = lqg.from_string(skey)
        self.systrekey.is_systrekey = True 
        mapping = result["mapping"]
        self.sk_mapping = {}
        for k in mapping:
            self.sk_mapping[int(k)-1] = mapping[k]-1
        # now let us try to map the edges from the
        self.sk_edge_mapping = {}
        for i, e in enumerate(self.edges):
            si = self.sk_mapping[e[0]]
            sj = self.sk_mapping[e[1]]
            self.sk_edge_mapping[i] = self.systrekey.get_edge(si, sj)[0]
        return self.systrekey


if __name__=="__main__":
    edges = [[0,0], [0,0], [0,0]]
    labels = [[1,0,0], [0,1,0], [0,0,1]]
    g = lqg(edges, labels)
    print (g.get_systrekey())








