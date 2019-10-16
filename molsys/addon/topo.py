"""topo is an addon for topo mol objects (we assert this at the init)

By adding this addon to a topo mol object you get a number of additional features relevant for working with topologies

Features:
    - lqg: adds the labeled quotient graph object to the mol object
        -> generate a barycentreic embedding
        -> compute the systre key and detect the name of the topo

TBI: we could think of a way to make a mol/topo object from the systrekey (or a lqg) which should directly convert the
     mol object into a topo and add the addon
"""

# from molsys.util.lqg import lqg
from molsys.util import systrekey

# make a systre key database 
skdb = systrekey.systre_db()

class topo:

    def __init__(self, mol):
        assert mol.is_topo
        self._mol = mol
        # we need pconn here in any case for many things so add it if it is not there, yet
        if not self._mol.use_pconn:
            self._mol.add_pconn()
        self._systrekey = None
        self._RCSRname   = None
        return
    
    @property
    def systrekey(self):
        """method calls javascript systreKey (must be installed) and computes the systreKey
        """
        if self._systrekey is None:
            edges = []
            labels = []
            for i in range(self._mol.get_natoms()):
                for j,v in enumerate(self._mol.conn[i]):
                    if v > i:
                        edges.append([i,v])
                        labels.append(list(self._mol.pconn[i][j].astype("int")))
            self._systrekey, mapping = systrekey.run_systrekey(edges, labels)
            self._skey_mapping = [-1 for i in range(self._mol.get_natoms())]
            for k in mapping:
                self._skey_mapping[int(k)-1] = mapping[k]
        return self._systrekey

    @property
    def RCSRname(self):
        if self._RCSRname is None:
            self._RCSRname = skdb.get_name(self.systrekey)
        return self._RCSRname