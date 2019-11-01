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
from molsys.util.images import idx2arr,arr2idx

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
                        # NOTE: in some cases pconn contains FLOAT numbers which is wrong!!! find out where thsi comes from!! and who did it?
                        labels.append(self._mol.pconn[i][j].astype("int32").tolist())
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

    def extend_images(self):
        """helper method to generate a new mol object where the periodic bonds are extended
        to the next image as defined by the pconn
        """
        nm = self._mol.clone()
        nm.unmake_topo()
        del_bonds = []
        for i in range(self._mol.get_natoms()):
            for ij,j in enumerate(self._mol.conn[i]):
                if arr2idx[self._mol.pconn[i][ij]] != 13:
                    # this is a periodic bond .. add image atom
                    xyz_image_j = self._mol.xyz[j]+(self._mol.cell*(self._mol.pconn[i][ij][:,None])).sum(axis=0)
                    new_atom = nm.add_atom("x", nm.get_elems()[j], xyz_image_j)
                    if j>i:
                        del_bonds.append((i,j))
                    nm.add_bond(i, new_atom)
        for b in del_bonds:
            print ("delete bond %s " % str(b))
            nm.delete_bond(b[0], b[1])
        return nm