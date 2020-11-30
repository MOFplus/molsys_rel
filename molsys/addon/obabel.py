"""obabel addon

   This addon allows access to various features of the openbabel library for molsys mol objects.
   You must have openbabel V3.X installed and currently only non-periodic molecules are supported

   Current focus of the addon: SMILES and canonical smiles etc
   TBI: FF optimization and conformer search

"""

from openbabel import openbabel as ob
from openbabel import pybel

ob_log_handler = pybel.ob.OBMessageHandler()

# helper class for rings
class ring:

    def __init__(self, pyb_ring, mol):
        self.mol = mol
        self.pyb_ring = pyb_ring
        self.atoms = set([i for i in range(self.mol.natoms) if self.pyb_ring.IsInRing(i+1)])
        self.size = self.pyb_ring.Size()
        # test for aromaticity
        print ("detect arom in ring %s" % str(self.atoms))
        self.arom = self.pyb_ring.IsAromatic()
        print ("openabel says %s" % self.arom)
        if not self.arom:
            # still test atomtypes
            arom = True
            for a in self.atoms:
                at = self.mol.atypes[a].split("_")[0]
                print ("atom %d is %s" % (a, at))
                if not at in ("c3", "n2", "s2", "o2"):
                    arom = False
                    print ("setting arom to False")
            if arom:
                self.arom = True
        print ("Final result %s" % self.arom)
        self.ringsys = None
        return
    
class obabel:

    def __init__(self, mol, loglevel = 0):
        # set log level of global logghandler (this is a global setting!)
        ob_log_handler.SetOutputLevel(0)
        assert mol.periodic == False and mol.bcond == 0
        self._mol = mol
        # generate the pybel object 
        molstring = mol.to_string(ftype="txyz", plain=True)
        self.pybmol = pybel.readstring("txyz", molstring)
        # defaults
        self._smiles = None
        self._cansmiles = None
        return

    @property
    def smiles(self):
        if self._smiles == None:
            self._smiles = self.pybmol.write("smi")[:-2]
        return self._smiles

    @property
    def cansmiles(self):
        if self._cansmiles == None:
            self._cansmiles = self.pybmol.write("can")[:-2]
        return self._cansmiles

    def smiles2filename(self, smi):
        """converts a smiles string to smoething useable as a filename

        we replace
        * -> x
        / -> y
        \ -> z
        # -> ^
        (do we need to replace parantheses?)

        Args:
            smi (string): smiles or canonical smiles
        """
        rep = (
            ("*", "x"),
            ("/", "y"),
            ('\\', 'z'),
            ('#', '^'),
        )
        for r in rep:
            smi = smi.replace(*r)
        return smi

    def get_aromatic(self):
        """get indices of all aromatic ring atoms

        NOTE: openbabel seems not to properly detect all heteraromatic rings ... not clear to me why.
        
        Returns:
            list: atom indices
        """
        aromatic = [a.OBAtom.GetIndex() for a in self.pybmol.atoms if a.OBAtom.IsAromatic()]
        return aromatic

    def determine_rings(self):
        self.rings = []
        for r in self.pybmol.OBMol.GetSSSR():
            self.rings.append(ring(r, self._mol))
        # now determine aromatic ringsystems
        links = []
        self.arom_rings = [r for r in self.rings if r.arom]
        for r in self.arom_rings:
            print ("aromatic ring %s" % r.atoms)
        self.arom_ringsys = []
        arng_idx = 0
        largest_j = -1
        for i,r1 in enumerate(self.arom_rings):
            for ii,r2 in enumerate(self.arom_rings[i+1:]):
                j = ii+i+1
                if not r1.atoms.isdisjoint(r2.atoms):
                    # print ("link between %d and %d" % (i,j))
                    if r1.ringsys is None:
                        # this is a new ringsystem
                        new_ringsys = [i, j]
                        r1.ringsys = arng_idx
                        r2.ringsys = arng_idx
                        self.arom_ringsys.append(new_ringsys)
                        arng_idx += 1
                    else:
                        # this is an exisiting ringsystem .. find it
                        for k, rgsys in enumerate(self.arom_ringsys):
                            if i in rgsys:
                                break
                        # add j to rgsys k
                        rgsys.append(j)
                        r2.ringsys = k
                        if j > largest_j:
                            largest_j = j
        # print ("linked arom rings %s" % self.arom_ringsys)
        # now add all remaining non-connected rings to individual ringsystems
        k = largest_j+1
        for r in self.arom_rings:
            if r.ringsys is None:
                self.arom_ringsys.append([k])
                k += 1
                r.ringsys = k
        print ("all aromatic rings %s" % self.arom_ringsys)
        return

    def get_aromatic_ringsystems(self):
        self.determine_rings()
        rgsys = []
        for s in self.arom_ringsys:
            aidx = set()
            for ri in s:
                aidx |= self.arom_rings[ri].atoms
            rgsys.append(aidx)
        return rgsys 

    def plot_svg(self,fname):
        self.pybmol.write(format="svg",filename=fname,overwrite=True) 
        return
   
    def check_chirality(self):
        centers = []
        is_chiral = False
        m = self.pybmol.OBMol 
        facade = ob.OBStereoFacade(m)
        for iat in range(1,m.NumAtoms()+1):
            tetstereo = facade.GetTetrahedralStereo(m.GetAtom(iat).GetId())
            if tetstereo is not None:
                config = tetstereo.GetConfig()
            local_check = m.GetAtom(iat).IsChiral()
            centers.append(local_check)
            is_chiral = is_chiral or (local_check)
        return is_chiral, centers 


