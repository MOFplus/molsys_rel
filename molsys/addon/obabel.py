"""obabel addon

   This addon allows access to various features of the openbabel library for molsys mol objects.
   You must have openbabel V3.X installed and currently only non-periodic molecules are supported

   Current focus of the addon: SMILES and canonical smiles etc
   TBI: FF optimization and conformer search

"""

from openbabel import openbabel as ob
from openbabel import pybel


class obabel:

    def __init__(self, mol):
        assert mol.periodic == False and mol.bcond == 0
        self._mol = mol
        # generate the pybel object 
        molstring = mol.to_string(ftype="txyz", plain=True)
        self.pybmol = pybel.readstring("txyz", molstring)
        return

    def get_SMILES(self):
        return self.pybmol.write("smi")

    def get_CANSMILES(self):
        return self.pybmol.write("can")