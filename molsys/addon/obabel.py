"""obabel addon

   This addon allows access to various features of the openbabel library for molsys mol objects.
   You must have openbabel V3.X installed and currently only non-periodic molecules are supported

   Current focus of the addon: SMILES and canonical smiles etc
   TBI: FF optimization and conformer search

"""

from openbabel import openbabel as ob
from openbabel import pybel

ob_log_handler = pybel.ob.OBMessageHandler()

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
            self._smiles = self.pybmol.write("smi")
        return self._smiles

    @property
    def cansmiles(self):
        if self._cansmiles == None:
            self._cansmiles = self.pybmol.write("can")
        return self._cansmiles