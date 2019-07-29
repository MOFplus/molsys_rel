
import string
import numpy
import copy
import molsys.mol
from molsys.addon import base
# molecules module

# atm it is just those routines that are removed from the old molsys stored here as backup
import logging
logger = logging.getLogger("molsys.molecules")

class mgroup:
    """class to keep a set of groups/molecules

    any atom must be in exactly one group/molecule
    Note that we use the names molecules and groups synonymously
    """

    def __init__(self, mol):
        self.whichmol = [] # a listof length natoms -> index of the molecule
        self.mols     = [] # a list of lists (length nmols) with the atom indices
        self.moltypes = [] # a list of length nmols with indices of molnames
        self.molnames = [] # a list of molnames
        self.nmols = 0     # len of moltypes/mols
        self.setup = False
        self.parent_mol = mol
        return

    def detect_molecules(self):
        ''' Detects independent (not connected) fragments and stores them as 
        '''
        # do not do this if already setup
        assert self.setup is False
        # the default moleculename is taken from the parent molfile
        self.molnames = [self.parent_mol.name]
        atoms = range(self.parent_mol.natoms)
        self.whichmol = self.parent_mol.natoms * [0]
        nmol = 0
        while len(atoms) > 0:
            # localize one molecule starting from the first available atom
            leafs = [atoms[0]]
            curr_mol = []
            while len(leafs) > 0:
                new_leafs = []
                # add all to curr_mol, remove from atoms and generate new_leafs
                for l in leafs:
                    atoms.remove(l)
                    curr_mol.append(l)
                    new_leafs += self.parent_mol.conn[l]
                # first remove duplicates in new_leafs
                for l in copy.copy(new_leafs):
                    i = new_leafs.count(l)
                    if i > 1:
                        for j in range(i-1):
                            new_leafs.remove(l)
                # now cut new_leafs (remove all those we already have in curr_mol)
                for l in copy.copy(new_leafs):
                    if curr_mol.count(l):
                        new_leafs.remove(l)
                # now make new_leafs to leafs and continue
                leafs = new_leafs
            # at this point the molecule is complete
            curr_mol.sort()
            self.mols.append(curr_mol)
            for i in curr_mol:
                self.whichmol[i] = nmol
            # at this point all molecules found get the type 0 = "xyz"
            # if len(self._molecules.keys()) != 0:
            #     #for latest GCMD version this needs to be done here
            #     # not at its final beazty here ... 
            #     if nmol != 0:
            #         self.moltypes.append(1)
            #     else:
            #         self.moltypes.append(0)
            # else:
            self.moltypes.append(0)
            nmol += 1
        # all atoms are assigned
        self.nmols = nmol
        # self.molnames += self.get_names()
        return

    def add_molecules(self, name, n, natoms):
        """add a number of molecules to this mgroup

        this implies that the molecules have been added to the parent_mol already
        
        Args:
            name (string): name of the molecule
            n (integer): number of molecules to be added
            natoms (integer): number of atoms in one molecule
        """
        orig_natoms = len(self.whichmol)
        self.molnames.append(name)
        mtype = len(self.molnames)-1
        for i in xrange(n):
            # add n molecules
            self.moltypes.append(mtype)
            self.mols.append(range(orig_natoms+i*natoms,orig_natoms+(i+1)*natoms))
            self.whichmol += natoms*[self.nmols+i]
        self.nmols = len(self.mols)
        assert len(self.whichmol) == self.parent_mol.natoms
        return


class molecules(base):

    def __init__(self, mol):
        super(molecules,self).__init__(mol)
        self.mgroups = {}
        self.mgroups["molecules"] = mgroup(self._mol)
        self.mgroups["molecules"].detect_molecules()
        self.default_mgroup = "molecules"
        return

    def __call__(self):
        '''
            keep deprecated call method for legacy reasons and map it onto detect_molecules'
        '''
        return


    def add_molecule(self, newmol, nmols=1):
        """Adds a molecule to the parent mol object

        A molecules is a sub mol object that can be added to the parent system many times

        Args:
            - mol (molsys obeject): nonperiodic molecule to be added
        """
        from molsys.addon.ff import ic
        # add newmol to parent mol
        offset = self._mol.natoms -1
        for i in range(nmols):
            rndxyz = self._mol.get_xyz_from_frac(numpy.random.uniform(0,1,(3,)))
            self._mol.add_mol(newmol,translate=rndxyz)
        # now add to mgroups
        for mg in self.mgroups:
            self.mgroups[mg].add_molecules(newmol.name, nmols, newmol.natoms)
        # if parent mol has a ff then the newmol needs it too
        if "ff" in self._mol.loaded_addons:
            assert "ff" in newmol.loaded_addons, "Added molecule object needs to have a ff setup"
            ### par update
            for k in newmol.ff.par.keys():
                self._mol.ff.par[k].update(newmol.ff.par[k])
            logger.info('ff information of molecule %s red' % (newmol.name,))
            ### ric_type update
            temp_ric_type = copy.copy(newmol.ff.ric_type)
            for k in temp_ric_type.keys():
                for moli in range(nmols):    
                    rictype = []
                    for i in range(len(temp_ric_type[k])):
                        for j in range(len(temp_ric_type[k][i])):
                            rictype.append(temp_ric_type[k][i][j] + offset + moli+1)
                    if rictype != []:
                        self._mol.ff.ric_type[k] += [ic(rictype)]
                        self._mol.ff.parind[k] += newmol.ff.parind[k]
        return




 