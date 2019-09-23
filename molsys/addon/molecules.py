
import string
import numpy
import copy
import molsys.mol
from molsys.addon import base

import os
import tempfile
#from molsys.molsys_mpi import mpiobject
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
        self.mols     = [] # a list of lists (length nmols) with the atom indices // OPTIONAL
        self.moltypes = [] # a list of length nmols with indices of molnames
        self.molnames = [] # a list of molnames
        self.nmols = 0     # len of moltypes/mols
        self.setup = False
        self.parent_mol = mol
        return

    def detect_molecules(self):
        """ Detect independent (not connected) molecules and store them

        We try to use first the molid information if it is present
        This has been generated by label_components from graphtools during assignement

        If this information is not present then we use the simple python fallback routine which is very slow and
        takes long for really large (>10000 atoms) systems
        """
        assert self.setup is False # in the future just return ... just to see if it is called twice errouneously
        # check if we have a molid
        if self.parent_mol.molid is not None:
            # ok, we can go and set things up
            self.whichmol = list(self.parent_mol.molid)
            # now we need to construct the other data from whichmol
            self.nmols = self.parent_mol.molid.max()+1
            # since we do not know molnames if this is a pure mfpx file we give them for the time being all the name of the parent mol obejct
            self.molnames = [self.parent_mol.name]
            self.moltypes = self.nmols*[0]
            # generate mols (list of lists)
            for i in range(self.nmols):
                self.mols.append([])  # nmols*[[]] does not work -> side effect
            for i, m in enumerate(self.whichmol):
                self.mols[m].append(i)
        else:
            self.detect_molecules_nograph()
        return

    def detect_molecules_nograph(self):
        ''' Detects independent (not connected) fragments and stores them

        This is the fallback procedure if the mol object does not already contain a molid

        '''
        # the default moleculename is taken from the parent molfile
        self.molnames = [self.parent_mol.name]
        atoms = list(range(self.parent_mol.natoms))
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
        # this has been called since the parents molid is empty => write whichmol there as a numpy array
        self.parent_mol.molid = numpy.array(self.whichmol)
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
        for i in range(n):
            # add n molecules
            self.moltypes.append(mtype)
            self.mols.append(list(range(orig_natoms+i*natoms,orig_natoms+(i+1)*natoms)))
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

    def __getattr__(self,name):
        ''' Overridden for temporary back compatibility
        
        Since both, pydlpoly and pylmps rely on mol.molecules.[lots of attributes], we map here
        the attributes of the default molecule group to the molecules instanace.
        may interfere at some point when mgroups[molecules] and this class have same named attributes
        currently this is not the case.

        Args:
            self (base): molecules instance
            name (str): attribute to be obtainde from the underlying mgroup
        
        Returns:
            varying: whatever the value of name is
        '''
        try:
            return object.__getattribute__(self,name)
        except:
            return getattr(self.mgroups[self.default_mgroup],name)

    def add_molecule(self, newmol, nmols=1, pack=False, packbound=0.5,ff=None):
        """Adds a molecule to the parent mol object

        A molecules is a sub mol object that can be added to the parent system many times

        Args:
            - mol (molsys obeject): nonperiodic molecule to be added (must be ONE connected molecule)
            - pack (boolean): defaults to False: if True use packmol to pack (needs to be installed and in the path)
            - packbound (float): defaults to 0.5: distance in Angstrom to reduce the box for packmol filling
            - ff (None):  dummy variable added for legacy reasons
        """
        
        # temporary hook to get legacy scripts running (JK)
        if type(newmol) == type('string'):
            name = newmol
            newmol = molsys.mol.from_file(name+'.mfpx')
            newmol.addon('ff')
            newmol.ff.read(name)
        
        # if pack is True we first try to use packmol to generate proper xyz coords and keep them.        
        if pack:
            assert self._mol.get_bcond() < 3
            # assume pydlpoly boundary conditions (orig in the center of box) -- this is what we get from MOF+ 
            cell = self._mol.get_cell().diagonal()
            cellh = (cell*0.5)-packbound
            box = (-cellh).tolist()
            box += cellh.tolist()
            # make a temp file and go there
            tmpd = tempfile.mkdtemp()
            cwd  = os.getcwd()
            os.chdir(tmpd)
            # write parent periodic mol as xyz
            self._mol.write("self.xyz")
            # write to be added mol
            newmol.write("newmol.xyz")
            packmolf = """
                tolerance 2.0
                output final.xyz
                filetype xyz

                structure self.xyz
                    filetype xyz
                    number 1
                    fixed 0.0 0.0 0.0 0.0 0.0 0.0
                end structure

                structure newmol.xyz
                    filetype xyz
                    number %d
                    inside box %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f
                end structure
            """ % tuple([nmols]+box)
            with open("pack.inp", "w") as packf:
                packf.write(packmolf)
            # now try to execute packmol
            os.system("packmol <pack.inp > pack.out")
            # check if final.xyz has been generated ... if not it probably failed becasue packmol is not installed or other reasons
            if os.path.isfile("final.xyz"):
                tempm = molsys.mol.from_file("final.xyz")
                pack_xyz = tempm.get_xyz()
                del(tempm)
                remove_files = ["self.xyz", "newmol.xyz", "pack.inp", "pack.out", "final.xyz"]
                for f in remove_files:
                    os.remove(tmpd+"/"+f)
                os.rmdir(tmpd)
                logger.info("Molecules packed with packmol")
            else:
                self._mol.pprint("packing with packmol failed. Temporary files at %s" % tmpd)
                logger.warning("packing with packmol failed. Temporary files at %s" % tmpd)
                pack = False
            os.chdir(cwd)
        # TBI verify that the added molecule is really only one molecule ... 
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
            logger.info('ff information of molecule %s read in' % (newmol.name,))
            ### ric_type update
            temp_ric_type = copy.copy(newmol.ff.ric_type)
            for k in temp_ric_type.keys():
                for moli in range(nmols):    
                    for i in range(len(temp_ric_type[k])):
                        # get attributes of ric
                        atrs = temp_ric_type[k][i].__dict__
                        rictype = []
                        for j in range(len(temp_ric_type[k][i])):
                            rictype.append(temp_ric_type[k][i][j] + offset + (moli*newmol.natoms) +1)
                        if rictype != []:
                            new_ric = ic(rictype)
                            for a in atrs.keys():
                                if (a != "used") and (a != "type"):
                                    new_ric.__dict__[a] = atrs[a]
                            self._mol.ff.ric_type[k] += [new_ric]
                        self._mol.ff.parind[k] += newmol.ff.parind[k]
            ### molid has been updated using self._mol.add_mol => update the molid in ff
            self._mol.ff.update_molid()
        # finish up packing
        if pack:
            self._mol.set_xyz(pack_xyz)
        return

    def get_names(self,legacy=True):
        """get a list of the molecule names stored in the default molecules dictionary
        
        Args:
            legacy (bool, optional): Defaults to True. if True, remove the name of the parent mol
        
        Returns:
            list: list of strings of molecule names
        """
        namelist = list(set(self.mgroups[self.default_mgroup].molnames))
        if legacy is True:
            namelist.remove(self.mgroups[self.default_mgroup].parent_mol.name)
        return namelist








 
