import copy
import string
import numpy
import molsys.mol
# molecules module

# atm it is just those routines that are removed from the old molsys stored here as backup
import logging
logger = logging.getLogger("molsys.molecules")


class molecules(object):

    def __init__(self, mol):
        self._mol = mol
        self._molecules = {}
        self._nmols = {}
        # self.detect_molecules()
        return

    def get_names(self):
        return self._molecules.keys()

    def add_molecule(self, molname, nmols=1, ff=False):
        """Adds a molecule to the parent mol object

        A molecules is a sub mol object that can be added to the parent system many times

        Args:
            molname (string): name of the molecule
            fname (string, optional): Defaults to None. if the filename is different to the name of the molecule,

        """
        from molsys.addon.ff import ic
        newmol = molsys.mol.from_file(molname)
        self._molecules.update({molname: newmol})
        self._nmols.update({molname:nmols})

        offset = self._mol.natoms -1
        for i in range(nmols):
            rndxyz = self._mol.get_xyz_from_frac(numpy.random.uniform(0,1,(3,)))
            self._mol.add_mol(newmol,translate=rndxyz)
        #self.nmols += nmols

        if ff is not False:
            if not hasattr(self._mol,'ff'): raise AssertionError('mol instance has no ff loaded!')
            newmol.addon('ff')
            newmol.ff.read(molname)
            ### par update
            for k in newmol.ff.par.keys():
                self._mol.ff.par[k].update(newmol.ff.par[k])
            logger.info('ff information of molecule %s red' % (molname,))
            ### ric_type update
            #import pdb; pdb.set_trace()
            temp_ric_type = copy.copy(newmol.ff.ric_type)
            for k in temp_ric_type.keys():
                for moli in range(nmols):    
                    rictype = []
                    for i in range(len(temp_ric_type[k])):
                        for j in range(len(temp_ric_type[k][i])):
                            rictype.append(temp_ric_type[k][i][j] + offset + moli+1)
                    #print k, moli, rictype
                    if rictype != []:
                        #rictype = [rictype]
                        self._mol.ff.ric_type[k] += [ic(rictype)]
                        self._mol.ff.parind[k] += newmol.ff.parind[k]
            ### parind update 
            for k in newmol.ff.parind.keys():
                pass
        return

    def detect_molecules_new(self):
        """
            Detect molecules by connectivity
            atm only what's absolutely necessary
        """
        self.molecules = []
        self.whichmol = []
        self.molnames = []
        self.moltypes = []
        self.nmols = 0
        for ic, c in enumerate(self._mol.get_conn()):
            if c == []:
                self.nmols += 1
                self.molecules.append('')
            else:
                pass
            self.whichmol.append(self.nmols)
            # WORK HERE!

    def __call__(self):
        '''
            keep deprecated call method for legacy reasons and map it onto detect_molecules'
        '''
        self.detect_molecules()
        return

    def detect_molecules(self):
        ''' Detects independent (not connected) fragments and stores them as 
            - mols     :
            - moltypes :
            - whichmol :
        '''
        self.mols = []
        self.moltypes = []
        # the default moleculename is "xyz" -> molecules from the xyz file
        self.molnames = ["xyz"]
        atoms = range(self._mol.natoms)
        self.whichmol = self._mol.natoms * [0]
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
                    new_leafs += self._mol.conn[l]
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
            # quick hack to get one type of molecule to run
            if nmol != 0:
                self.moltypes.append(1)
            else:
                self.moltypes.append(0)
            nmol += 1
        # all atoms are assigned
        # if mol.verbose:
        #print("$$ -- found %d independent molecules from connectivity" % nmol)
        self.nmols = nmol
        self.molnames += self.get_names()
        return

    ### for use with pydlpoly's pdlpmol
    def to_pdlpmol(self,pd):
        """fill the pdlpmol instance with the data from the molecuels addon
        
        Fills the date from the molecules addon into the pdlpmol instance from pydlpoly

        
        Args:
            pdlpmol (pdlpmol.pdlpmol): pdlpmol instance
        """
        pdmol = pd.mol
        pdmolecules = pd.mol.molecules
        for i,name in self.get_names():
            moltype = len(pdmol.molnames)
            pdmol.molnames += [name]
            pdmol.moltypes.append([moltype for i in self.nmols[name]])
        
        
        return

    # that one is from assign_FF!
    def add_mol_tinker_xyz(self, fname, N, molname, offset, scale, rotate=True):

        moltype = len(self.molnames)
        self.molnames.append(molname)
        # read it in first into some temp arrays
        #f = open(fname, "r")
        #lbuffer = string.split(f.readline())
        #mna = string.atoi(lbuffer[0])
        #mxyz = []
        #melems = []
        #mtypes = []
        #mcnct = []
        #for i in xrange(mna):
        #    lbuffer = string.split(f.readline())
        #    mxyz.append(map(string.atof, lbuffer[2:5]))
        #    melems.append(string.lower(lbuffer[1]))
        #    t = lbuffer[5]
        #    mtypes.append(t)
        #    if not self.typedata.has_key(t): self.typedata[t] = None
        #    mcnct.append(num.array(map(string.atoi, lbuffer[6:]))-1)
        #f.close()
        #self.ntypes = len(self.typedata.keys())
        # center
        #amass = []
        #for e in melems: amass.append(atomicmass[e])
        #amass = num.array(amass)
        #mxyz = num.array(mxyz)
        #com = sum(mxyz*amass[:,num.newaxis],0)/sum(amass)
        #mxyz -= com
        # now generate N translated/rotated copies
        #if self.verbose: print  ("$$ -- generating %d copies of the molecule in the box"%N)
        #cell = num.array(self.cell)
        for i in xrange(N):
        #    act_mxyz = mxyz.copy()
        #    if rotate:
        #        # rotate molecule by random quaternion
        #        act_mxyz = rotate_random(act_mxyz)
        #    # translate
        #    #act_gxyz += cell*num.random.random(3)
        #    rand_vect = offset+(num.random.random(3)*scale)
        #    act_mxyz += num.sum(cell*rand_vect,axis=0)
            # register moltype (per molecule) and to which mol the atom belongs (per atom)
            self.moltypes.append(moltype)
            current_mol = self.nmols+i
            for j in xrange(mna):
                self.xyz.append(act_mxyz[j].tolist())
                self.elems.append(melems[j])
                self.types.append(mtypes[j])
                self.cnct.append((mcnct[j]+self.natoms).tolist())
                self.whichmol.append(current_mol)
            self.mols.append(range(self.natoms,self.natoms+mna))
            self.natoms += mna
        self.nmols += N
        # update exclude list
        old_natoms = len(self.exclude)
        self.exclude += (self.natoms-old_natoms)*[False]
        self.frozen  += (self.natoms-old_natoms)*[0]        
        return
    
