#### molecules module

###### atm it is just those routines that are removed from the old molsys stored here as backup


def find_molecules(mol):
    ''' Detects independent (not connected) fragments and stores them as 
        - mols     :
        - moltypes :
        - whichmol :
        '''
    mol.mols = []
    mol.moltypes = []
    # the default moleculename is "xyz" -> molecules from the xyz file
    mol.molnames = ["xyz"]
    atoms = range(mol.natoms)
    mol.whichmol = mol.natoms * [0]
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
                new_leafs += mol.conn[l]
            # first remove duplicates in new_leafs
            for l in copy.copy(new_leafs):
                i = new_leafs.count(l)
                if i > 1:
                    for j in xrange(i-1):
                        new_leafs.remove(l)
            # now cut new_leafs (remove all those we already have in curr_mol)
            for l in copy.copy(new_leafs):
                if curr_mol.count(l): new_leafs.remove(l)
            # now make new_leafs to leafs and continue
            leafs = new_leafs
        # at this point the molecule is complete
        curr_mol.sort()
        mol.mols.append(curr_mol)
        for i in curr_mol: mol.whichmol[i] = nmol
        # at this point all molecules found get the type 0 = "xyz"
        mol.moltypes.append(0)
        nmol += 1
    # all atoms are assigned
    #if mol.verbose:
    #print "$$ -- found %d independent molecules from connectivity" % nmol
    mol.nmols=nmol
    return