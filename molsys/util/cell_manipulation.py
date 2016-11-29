#### cell manipulations

def extend_cell(mol,offset):
    ''' Atoms as close as offset to the box boundaries are selected to be copied.
        They are then added at the other side of the cell to "extend" the system periodically 
        Mainly for visualization purposes
        WARNING: Connectivity is destroyed afterwards 
        :Params: 
            - offset: The distance (in Angstroms) from the box boundary at which to duplicate the atoms
    '''
    logging.warning('connectivity is destroyed')
    frac_xyz = mol.get_frac_xyz()
    wherexp = np.where(np.less(frac_xyz[:,0], offset))
    wherexm = np.where(np.greater(frac_xyz[:,0], 1.0-offset))
    whereyp = np.where(np.less(frac_xyz[:,1], offset))
    whereym = np.where(np.greater(frac_xyz[:,1], 1.0-offset))
    wherezp = np.where(np.less(frac_xyz[:,2], offset))
    wherezm = np.where(np.greater(frac_xyz[:,2], 1.0-offset))
    new_xyz = frac_xyz
    #print new_xyz.shape
    new_xyz = np.append(new_xyz, frac_xyz[wherexp[0]]+[1.0,0.0,0.0],0)
    new_xyz = np.append(new_xyz, frac_xyz[whereyp[0]]+[0.0,1.0,0.0],0)
    new_xyz = np.append(new_xyz, frac_xyz[wherezp[0]]+[0.0,0.0,1.0],0)
    new_xyz = np.append(new_xyz, frac_xyz[wherexm[0]]-[1.0,0.0,0.0],0)
    new_xyz = np.append(new_xyz, frac_xyz[whereym[0]]-[0.0,1.0,0.0],0)
    new_xyz = np.append(new_xyz, frac_xyz[wherezm[0]]-[0.0,0.0,1.0],0)
    #print new_xyz
    #print new_xyz.shape
    mol.set_xyz_from_frac(new_xyz)
    for i in range(len(wherexp[0])):
        mol.elems.append(mol.elems[wherexp[0][i]])
    for i in range(len(whereyp[0])):
        mol.elems.append(mol.elems[whereyp[0][i]])
    for i in range(len(wherezp[0])):
        mol.elems.append(mol.elems[wherezp[0][i]])
    for i in range(len(wherexm[0])):
        mol.elems.append(mol.elems[wherexm[0][i]])
    for i in range(len(whereym[0])):
        mol.elems.append(mol.elems[whereym[0][i]])
    for i in range(len(wherezm[0])):
        mol.elems.append(mol.elems[wherezm[0][i]])
    #print new_xyz
    mol.natoms = len(mol.xyz)
    #logging.info('Cell was extended by %8.4f AA in each direction' % (offset))
    return