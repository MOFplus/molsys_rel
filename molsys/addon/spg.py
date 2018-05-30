# -*- coding: utf-8 -*-
"""

    spg

    implements an addon to access the features of the spglib within molsys
    https://atztogo.github.io/spglib/

    you need a recent spglib (>= 1.9.0) because the python import has changed at some point
    it can be installed via pip (pip install spglib)

    comment by JK: i've added spacegroups.py in util containing a (an incomplete) list of spacegroup
    strings with the corresponding spacegroup numbers.
    molsys.util.spacegroups (herein imported as spacegroup) call via:
    spacegroup.get_spacegroup_number(sgname). returns None if not in dict.

Created on Wed Dec  7 15:44:36 2016

@author: rochus
"""

import spglib
import numpy
from molsys.util import elems
from molsys.util import spacegroups
import molsys
import sys
import numpy as np
import copy

import logging
logger = logging.getLogger("molsys.spg")

class spg:

    def __init__(self, mol):
        """
        generate a spg object

        :Parameters:

            - mol: mol object to be kept as a parent ref
        """
        self.mol = mol
        self.spgcell = None # tuple of (lattice, position, numbers) as used in spglib
        self.spg_version = spglib.get_version()
        self.symprec = 1.0e-2
        logger.info("Addon spg loaded (version %d.%d.%d)" % self.spg_version)
        return

    def set_symprec(self, thresh):
        """
        set the symmetry threshold
        """
        self.symprec = thresh
        return

    def get_symprec(self):
        """
        get the symmetry threshold
        """
        return self.symprec

    def generate_spgcell(self, omit=[]):
        """
        Generate the spglib specific representation of the structure
        (Needs to be called before any other call to spg methods)
         :Parameters:

             - omit : a list of either integers (atomic numbers) or element strings to be omited [optional]
        """
        # convert element symbols into atomic numbers
        if len(omit)>0:
            new_omit = []
            for e in omit:
                if type(e)==str and len(e)==0: continue
                if type(e) != type(1):
                    new_omit.append(elems.number[e.lower()])
                else:
                    new_omit.append(e)
            omit = new_omit
        lattice = numpy.array(self.mol.get_cell(), order="C", dtype="double")
        pos     = self.mol.get_frac_xyz()
        pos     = pos%1.0
        # pos     = pos%1.0
        num     = self.mol.get_elems_number()
        pos_rem = []
        num_rem = []
        for i in range(self.mol.natoms):
            if num[i] not in omit:
                pos_rem.append(pos[i])
                num_rem.append(num[i])
        pos_rem = numpy.array(pos_rem, order="C", dtype="double")
        num_rem = numpy.array(num_rem, dtype="intc")
        self.spgcell = (lattice, pos_rem, num_rem)
        return

    def get_spacegroup(self):
        """
        determine the space group of the current system
        returns a tuple with the symbol and the integer number
        """
        try:
            assert self.spgcell != None
        except:
            self.generate_spgcell()
        #print(self.spgcell)
        result = spglib.get_spacegroup(self.spgcell, symprec=self.symprec)
        result = result.split()
        symbol = result[0]
        number = int(result[1][1:-1])
        if number == 1:
            logger.warning("symmetry detection claims it's P1")
        else:
            logger.info('detected spacegroup %s %i with symprec=%5.4f' % (symbol, number, self.symprec))
        return (symbol, number)

    def make_P1(self, spgnum=-1, sg_setting=1):
        """
        to be implemented by Julian from his topo tools

        :Parameters:

            - spgnum : integer space group number

        :KNOWN BUGS:
            - scaled_positions could be equivalent from a cif file, so it fails to make_P1
        """
        # how to convert international spgnum to hall number
        # apply operations to self.mol and generate a new mol object
        # use detect_conn etc to complete it.
        
        #Okay, what i did was to use ASE as:
        try: 
            from ase.lattice.spacegroup import Spacegroup
        except:
            logger.error('make_P1 requires ASE (i.e. ase.lattice.spacegroup) to function properly')
            return
        
        # 1) if provided, use the spacegroup set by the user
        # 2) the spacegroup number is supplied with the cif
        # 3) there is at least the H-M symbol of it, try to find it via the dictionary!
        
        if spgnum == -1:
            try:
                spgnum_cif = int(self.mol.cifdata['_symmetry_int_tables_number'])
                try:
                    spgsym_cif = self.mol.cifdata['_symmetry_space_group_name_H-M'].replace(' ','')
                except:
                    sgs = spacegroups.spacegroups
                    spgsym_cif = str([i for i in sgs.keys() if sgs[i] == spgsym_cif][0])
                logger.info('using spacegroup number from cif file: %i (%s)' % (spgnum_cif,spgsym_cif))
                spgnum=spgnum_cif
            except:
                try:
                    spgsym_cif = self.mol.cifdata['_symmetry_space_group_name_H-M'].replace(' ','')
                    spgnum_cif = spacegroups.get_spacegroup_number(spgsym_cif)
                    if spgnum_cif ==None: 
                        logger.error('spacegroup %s could not be found, add it to spacegroups.py ?!' %spgsym_cif)
                        logger.error('make_P1 failed')
                        return False
                    logger.info('using spacegroup symbol from cif file, sgnumber looked up: %i (%s)' % (spgnum_cif,spgsym_cif))
                    spgnum = spgnum_cif
                except:
                    logger.error('I do not have any spacegroup informations, make_P1 failed!')
                    return False
        
        dataset = spglib.get_symmetry_from_database(spgnum)
        #print(dataset)
        
        #self.sg = Spacegroup(spgnum,setting=sg_setting)#,sprec = 1e-3) 
        self.sg = Spacegroup(spgnum,setting=sg_setting)#,sprec = 1e-3) 
        
        new_xyz = []
        new_elems = []
        new_atypes = []
        frac_xyz = self.mol.get_frac_xyz()
        #new_xyz,kinds =self.sg.equivalent_sites(frac_xyz,symprec=self.symprec)
        try:
            new_xyz,kinds =self.sg.equivalent_sites(frac_xyz,symprec=1.0e-6)
        except:
            import sys
            logger.error('could not get equivalent sites, '+str(sys.exc_info()[1]))
            return False
        #now do the new elems and stuff:
        for i,k in enumerate(kinds):
            new_elems.append(self.mol.elems[k])
            new_atypes.append(self.mol.atypes[k])
            # fragtypes = self.mol.fragtypes[k]  #### Q: work on this here? they there?
        ## now we try to get the connectivity right and find duplicates during the search
        self.mol.set_natoms(len(new_xyz))
        self.mol.set_elems(new_elems)
        self.mol.set_atypes(new_atypes)
        self.mol.set_xyz_from_frac(new_xyz)
        self.mol.set_nofrags()
        
        #self.mol.elems  = new_elems
        #self.mol.atypes = new_atypes
        self.mol.detect_conn(tresh = 0.1,remove_duplicates = True)
        return True

    def get_primitive_cell(self):
        """
        get the primitve cell as a new mol object
        """
        assert self.spgcell != None
        new_spgcell = spglib.find_primitive(self.spgcell)
        if new_spgcell is None:
            logger.error("Search for primitive cell failed with symprec %f" % self.symprec)
            return
        print(new_spgcell[0])
        print(new_spgcell[2])
        new_mol = molsys.mol()
        new_mol.set_natoms(len(new_spgcell[2]))
        new_mol.set_cell(new_spgcell[0])
        new_mol.set_xyz(new_mol.get_real_from_frac(new_spgcell[1]))
        new_mol.set_elems_number(new_spgcell[2])
        # now add the connectivity
        new_mol.detect_conn()
        # RS: we could do atomtyping ... but this would have to be a method of mol ...
        new_mol.set_atypes(["0"]*new_mol.get_natoms())
        new_mol.set_nofrags()
        return new_mol

    def get_symmetry(self):
        """
        returns lists of rotations, translations and equivalent atoms according to the spgcell
        n.b.: spgcell must be generated with generate_spgcell
        example:
        >>> import molsys
        >>> import numpy as np
        >>> m = molsys.mol()
        >>> m.read(filename)
        >>> m.addon("spg")
        >>> m.spg.generate_spgcell()
        >>> sym = m.spg.get_symmetry()
        >>> n=0 #just an example, n could be any btw. 0 and len(sym)-1
        >>> rota, tran = sym['rotations'][n], sym['translations'][n]
        >>> new_vector = rota*old_vector[:,np.newaxis] + tran
        """
        logger.info("Get symmetries")
        sym = spglib.get_symmetry(self.spgcell)
        logger.info("Found %s symmetry/ies and %s equivalent atom/s" % \
            (len(sym['rotations']), len(sym['equivalent_atoms'])))
        return sym['rotations'], sym['translations'], sym['equivalent_atoms']

    def generate_symmetries(self):
        """
        Generate list of coordinates by symmetries
        scale (same scale as per supercell) ###TBI: non-orthorombic cells
        """
        logger.info("Generating symmetries")
        self.generate_spgcell()
        lrot, ltra, leqa = self.get_symmetry() ###TBI equivalent atoms
        nsym = len(lrot)
        supercell = np.diagonal(self.spgcell[0])
        superinv = 1./supercell
        self.mol.scale_cell(superinv)
        xyzsymlist = []
        fracsymlist = []
        for i in range(nsym):
            frac = np.tensordot(self.mol.xyz, lrot[i], axes=1)+ltra[i]
            frac[np.isclose(frac,0)]=0. ###avoid negative zeros for floor
            frac[np.isclose(frac,1)]=0. ###ones are zeros
            frac[np.isclose(frac,-1)]=0. ###minus ones are zeros
            ### TBI: TO BE TRIED, WORKS FOR ANY INTEGER
            ### fround = np.rint(frac)
            ### intindex = np.isclose(frac,fround)
            ### frac[intindex]=fround[intindex]
            frac -= np.floor(frac)
            fracsymlist.append(frac)
            xyzsym = frac*supercell
            xyzsymlist.append(xyzsym)
        self.mol.scale_cell(supercell)
        self.syms = xyzsymlist
        self.fracsyms = fracsymlist

    def generate_symperms(self):
        """
        Each symmetry permutation stores the indices that would sort an array
        according to each symmetry operation in the symmetry space group.

        >>> m.addon('spg')
        >>> m.spg.generate_spgcell()
        >>> m.spg.generate_symmetries()
        >>> m.spg.generate_symperms()
        """
        logger.info("Generating symmetry permutations")
        xyzfrac = self.mol.get_frac_xyz()
        symperms = []
        for i,isym in enumerate(self.fracsyms):
            symperm = []
            for c in isym:
                frac = xyzfrac-c
                frac[np.isclose(frac,0)]=0. ###avoid negative zeros for floor
                frac[np.isclose(frac,1)]=0. ###ones are zeros
                frac[np.isclose(frac,-1)]=0. ###ones are zeros
                frac -= np.floor(frac)
                sype = np.where((frac < 1.5e-7).all(axis=1))[0]
                symperm.append(sype[0])
            self.syms[i] = frac[symperm] ###ensures equality, overcomes "float" uncertainty
            symperms.append(symperm)
        self.symperms = symperms

    def find_symmetry(self, xyzref):
        """
        If a match is found, return True. Else, return False.
        """
        logger.info("Seeking symmetry match")
        match = False
        for i,isp in enumerate(self.symperms):
            if np.isclose(self.mol.xyz[isp], xyzref).all():
                match = True
                logger.info("Find symmetry!\nIndex: %d\nPermutation: %s" % (i,isp))
                return i, isp
        if match == False:
            logger.info("No symmetry found")
            raise ValueError("No symmetry found")

    def find_symmetry_from_frac(self, fracref):
        """
        If a match is found, return True. Else, return False.
        """
        logger.info("Seeking symmetry match")
        match = False
        frac = self.mol.get_frac_xyz()
        for i,isp in enumerate(self.symperms):
            if np.isclose(frac[isp], fracref).all():
                match = True
                logger.info("Find symmetry!\nIndex: %d\nPermutation: %s" % (i,isp))
                return i, isp
        if match == False:
            logger.info("No symmetry found")
            raise ValueError("No symmetry found")

    def find_symmetry_from_colors(self, colref=None, symperms = None):
        """
        If a match is found, return True. Else, return False.
        """
        if symperms is None: symperms = self.symperms
        if colref is None: colref = self.colors
        logger.info("Seeking symmetry match")
        match = False
        col = self.mol.colors ###???
        for i,isp in enumerate(symperms): ###???
            if np.isclose(col[isp], colref).all():
                match = True
                logger.info("Find symmetry!\nIndex: %d\nPermutation: %s" % (i,isp))
                return i, isp
        if match == False:
            logger.info("No symmetry found")
            raise ValueError("No symmetry found")

    def get_frac_match(frac, sym, thresh=5e-6, eps=1e-8):
        """retrieve equivalent 
        
        [description]
        
        Arguments:
            frac {(N,3) numpy array of floats} -- 
            sym  {(N,3) numpy array of floats} -- 
        
        Keyword Arguments:
            thresh {[type]} -- [description] (default: {5e-6})
            eps {[type]} -- [description] (default: {1e-8})
        
        Returns:
            [type] -- [description]
        """
        symperm = []
        x = frac[np.newaxis,:]-sym[:,np.newaxis]
        whereint = np.where(np.isclose(np.round(x), x, atol=eps))
        x[whereint] = np.round(x[whereint]) + eps
        x -= np.floor(x) ### np.round does not work! [?RA]
        symperm = np.where((x<thresh).all(axis=-1))[-1]
        return symperm

    def generate_symmetry_dataset(self,eps=1e-13):
        """
            Set up the data necessary to exploit symmetry within weaver.
            transformations in the end contains a set of rotations and translations
            for every vertex.

        """
        self.dataset = spglib.get_symmetry_dataset(self.spgcell)
        self.RT = [(r, t) for r, t in zip(self.dataset['rotations'], self.dataset['translations'])]
        self.R = self.dataset['rotations']
        self.T = self.dataset['translations']
        equiv = self.dataset['equivalent_atoms'].tolist()
        equiv_set = list(set(equiv))
        base_indices = [equiv.index(i) for i in equiv_set]
        derived_indices = [[i for i,e in enumerate(equiv) if ((e == j) and (i != -1))] for j in equiv_set]
        
        transformations = {}
        frac = copy.copy(self.spgcell[1])
        for ie,e in enumerate(equiv_set):
            exyz = frac[e]
            M_exyz = (numpy.dot(self.R,exyz)+self.T) #% 1.0
            for id,d in enumerate(derived_indices[ie]):
                dxyz  = frac[d]
                distvects = M_exyz-dxyz
                whereint = np.where(np.isclose(np.round(distvects), distvects, atol=eps))
                distvects[whereint] = np.round(distvects[whereint]) + eps
                distvects -= numpy.floor(distvects)
                
                dists = numpy.linalg.norm(distvects,axis=1)
                idx_dmin  = [i for i,dist in enumerate(dists) if dist < 1e-5][0]
                #idx_dmin = numpy.argmin(dists)
                if dists[idx_dmin] > 1e-8:
                    #import pdb; pdb.set_trace()
                    raise ValueError('no transformation found! for %i ' % (d,))
                #print id,d,numpy.min(dists), numpy.argmin(dists)
                transformations[d] = self.RT[idx_dmin]
        self.transformations = transformations
        self.equiv = equiv
        self.equiv_set = equiv_set
        self.base_indices = base_indices
        self.derived_indices = derived_indices
        
        return transformations