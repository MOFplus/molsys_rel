# -*- coding: utf-8 -*-
import string as st
import numpy as np
import types
import copy
import string
import logging

#import util.elems as elements
#import util.unit_cell as unit_cell
#import util.rotations as rotations
#from molsys.io import formats
from io import formats
#import io.formats
import random



try:
    from ase import Atoms
    from pyspglib import spglib
except ImportError:
    spg = False
else:
    spg = True

"""

        ToDo list ...

        - Complete rework, main molsys does not have pconn any more
        - make_supercell_preserve_conn has to be rewritten not to work with pconn ;)
        - ...

"""

np.set_printoptions(threshold=20000)

deg2rad = np.pi/180.0
SMALL_DIST = 1.0e-3

class molsys:

    def __init__(self):
        self.natoms=0
        self.cell=None
        self.cellparams=None
        self.images_cellvec=None
        self.xyz=None
        self.elems=[]
        self.atypes=[]
        self.conn=[]
        self.fragment=[]
        self.fragment_number=[]
        self.periodic= None
        logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)
        return

    ######  I/O stuff ############################
    
    def read(self,fname,ftype='mfpx',**kwargs):
        formats.read[ftype](self,fname,**kwargs)
        return
    
    def write(self,fname,ftype='mfpx',**kwargs):
        formats.write[ftype](self,fname,**kwargs)
        return
        

    ###### helper functions #######################
        
    def get_elemlist(self):
        el = []
        for e in self.elems:
            if not el.count(e): el.append(e)
        return el
        
    def get_atypelist(self):
        if not self.atypes: return None
        at = []
        for a in self.atypes:
            if not at.count(a): at.append(a)
        return at
    
    def get_distvec(self, i, j):
        """ vector from i to j
        This is a tricky bit, because it is needed also for distance detection in the blueprint
        where there can be small cell params wrt to the vertex distances.
        In other words: i can be bonded to j multiple times (each in a different image)
        and i and j could be the same!! """
        ri = self.xyz[i]
        rj = self.xyz[j]
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - ri
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            d_sort = np.argsort(all_d)
            if i == j:
                # if this was requested for i==j then we have to eliminate the shortest 
                # distance which NOTE unfinished!!!!!!!!
                pass 
            closest = d_sort[0]
            closest=[closest]  # THIS IS A BIT OF A HACK BUT WE MAKE IT ALWAYS A LIST ....
            if (abs(all_d[closest[0]]-all_d[d_sort[1]]) < SMALL_DIST):
                # oops ... there is more then one image atom in the same distance
                #  this means the distance is larger then half the cell width
                # in this case we have to return a list of distances
                for k in d_sort[1:]:
                    if (abs(all_d[d_sort[0]]-all_d[k]) < SMALL_DIST):
                        closest.append(k)
            d = all_d[closest[0]]
            r = all_r[closest[0]]
        else:
            if i == j: return
            r = rj-ri
            d = np.sqrt(np.sum(r*r))
            closest=[0]
        return d, r, closest

    # thes following functions rely on an exisiting connectivity conn (and pconn)

    def get_neighb_coords(self, i, ci):
        """ returns coordinates of atom bonded to i which is ci'th in bond list """
        j = self.conn[i][ci]
        rj = self.xyz[j].copy()
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - self.xyz[i]
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            closest = np.argsort(all_d)[0]
            return all_rj[closest]
        return rj

    def get_neighb_dist(self, i, ci):
        """ returns coordinates of atom bonded to i which is ci'th in bond list """
        ri = self.xyz[i]
        j = self.conn[i][ci]
        rj = self.xyz[j].copy()
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - self.xyz[i]
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            closest = np.argsort(all_d)[0]
            return all_rj[closest]
        dr = ri-rj
        d = np.sqrt(np.sum(dr*dr))
        return d

    ######## manipulations in particular for blueprints

    def make_supercell(self, supercell,shift=False):
        supercell = np.array(supercell)
        cellfact = supercell.prod()
        new_xyz = []
        for a in xrange(supercell[0]):
            for b in xrange(supercell[1]):
                for c in xrange(supercell[2]):
                    if shift:
                        dispvect = np.sum(self.cell*np.array([a-1,b-1,c-1])[:,np.newaxis],axis=0)
                    else:
                        dispvect = np.sum(self.cell*np.array([a,b,c])[:,np.newaxis],axis=0)
                    new_xyz += (self.xyz+dispvect).tolist()
        self.xyz = np.array(new_xyz, "d")
        self.cell *= supercell[:,np.newaxis]
        self.cellparams[0:3] *= supercell
        self.natoms*=cellfact
        self.elems *=cellfact
        self.atypes*=cellfact
        self.images_cellvec = np.dot(images, self.cell)
        return
        
    def make_supercell_symmetric(self,supercell,shift=False):    
        supercell = np.array(supercell)
        cellfact = supercell.prod()
        new_xyz = []
        dispvect = np.sum(self.cell*np.array([0,0,0])[:,np.newaxis],axis=0)
        new_xyz += (self.xyz+dispvect).tolist()
        #print len(self.xyz)
        for a in xrange(supercell[0]):
            for b in xrange(supercell[1]):
                for c in xrange(supercell[2]):
                    if (a == 1) and (b == 1) and (c==1): continue
                    if shift:
                        dispvect = np.sum(self.cell*np.array([a-1,b-1,c-1])[:,np.newaxis],axis=0)
                    else:
                        dispvect = np.sum(self.cell*np.array([a,b,c])[:,np.newaxis],axis=0)
                    new_xyz += (self.xyz+dispvect).tolist()
        self.xyz = np.array(new_xyz, "d")
        #print len(self.xyz)
        self.cell *= supercell[:,np.newaxis]
        self.cellparams[0:3] *= supercell
        self.natoms*=cellfact
        self.elems *=cellfact
        self.atypes*=cellfact
        self.images_cellvec = np.dot(images, self.cell)
        return
        
    def remove_duplicates(self, thresh=SMALL_DIST):
        badlist = []
        for i in xrange(self.natoms):
            for j in xrange(i+1, self.natoms):
                d,r,imgi=self.get_distvec(i,j)
                if d < thresh:
                    badlist.append(j)
        new_xyz = []
        new_elems = []
        new_atypes = []
        for i in xrange(self.natoms):
            if not badlist.count(i):
                new_xyz.append(self.xyz[i].tolist())
                new_elems.append(self.elems[i])
                new_atypes.append(self.atypes[i])
        self.xyz = np.array(new_xyz, "d")
        self.elems = new_elems
        self.natoms = len(self.elems)
        self.atypes = new_atypes
        return

    def wrap_in_box(self, thresh=SMALL_DIST):
        if not self.periodic: return
        # puts all atoms into the box again
        frac_xyz = self.get_frac_xyz()
        # now add 1 where the frac coord is negative and subtract where it is larger then 1
        frac_xyz += np.where(np.less(frac_xyz, 0.0), 1.0, 0.0)
        #frac_xyz -= np.where(np.greater_equal(frac_xyz, 1.0-thresh), 1.0, 0.0)
        #frac_xyz -= np.where(np.greater_equal(frac_xyz, 1.0), 1.0, 0.0)
        frac_xyz -= np.where(np.greater_equal(frac_xyz, 1.0+thresh*0.1), 1.0, 0.0)
        # convert back
        self.set_xyz_from_frac(frac_xyz)
        return
    
    def extend_cell(self,offset,bothsides=True):
        print 'warning, connectivity is destroyed!'
        frac_xyz = self.get_frac_xyz()
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
        self.set_xyz_from_frac(new_xyz)
        for i in range(len(wherexp[0])):
            self.elems.append(self.elems[wherexp[0][i]])
        for i in range(len(whereyp[0])):
            self.elems.append(self.elems[whereyp[0][i]])
        for i in range(len(wherezp[0])):
            self.elems.append(self.elems[wherezp[0][i]])
        for i in range(len(wherexm[0])):
            self.elems.append(self.elems[wherexm[0][i]])
        for i in range(len(whereym[0])):
            self.elems.append(self.elems[whereym[0][i]])
        for i in range(len(wherezm[0])):
            self.elems.append(self.elems[wherezm[0][i]])
        #print new_xyz
        self.natoms = len(self.xyz)
        

    def detect_conn(self, fixed_cutoff=None, exclude_pairs=None, cov_rad_buffer=0.1):
        self.conn = []
        for i in xrange(self.natoms):
            self.conn.append([])
        for i in xrange(self.natoms):
            for j in xrange(i+1,self.natoms):
                d,r,imgi=self.get_distvec(i,j)
                bond = False
                if fixed_cutoff:
                    if d<fixed_cutoff: bond = True
                else:
                    covradi = elements.cov_radii[self.elems[i]]
                    covradj = elements.cov_radii[self.elems[j]]
                    if d<(covradi+covradj+cov_rad_buffer) : bond = True
                # exclude pairs testing
                if exclude_pairs and bond:
                    el_p1,el_p2 = (self.elems[i], self.elems[j]),(self.elems[j], self.elems[i])
                    for expair in exclude_pairs:
                        if (expair == el_p1) or (expair == el_p2):
                            bond= False
                            break 
                if bond:
                    if len(imgi)>1:
                        raise ValueError, "Error in connectivity detection: use pconn!!!"
                    for ii in imgi:
                        self.conn[i].append(j)
                        self.conn[j].append(i)
        return
    
    def report_conn(self):
        print "REPORTING CONNECTIVITY"
        for i in xrange(self.natoms):
            conn = self.conn[i]
            print "atom %3d   %2s coordination number: %3d" % (i, self.elems[i], len(conn))
            for j in xrange(len(conn)):
                d = self.get_neighb_dist(i,j)
                print "   -> %3d %2s : dist %10.5f " % (conn[j], self.elems[conn[j]], d)
        return
    
    ###  specific to periodic systems .. cell manipulation ############
    
    def get_frac_xyz(self):
        if not self.periodic: return None
        cell_inv = np.linalg.inv(self.cell)
        return np.dot(self.xyz, cell_inv)
    
    def get_frac_from_real(self,real_xyz):
        if not self.periodic: return None
        cell_inv = np.linalg.inv(self.cell)
        return np.dot(real_xyz, cell_inv)
    
    def get_real_from_frac(self,frac_xyz):
        return np.dot(np.array(frac_xyz),self.cell)
        
    def set_xyz_from_frac(self, frac_xyz):
        if not self.periodic: return
        self.xyz = np.dot(frac_xyz,self.cell)

    def get_image(self,xyz, img):
        xyz = np.array(xyz)
        try:
            l = len(img)
            dispvec = np.sum(self.cell*np.array(img)[:,np.newaxis],axis=0)
        except TypeError:
            dispvec = np.sum(self.cell*np.array(images[img])[:,np.newaxis],axis=0)
        return xyz + dispvec
        
    def change_cell(self, new_cell):
        if not self.periodic: return
        frac_xyz = self.get_frac_xyz()
        if len(new_cell) == 6:
            self.cellparams = new_cell
            self.cell = unit_cell.vectors_from_abc(self.cellparams)
        elif len(new_cell) == 9:
            self.cell = new_cell
            self.cellparams = unit_cell.abc_from_vectors(self.cell)
        else:
            raise ValueError
        self.images_cellvec = np.dot(images, self.cell)
        self.set_xyz_from_frac(frac_xyz)
        return
        
    def scale_cell(self, scale):
        if not type(scale) == types.ListType:
            scale = 3*[scale]
        self.cellparams *= np.array(scale+[1,1,1])
        frac_xyz = self.get_frac_xyz()
        self.cell *= np.array(scale)[:,np.newaxis]
        self.images_cellvec = np.dot(images, self.cell)
        self.set_xyz_from_frac(frac_xyz)
        return
    
    ###  molecular manipulations #######################################
    
    def translate(self, vec):
        self.xyz += vec
        return
        
    def translate_frac(self, vec):
        if not self.periodic: return
        self.xyz += np.sum(self.cell*vec, axis=0)
        return
    
    def rotate_euler(self, euler):
        self.xyz = rotations.rotate_by_euler(self.xyz, euler)
        return
        
    def rotate_triple(self, triple):
        self.xyz = rotations.rotate_by_triple(self.xyz, triple)
        return

    def center_com(self):
        if self.periodic: return
        amass = []
        for e in self.elems: amass.append(elements.mass[e])
        amass = np.array(amass)
        center = np.sum((self.xyz*amass[:,np.newaxis]),axis=0)/np.sum(amass)
        self.translate(-center)
        return
    
    def get_com(self):
        amass = []
        for e in self.elems: amass.append(elements.mass[e])
        amass = np.array(amass,dtype='float64')
        #print amass, amass[:,np.newaxis].shape,self.xyz.shape
        center = np.sum((amass[:,np.newaxis]*self.xyz),axis=0)/np.sum(amass)
        #center = np.sum((amass[:,np.newaxis]*self.xyz),axis=0)/np.sum(amass)
        return center
    
    def get_coc(self,conns):
        amass = []
        for e in conns: amass.append(1.0)
        amass = np.array(amass,dtype='float64')
        #print amass, amass[:,np.newaxis].shape,self.xyz.shape
        conns=np.array(conns,dtype='int')
        center = np.sum((amass[:,np.newaxis]*self.xyz[conns,:]),axis=0)/np.sum(amass)
        #center = np.sum((amass[:,np.newaxis]*self.xyz),axis=0)/np.sum(amass)
        return center
    ###  system manipulations ##########################################
    
    def copy(self):
        return copy.deepcopy(self)
         
    def add_molsys(self, other, translate=None,rotate=None, scale=None, roteuler=None):
        # adds only nonperiodic molsys ... self can be both
        if other.periodic:
            if not (self.cell==other.cell).all():
                print "can not add periodic systems with unequal cells!!"
                return
        other_xyz = other.xyz.copy()
        # NOTE: it is important ot keep the order of operations
        #       1) scale
        #       2) rotate by euler angles
        #       3) rotate by orientation triple
        #       4) translate
        if scale    !=None:
            other_xyz *= np.array(scale)
        if roteuler != None:
            other_xyz = rotations.rotate_by_euler(other_xyz, roteuler)
        if rotate   !=None:
            other_xyz = rotations.rotate_by_triple(other_xyz, rotate,use_new=False)
        if translate!=None:
            other_xyz += translate
        if self.natoms==0:
            self.xyz = other_xyz
        else:
            self.xyz = np.array(self.xyz.tolist()+other_xyz.tolist(),"d")
        self.elems += other.elems
        self.atypes+= other.atypes
        for c in other.conn:
            cn = (np.array(c)+self.natoms).tolist()
            self.conn.append(cn)
        self.natoms += other.natoms
        return
        
    def make_periodic(self, periodic, cell):
        # makes system periodic without shifting, wrapping or changing connectivity!!
        self.periodic = periodic
        if type(cell) == types.ListType:
            print len(cell)
            if len(cell) == 6:
                self.cellparams = cell
                self.cell = unit_cell.vectors_from_abc(self.cellparams)
            elif len(cell) == 9:
                self.cell = np.array(cell).reshape([3,3])
                self.cellparams = unit_cell.abc_from_vectors(self.cell)
            else:
                raise ValueError
        else:
            if len(cell) == 6:
                # these are cellparams as numpy array
                self.cellparams = cell.tolist()
                self.cell = unit_cell.vectors_from_abc(self.cellparams)
            elif len(cell) == 3:
                # these are cellvectors as numpy array
                self.cell = cell
                self.cellparams = unit_cell.abc_from_vectors(self.cell)
        self.images_cellvec = np.transpose(np.dot(self.cell, np.transpose(images)))
        return

    def insert_atom(self, lab, aty, xyz, i, j):
        xyz.shape=(1,3)
        self.xyz = np.concatenate((self.xyz, xyz))
        self.elems.append(lab)
        self.atypes.append(aty)
        ci = self.conn[i]
        cj = self.conn[j]
        self.natoms += 1
        if ((i <= -1) or (j <= -1)): 
            self.conn.append([])
            return
        ci.remove(j)
        cj.remove(i)
        ci.append(self.natoms)
        cj.append(self.natoms)
        self.conn.append([i,j])
        return
        
    def add_bond(self, a1, a2):
        """ add a connection between a1 and a2 (in both directions)
        """
        if self.use_pconn:
            raise ValueError, "Can not add bonds to systems with pconn"
        self.conn[a1].append(a2)
        self.conn[a2].append(a1)
        return
        
    ### logical operations #####################################

    def is_superpose(self, other, thresh=1.0e-1):
        """ we test if two molecular systems are equal (superimpose)
        """
        if self.natoms != other.natoms: return False
        rmsd = 0.0
        for i in xrange(self.natoms):
            sxyz = self.xyz[i]
            r = other.xyz-sxyz
            d = np.sqrt(np.sum(r*r, axis=1))
            closest = np.argsort(d)[0]
            if d[closest] > thresh: return False, 0.0
            if self.elems[i] != other.elems[closest]: return False, 0.0
            rmsd += d[closest]
        rmsd = np.sqrt(np.sum(rmsd*rmsd))/self.natoms
        return True, rmsd
            
    
### additional stuff needed for the DEMOF version of weaver ####################################

    def find_molecules(self):
        self.mols = []
        self.moltypes = []
        # the default moleculename is "xyz" -> molecules from the xyz file
        self.molnames = ["xyz"]
        atoms = range(self.natoms)
        self.whichmol = self.natoms * [0]
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
                    new_leafs += self.conn[l]
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
            self.mols.append(curr_mol)
            for i in curr_mol: self.whichmol[i] = nmol
            # at this point all molecules found get the type 0 = "xyz"
            self.moltypes.append(0)
            nmol += 1
        # all atoms are assigned
        #if self.verbose: 
        #print "$$ -- found %d independent molecules from connectivity" % nmol
        self.nmols=nmol
        return
    
    def delete_atom(self,bad):
        new_xyz = [] 
        new_elems = []
        new_atypes = []
        new_conn = []
        for i in xrange(self.natoms):
            if i != bad:
                new_xyz.append(self.xyz[i].tolist())
                new_elems.append(self.elems[i])
                new_atypes.append(self.atypes[i])
                new_conn.append(self.conn[i])
                for j in xrange(len(new_conn[-1])):
                    if new_conn[-1].count(bad) != 0:
                        new_conn[-1].pop(new_conn[-1].index(bad))
        self.xyz = np.array(new_xyz, "d")
        self.elems = new_elems
        self.natoms = len(self.elems)
        self.atypes = new_atypes
        for i in range(len(new_conn)):
            #try: 
                #len(new_conn[i])
            #except:
                #new_conn[i] = [new_conn[i]]
            for j in range(len(new_conn[i])):
                if new_conn[i][j] >= bad:
                    new_conn[i][j]=new_conn[i][j]-1
        self.conn = new_conn
        return 
    

    def delete_atom_only(self,bad):
        new_xyz = [] 
        new_elems = []
        new_atypes = []
        new_conn = []
        for i in xrange(self.natoms):
            if i != bad:
                new_xyz.append(self.xyz[i].tolist())
                new_elems.append(self.elems[i])
                new_atypes.append(self.atypes[i])
                new_conn.append(self.conn[i])
                #for j in xrange(len(new_conn[-1])):
                    #if new_conn[-1].count(bad) != 0:
                        #new_conn[-1].pop(new_conn[-1].index(bad))
        self.xyz = np.array(new_xyz, "d")
        self.elems = new_elems
        self.natoms = len(self.elems)
        self.atypes = new_atypes
        for i in range(len(new_conn)):
            #try: 
                #len(new_conn[i])
            #except:
                #new_conn[i] = [new_conn[i]]
            for j in range(len(new_conn[i])):
                if new_conn[i][j] >= bad:
                    new_conn[i][j]=new_conn[i][j]-1
        self.conn = new_conn
        return 

    def delete_conn(self,el1,el2):  
        self.conn[el1].remove(el2)
        self.conn[el2].remove(el1)
        return
    
    def get_comdist(self,com,i):
        ri = self.xyz[i]
        rj = com
        if self.periodic:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - ri
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            d_sort = np.argsort(all_d)
            closest = d_sort[0]
            closest=[closest]  # THIS IS A BIT OF A HACK BUT WE MAKE IT ALWAYS A LIST ....
            if (abs(all_d[closest[0]]-all_d[d_sort[1]]) < SMALL_DIST):
                # oops ... there is more then one image atom in the same distance
                #  this means the distance is larger then half the cell width
                # in this case we have to return a list of distances
                for k in d_sort[1:]:
                    if (abs(all_d[d_sort[0]]-all_d[k]) < SMALL_DIST):
                        closest.append(k)
            d = all_d[closest[0]]
            r = all_r[closest[0]]
        else:
            r = rj-ri
            d = np.sqrt(np.sum(r*r))
            closest=[0]
        return d, r, closest
                
############# Plotting
            
    def plot(self,scell=False,bonds=False,labels=False):
        col = ['r','g','b','m','c','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k']+['k']*200
        fig = plt.figure(figsize=plt.figaspect(1.0)*1.5)
        ax = fig.add_subplot(111, projection='3d')
        atd = {}
        for i,aa in enumerate(list(set(self.atypes))):
            atd.update({aa:col[i]})
        print atd
        if bonds:
            for i in range(self.natoms):
                conn = self.conn[i]
                for j in range(len(conn)):
                    ax.plot([self.xyz[i][0],self.xyz[conn[j]][0]],[self.xyz[i][1],self.xyz[conn[j]][1]],[self.xyz[i][2],self.xyz[conn[j]][2]],color=atd[self.atypes[i]])
                        
        if labels:
            for i in range(self.natoms):
                label = str(i)+'-'+str(self.atypes[i]) +'-'+str(len(self.conn[i]))
                ax.text(self.xyz[i][0], self.xyz[i][1], self.xyz[i][2]+0.005, label, color='k',fontsize=9)
        if scell:
            xyz3 = self.make_333(out=True) 
            xyz3 =  np.array(xyz3)
            ax.scatter(xyz3[:,0],xyz3[:,1],xyz3[:,2],color='r',alpha=0.5)
        xyz=np.array(self.xyz)
        for i,xx in enumerate(xyz):
            ax.scatter(xx[0],xx[1],xx[2],color=atd[self.atypes[i]])
        minbound = np.min([np.min(xyz[:,0]),np.min(xyz[:,1]),np.min(xyz[:,2])])
        maxbound = np.max([np.max(xyz[:,0]),np.max(xyz[:,1]),np.max(xyz[:,2])])
        ax.auto_scale_xyz([0.0, maxbound], [0.0, maxbound], [0.0, maxbound])
        #ax.scatter(xyz1[:,0],xyz1[:,1],xyz1[:,2],color='k')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
            