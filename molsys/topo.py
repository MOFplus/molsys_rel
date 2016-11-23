# -*- coding: utf-8 -*-
import string as st
import numpy as np
import types
import copy
import string
import logging

import util.elems as elements
import util.unit_cell as unit_cell
import util.rotations as rotations
import random
import molsys

try:
    from ase import Atoms
    from pyspglib import spglib
except ImportError:
    spg = False
else:
    spg = True

"""

        ToDo list ...

        - if use_pconn: write tinker xyz file with images (add a keyword in the file)
        ---> This is done, the file is now called topo file and can be i/o'd with read/write topo
        - if this keyword is found in a tinker file read in with pconn

"""

np.set_printoptions(threshold=20000)



deg2rad = np.pi/180.0
SMALL_DIST = 1.0e-3


class topo(molsys.molsys):
    def __init__(self):
        molsys.molsys.__init__(self)
        self.periodic= True
        self.use_pconn = True   # flag to use pconn: keeps the image number along with the bond
        self.pconn=[]
        # extra deafult sfor pyspglib
        if spg:
            self.symprec = SMALL_DIST        # precsion in symmetry detection .. pyspglib default of 1.0e-5 seems to be way too small for large systems (in Angstrom)
            self.nonhydrogen = False  # use only non-hydrogen atoms if True in symmetry detection or any operation
        return

    ####  I/O stuff ############################


    def write_topo(self,filename):
        ff=open(filename,'w')
        ff.write('%i %16.10f %16.10f %16.10f %10.6f %10.6f %10.6f\n' % (len(self.xyz),self.cellparams[0],self.cellparams[1],self.cellparams[2],self.cellparams[3],self.cellparams[4],self.cellparams[5]) )
        for i in xrange(self.natoms):
            line = ("%3d %-3s" + 3*"%12.6f" + " %5s") % \
               tuple([i+1]+[self.elems[i]]+ self.xyz[i].tolist() + [self.atypes[i]])
            conn = (np.array(self.conn[i])+1).tolist()
            pconn = self.pconn[i]
            pimg = []
            for pc in pconn:
                for ii,img in enumerate(images):
                    if all(img==pc): 
                        pimg.append(ii)
                        break
            if len(conn) != 0:
                for cc,pp in zip(conn,pimg):
                    if pp < 10:
                        line += "%8d/%1d" % (cc,pp)
                    else:
                        line += "%7d/%2d" % (cc,pp)
                #line += (len(conn)*"%7d") % tuple(conn)
            ff.write("%s \n" % line)
        ff.close()

    def write_v1(self,fname):
        if fname.split('.')[-1] != 'v1':
            fname = fname + '.v1'
        ff = open(fname,'w')
        ff.write(("%s\n") % ("unit cell vectors:"))
        ff.write("%s %12.6f %12.6f %12.6f\n" % ("va= ", self.cell[0,0], self.cell[1,0], self.cell[2,0]))
        ff.write("%s %12.6f %12.6f %12.6f\n" % ("va= ", self.cell[0,1], self.cell[1,1], self.cell[2,1]))
        ff.write("%s %12.6f %12.6f %12.6f\n" % ("va= ", self.cell[0,2], self.cell[1,2], self.cell[2,2]))
        ff.write("%i\n" % self.natoms)

        for i in range(self.natoms):
            s=self.elems[i]
            s = string.capitalize(s[0])+s[1:]
            ff.write("%s %12.6f %12.6f %12.6f\n" % (s,self.xyz[i][0], self.xyz[i][1], self.xyz[i][2]))
        ff.close()
        return
        
    ###### helper functions #######################

    def get_distvec2(self, i, j,exclude_self=True):
        """ vector from i to j
        This is a tricky bit, because it is needed also for distance detection in the blueprint
        where there can be small cell params wrt to the vertex distances.
        In other words: i can be bonded to j multiple times (each in a different image)
        and i and j could be the same!! """
        leni = True
        lenj = True
        try:
            l=len(i)
            if l > 1:
                ri = np.array(i)
            else:
                leni = False
                ri = self.xyz[i]
        except:
            ri = self.xyz[i]
        try:
            l=len(j)
            if l > 1:
                rj = np.array(j)
            else:
                rj = self.xyz[j]
        except:
            rj = self.xyz[j]
            lenj = False 
        if 1:
            all_rj = rj + self.images_cellvec
            all_r = all_rj - ri
            all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
            d_sort = np.argsort(all_d)
            if (np.linalg.norm(ri-rj) <= 0.001) and exclude_self:
                d_sort = d_sort[1:]
            closest = d_sort[0]
            closest=[closest]
            if (abs(all_d[closest[0]]-all_d[d_sort[1]]) < SMALL_DIST):
                for k in d_sort[1:]:
                    if (abs(all_d[d_sort[0]]-all_d[k]) < SMALL_DIST):
                        closest.append(k)
            d = all_d[closest[0]]
            r = all_r[closest[0]]
        return d, r, closest
    
    # thes following functions rely on an exisiting connectivity conn (and pconn)

    def get_neighb_coords(self, i, ci):
        """ returns coordinates of atom bonded to i which is ci'th in bond list """
        j = self.conn[i][ci]
        rj = self.xyz[j].copy()
        if self.periodic:
            if self.use_pconn:
                img = self.pconn[i][ci]
                rj += np.dot(img, self.cell)
            else:
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
            if self.use_pconn:
                img = self.pconn[i][ci]
                rj += np.dot(img, self.cell)
            else:
                all_rj = rj + self.images_cellvec
                all_r = all_rj - self.xyz[i]
                all_d = np.sqrt(np.add.reduce(all_r*all_r,1))
                closest = np.argsort(all_d)[0]
                return all_rj[closest]
        dr = ri-rj
        d = np.sqrt(np.sum(dr*dr))
        return d

    ######## manipulations in particular for blueprints

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
        # if pconn was set this needs to be redone now (conn is unaffected)
        if self.use_pconn:
            self.add_pconn()
        return
                
    
    def make_supercell_preserve_conn(self,supercell):
        img = [np.array(i) for i in images.tolist()]
        ntot = np.prod(supercell)
        nat = copy.deepcopy(self.natoms)
        nx,ny,nz = supercell[0],supercell[1],supercell[2]
        pconn = [copy.deepcopy(self.pconn) for i in range(ntot)]
        conn =  [copy.deepcopy(self.conn) for i in range(ntot)]
        xyz =   [copy.deepcopy(self.xyz) for i in range(ntot)]
        elems = copy.deepcopy(self.elems)
        left,right,front,back,bot,top =  [],[],[],[],[],[]
        neighs = [[] for i in range(6)]
        iii = []
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    ixyz = ix+nx*iy+nx*ny*iz
                    iii.append(ixyz)
                    if ix == 0   : left.append(ixyz)
                    if ix == nx-1: right.append(ixyz)
                    if iy == 0   : bot.append(ixyz)
                    if iy == ny-1: top.append(ixyz)
                    if iz == 0   : front.append(ixyz)
                    if iz == nz-1: back.append(ixyz)
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    ixyz = ix+nx*iy+nx*ny*iz
                    dispvect = np.sum(self.cell*np.array([ix,iy,iz])[:,np.newaxis],axis=0)
                    xyz[ixyz] += dispvect
                    
                    i = copy.copy(ixyz)
                    for cc in range(len(conn[i])):
                        for c in range(len(conn[i][cc])):
                            if (img[13] == pconn[i][cc][c]).all():
                                #conn[i][cc][c] += ixyz*nat
                                conn[i][cc][c] = int( conn[i][cc][c] + ixyz*nat )
                                pconn[i][cc][c] = np.array([0,0,0])
                            else:
                                px,py,pz     = pconn[i][cc][c][0],pconn[i][cc][c][1],pconn[i][cc][c][2]
                                #print px,py,pz
                                iix,iiy,iiz  = (ix+px)%nx, (iy+py)%ny, (iz+pz)%nz
                                iixyz= iix+nx*iiy+nx*ny*iiz
                                conn[i][cc][c] = int( conn[i][cc][c] + iixyz*nat )
                                pconn[i][cc][c] = np.array([0,0,0])
                                if ((px == -1) and (left.count(ixyz)  != 0)): pconn[i][cc][c][0] = -1
                                if ((px ==  1) and (right.count(ixyz) != 0)): pconn[i][cc][c][0] =  1   
                                if ((py == -1) and (bot.count(ixyz)   != 0)): pconn[i][cc][c][1] = -1
                                if ((py ==  1) and (top.count(ixyz)   != 0)): pconn[i][cc][c][1] =  1  
                                if ((pz == -1) and (front.count(ixyz) != 0)): pconn[i][cc][c][2] = -1
                                if ((pz ==  1) and (back.count(ixyz)  != 0)): pconn[i][cc][c][2] =  1   
                                #print px,py,pz
        self.conn, self.pconn, self.xyz = [],[],[]
        for cc in conn:
            for c in cc:
                self.conn.append(c)
        for pp in pconn:
            for p in pp:
                self.pconn.append(p)
        self.natoms = nat*ntot
        self.xyz = np.array(xyz).reshape(nat*ntot,3)
        self.cellparams[0:3] *= np.array(supercell)
        self.cell *= np.array(supercell)[:,np.newaxis]
        self.elems *= ntot
        self.atypes*=ntot
        self.images_cellvec = np.dot(images, self.cell)
        #print xyz
        return xyz,conn,pconn
    
    ######### connectivity things #################################

    def detect_conn(self, fixed_cutoff=None, pconn=False, exclude_pairs=None, cov_rad_buffer=0.1):
        self.conn = []
        if pconn: self.use_pconn=True
        for i in xrange(self.natoms):
            self.conn.append([])
            if self.use_pconn: self.pconn.append([])
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
                    if len(imgi)>1 and not self.use_pconn:
                        raise ValueError, "Error in connectivity detection: use pconn!!!"
                    for ii in imgi:
                        self.conn[i].append(j)
                        self.conn[j].append(i)
                        if self.use_pconn:
                            image = images[ii]
                            self.pconn[i].append(image)
                            self.pconn[j].append(image*-1)
        return
    
    def add_pconn(self):
        """ with the method detect_conn the connectivity is detected from a distance search
            if a connectivity is read via a tinker file there is no pconn present.
            with this metod it is added for the use with weaver2 """
        self.use_pconn= True
        self.pconn = []
        for i,c in enumerate(self.conn):
            atoms_pconn = []
            atoms_image = []
            for ji, j in enumerate(c):
                d,r,imgi = self.get_distvec(i,j)
                if len(imgi) == 1:
                    # only one neighbor .. all is fine
                    atoms_pconn.append(images[imgi[0]])
                    atoms_image.append(imgi[0])
                else:
                    # we need to assign an image to each connection
                    # if an atom is connected to another atom twice this means it must be another
                    # image
                    for ii in imgi:
                        # test if this image is not used for this atom .. then we can use it
                        if atoms_image.count(ii)==0:
                            atoms_image.append(ii)
                            atoms_pconn.append(images[ii])
                        else:
                            # ok, we have this image already
                            use_it = True
                            #print c
                            #print atoms_image
                            for k, iii in enumerate(atoms_image):
                                #print 'k',k
                                if (iii == ii) and (c[k] == j): use_it=False
                            if use_it:
                                atoms_image.append(ii)
                                atoms_pconn.append(images[ii])
            self.pconn.append(atoms_pconn)
        return

    def insert_atom(self, lab, aty, xyz, i, j):
        xyz.shape=(1,3)
        self.xyz = np.concatenate((self.xyz, xyz))
        self.elems.append(lab)
        self.atypes.append(aty)
        ci = self.conn[i]
        cj = self.conn[j]
        self.natoms += 1
        #print i, ci
        #print j, cj
        if ((i <= -1) or (j <= -1)): 
            self.conn.append([])
            #self.pconn.append([])
            return
        if self.use_pconn:
            pci = self.pconn[i].pop(ci.index(j))
            pcj = self.pconn[j].pop(cj.index(i))
        ci.remove(j)
        cj.remove(i)
        ci.append(self.natoms)
        cj.append(self.natoms)
        if self.use_pconn:
            self.pconn[i].append(np.zeros([3]))
            self.pconn[j].append(pcj)
        self.conn.append([i,j])
        if self.use_pconn:
            self.pconn.append([np.zeros([3]),pci])
        #print "end of insert .. conn:"
        #print self.conn
        return
        
    def add_bond(self, a1, a2):
        """ add a connection between a1 and a2 (in both directions)
        """
        if self.use_pconn:
            raise ValueError, "Can not add bonds to systems with pconn - well, we can fix this ;) "
        self.conn[a1].append(a2)
        self.conn[a2].append(a1)
        return
        

                
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
                    if self.pconn:
                        if np.sum(np.abs(self.pconn[i][j])) == 0:
                            ax.plot([self.xyz[i][0],self.xyz[conn[j]][0]],[self.xyz[i][1],self.xyz[conn[j]][1]],[self.xyz[i][2],self.xyz[conn[j]][2]],color='black')
                        else:
                            xyznew = self.get_image(self.xyz[conn[j]],self.pconn[i][j])
                            ax.scatter(xyznew[0],xyznew[1],xyznew[2],color='orange')
                            ax.plot([self.xyz[i][0],xyznew[0]],[self.xyz[i][1],xyznew[1]],[self.xyz[i][2],xyznew[2]],color='green')
                    else:
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

# ########## additional stuff for edge coloring ############################################

    def color_edges(self, proportions, maxiter=100, maxstep=100000, penref=0.3):
        """
        wrapper to search for a zero penalty solution of edge coloring. An initial coloring is generated
        randomly and a flipMC run is started ... if no zero penalty is found after maxsteps it is repeated ...
        """
        converged = False
        niter = 0
        while not (converged and (niter<maxiter)):
            self.init_color_edges(proportions)
            result = self.run_flip(maxstep, penref=penref)
            if result[0] : converged=True
            niter += 1
        print "**********************************************************************"
        print "edge coloring is convereged !!"
        print "final penalty is %12.6f" % result[2]
        return
                       

    def init_color_edges(self, proportions):
        """
        generate datastructures for edge coloring and set up random
        the proportions are a list or tuple of integer.
        the length determines the number of colors, the numbers are the multiples
        one integer should be always 1.
        so for example proportions=[2,1] generates a two color system (red, blue)
        with 2 red and 1 blue ... this means the total number of edges inthe periodic net must
        be a multiple of 3
        """
        self.prop = proportions
        self.ncolors = len(proportions)
        self.nprop = np.array(proportions).sum()
        # generate the list of bonds only "upwards" bonds (i1<i2) are stored
        blist = []
        for i, ci in enumerate(self.conn):
            for j in ci:
                if i<j: 
                    blist.append([i,j])
        # convert blist into numpy array
        self.blist = np.array(blist)
        self.nbonds = len(self.blist)
        assert self.nbonds%self.nprop==0,  "these proportions do not work"
        colors = []
        nc = self.nbonds/self.nprop
        for c, p in enumerate(self.prop):
            col = []
            for i in xrange(p*nc): col.append(c)
            colors += col
        self.colors = np.array(colors)
        numrand.shuffle(self.colors)
        # generate the bcolors table (same as self.conn but with colors)
        self.bcolors = copy.deepcopy(self.conn)
        for i in xrange(self.nbonds):
            self.set_bcol(i)
        # set up penalty table (per vertex)
        # defaults
        self.colpen_sum_fact    = 1.0
        self.colpen_orient_fact = 0.5
        # self.colpen_sumrule = self.natoms*[4.0]
        self.penalty = np.zeros([self.natoms],dtype="float64")
        # set up the scalmat array for all vertices
        self.setup_scalprod_mats()
        for i in xrange(self.natoms):
            self.penalty[i] = self.calc_colpen(i)
        self.totpen = self.penalty.sum()
        return
    
    def flip_color(self):
        """
        does one color flip and computes the new penalty
        the old settings are kept
        """
        # get bond A
        self.bondA = random.randint(0, self.nbonds-1)
        self.colA  = self.colors[self.bondA]
        vert = self.blist[self.bondA].tolist()
        peninit = [self.penalty[vert[0]], self.penalty[vert[1]]]
        # get bond B
        self.colB = self.colA
        while self.colA == self.colB:
            self.bondB = random.randint(0, self.nbonds-1)
            self.colB = self.colors[self.bondB]
        vert += self.blist[self.bondB].tolist()
        peninit += [self.penalty[vert[2]], self.penalty[vert[3]]]
        peninit = np.array(peninit)
        # now flip the colors
        self.colors[self.bondA] = self.colB
        self.colors[self.bondB] = self.colA
        # now correct entries in bcolors
        self.set_bcol(self.bondA)
        self.set_bcol(self.bondB)
        # compute penalty of he four affected vertices only
        pennew = []
        for v in vert:
            pennew.append(self.calc_colpen(v))
        self.pennew = np.array(pennew)
        self.changed_vert = vert
        # compute delta in penalty
        delta_pen = self.pennew.sum()-peninit.sum()
        # print "flipped colors of bonds %3d and %3d (vertices: %20s) -- delta penalty: %10.3f" % (self.bondA, self.bondB, str(vert), delta_pen)
        return delta_pen
    
    def unflip_colors(self):
        """
        call this directly after a flip to put everything back
        """
        self.colors[self.bondA] = self.colA
        self.colors[self.bondB] = self.colB
        self.set_bcol(self.bondA)
        self.set_bcol(self.bondB)
        return
    
    def accept_flip(self):
        """
        call this directly after flip to keep the flip
        """
        for i, v in enumerate(self.changed_vert):
            self.penalty[v] = self.pennew[i]
        self.totpen = self.penalty.sum()
        return
    
    def run_flip(self, maxstep, nprint=1000, penref=0.2, thresh=1.0e-3):
        """
        run a MC with color flips as moves for maxstep or until the total penalty
        is below thresh. the virtual "temperature" or reference energy for the acceptance
        (kT) is in the same unit as the penalty and is given as penref 
        
        :Parameters:
            - maxstep : maximum number of MC steps
            - nprint  : number of steps after which a printout is made [100]
            - penref  : reference penalty for the MC aceptance criterion exp(-pen/penref) [0.2]
            - thresh  : threshold under which convergence is assumed (zero penalty is not always reached for orientation penalty) [1.0e-3]
        """
        
        step = 0
        while (step < maxstep) and (self.totpen>thresh):
            dpen = self.flip_color()
            accept = False
            if dpen <= 0.0:
                accept = True
            else:
                prob = np.exp(-dpen/penref)
                #print "dpen %10.5f prob %10.5f" %(dpen, prob)
                if random.random() < prob:
                    accept=True
            if accept:
                self.accept_flip()
            else:
                self.unflip_colors()
            if (step%nprint) == 0:
                print "step: %7d ; penalty %10.5f" % (step, self.totpen)
            step += 1
        if step<maxstep:
            print "Converged after %7d steps with penalty %10.5f" % (step, self.totpen)
            print "last delta_pen was %10.5f" % dpen
            converged = True
        else:
            print "Not converged!!!"
            converged = False
        return (converged, step, self.totpen)
    
    def add_vertex_on_color(self, col, lelem, laty):
        for i,b in enumerate(self.blist):
            if self.colors[i] == col:
                # yes .. add a vertex here
                i,j = b
                ci = self.conn[i].index(j)
                xyz_j = self.get_neighb_coords(b[0], ci)
                xyz = (self.xyz[i]+xyz_j)/2.0
                self.insert_atom(lelem, laty, xyz, i, j)
        return
        
    # utility functions
    def set_bcol(self, bond):
        """
        utility to set color values in bcolors for bond
        
        :Parameters:
            - bond : index of bond i self.blist (and colors)
            
        """
        i, j = self.blist[bond]
        c = self.colors[bond]
        #print "bond from %3d to %3d : color %3d" % (i, j, c)
        #print "i_conn:  ", self.conn[i]
        #print "j_conn:  ", self.conn[j]
        ### set for i
        j_ind = self.conn[i].index(j)
        i_ind = self.conn[j].index(i)
        #print "i_ind ", i_ind
        #print "j_ind ", j_ind
        self.bcolors[i][j_ind] = c
        self.bcolors[j][i_ind] = c
        return
    
    def calc_colpen(self, vert):
        # print "calculating penalty for vert %d (%s)  colors: %s" % (vert, self.elems[vert], str(self.bcolors[vert]))
        pen_sum = self.calc_colpen_sum(vert)
        if pen_sum == 0.0:
            # this vertex has the correct number of colors on the edges
            # now compute in addition the penalty on the orientation
            return self.calc_colpen_orient(vert)
        else:
            return pen_sum
        
    def calc_colpen_sum(self, vert):
        """ compute the color penalty for vertex vert
            rules are in list self.colpen_sumrule
        """
        # get rule for first color (currently only two colors are supported for testing)
        nc0 = self.colpen_sumrule[vert]
        nc = self.bcolors[vert].count(0)
        pen = abs(nc-nc0)*self.colpen_sum_fact
        return pen*pen
    
    def set_colpen_sumrule(self, vert_dict):
        """
        set the color penalty for the sum of colors
        :Paramteres:
        
            - vert_dict: dictionary of vertices with the number of expected edges for color 0
        """
        self.colpen_sumrule = np.zeros([self.natoms], dtype="int32")
        for i in xrange(self.natoms):
            self.colpen_sumrule[i] = vert_dict[self.elems[i]]
        return
    
    def setup_scalprod_mats(self):
        self.scalmat = []
        for i in xrange(self.natoms):
            v = []
            for ji in xrange(len(self.conn[i])):
                v.append(self.get_neighb_coords(i, ji)-self.xyz[i])
            vn = vector.normalize(np.array(v))
            mat = np.sum(vn[:,np.newaxis,:]*vn[np.newaxis,:,:], axis=2)
            self.scalmat.append(mat)
        #self.scalmat = np.array(self.scalmat)
        return
        
    def calc_colpen_orient(self, vert):
        # this is a HACK ... works only for vertices with two colors
        # if self.colpen_orientrule == None ignore
        # if not use the number to be added to scal
        if self.colpen_orientrule[vert] != None:
            col0_edges = []
            for i,c in enumerate(self.bcolors[vert]):
                if c == 0: col0_edges.append(i)
            scal = self.scalmat[vert][col0_edges[0], col0_edges[1]]
            # print "check orient for vert %d  : %s %103f" % (vert, str(col0_edges), scal)
            return self.colpen_orient_fact*(scal+self.colpen_orientrule[vert])
        else:
            return 0.0
        
    def set_colpen_orientrule(self, vert_dict):
        """
        set the color penalty for the orientation of colors
        currently this works only for color zero sum=2 (if more ... how to dadd up penalties?)
        
        :Paramteres:
        
            - vert_dict: dictionary of vertices with either None or what to add to skal
            
        example: for color 0 (sum=2) being 180deg set it to 1.0 (-1.0+1.0 = 0.0)
                                            90deg set it to 0.0 and the fact to -0.5
        """
        self.colpen_orientrule = []
        for i in xrange(self.natoms):
            self.colpen_orientrule.append(vert_dict[self.elems[i]])
        return
         
                        
            