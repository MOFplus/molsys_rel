# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 09:37:25 2017

@author: rochus


              ff2lammps
              
class to be instantiated either with a mfpx file name or a mol object
it will setup the force field params from MOF+ and then write a data and a lamps input file              

"""

import numpy as np
import string
import molsys
import molsys.util.elems as elements
import copy


mdyn2kcal = 143.88
angleunit = 0.02191418
rad2deg = 180.0/np.pi 

class ff2lammps(object):
    
    def __init__(self, mol):
        """
        setup system and get parameter 
        
        :Parameters:
        
            - mol: either a string (mfpx file) or a mol object
        """
        
        if type(mol)==type(""):
            self.mol = molsys.mol()
            self.mol.read(mol)
        else:
            self.mol = mol
        # generate the force field
        self.mol.addon("ff")
        self.mol.ff.assign_params("MOF-FF")
        self.mol.ff.setup_pair_potentials()
        # set up the molecules
        self.mol.addon("molecules")
        self.mol.molecules()
        # make lists of paramtypes and conenct to mol.ff obejcts as shortcuts
        self.ricnames = ["bnd", "ang", "dih", "oop", "cha", "vdw"]
        self.par_types = {}
        self.par = {}
        self.parind = {}
        self.rics = {}
        self.npar = {}
        for r in self.ricnames:
            self.par[r]       = self.mol.ff.par[r]
            self.parind[r]    = self.mol.ff.parind[r]
            self.rics[r]      = self.mol.ff.ric_type[r]
            # sort identical parameters (sorted) using a tuple to hash it into the dict par_types : value is a number starting from 1 
            par_types = {}
            i = 1
            for pil in self.parind[r]:
                if pil:
                    pil.sort()
                    tpil = tuple(pil)
                    if not tpil in par_types:
                        par_types[tpil] = i
                        i += 1
            self.par_types[r] = par_types
            self.npar[r] = i-1
        # we need to verify that the vdw types and the charge types match because the sigma needs to be in the pair_coeff for lammps
        # thus we build our own atomtypes list combining vdw and cha and use the mol.ff.vdwdata as a source for the combined vdw params
        # but add the combined 1.0/sigma_ij here
        self.plmps_atypes = []
        self.plmps_pair_data = {}
        self.plmps_mass = {} # mass from the element .. even if the vdw and cha type differ it is still the same atom
        for i in xrange(self.mol.get_natoms()):
            vdwt = self.parind["vdw"][i][0]
            chrt = self.parind["cha"][i][0]
            at = vdwt+"/"+chrt
            if not at in self.plmps_atypes:
                #print("new atomtype %s" % at)
                self.plmps_atypes.append(at)
                # extract the mass ...
                etup = vdwt.split("->")[1].split("|")[0]
                etup = etup[1:-2]
                e = etup.split("_")[0]
                e = filter(lambda x: x.isalpha(), e)
                self.plmps_mass[at] = elements.mass[e]
                #print "with mass %12.6f" % elements.mass[e]
        for i, ati in enumerate(self.plmps_atypes):
            for j, atj in enumerate(self.plmps_atypes[i:],i):
                vdwi, chai = ati.split("/")
                vdwj, chaj = atj.split("/")
                vdwpairdata = self.mol.ff.vdwdata[vdwi+":"+vdwj]
                sigma_i = self.par["cha"][chai][1][1]
                sigma_j = self.par["cha"][chaj][1][1]
                # compute sigma_ij
                sigma_ij = np.sqrt(sigma_i*sigma_i+sigma_j*sigma_j)
                # vdwpairdata is (pot, [rad, eps])
                pair_data = copy.copy(vdwpairdata[1])
                pair_data.append(1.0/sigma_ij)
                self.plmps_pair_data[(i+1,j+1)] = pair_data
        # general settings                
        self._settings = {}
        # set defaults
        self._settings["cutoff"] = 12.0
        self._settings["parformat"] = "%15.8f"
        self._settings["vdw_a"] = 1.84e5
        self._settings["vdw_b"] = 12.0
        self._settings["vdw_c"] = 2.25
        self._settings["vdw_dampfact"] = 0.25
        return

    def setting(self, s, val):
        if not s in self._settings:
            print("This settings %s is not allowed" % s)
            return
        else:
            self._settings[s] = val
            return
        
    def write_data(self, filename="tmp.data"):
        self.data_filename = filename
        f = open(filename, "w")
        # write header 
        header = "LAMMPS data file for mol object with MOF-FF params from www.mofplus.org\n\n"
        header += "%10d atoms\n"      % self.mol.get_natoms()
        header += "%10d bonds\n"      % len(self.rics["bnd"])
        header += "%10d angles\n"     % len(self.rics["ang"])
        header += "%10d dihedrals\n"  % len(self.rics["dih"])
        header += "%10d impropers\n"  % len(self.rics["oop"])
        # types are different paramtere types 
        header += "%10d atom types\n"       % len(self.plmps_atypes)
        header += "%10d bond types\n"       % len(self.par_types["bnd"]) 
        header += "%10d angle types\n"      % len(self.par_types["ang"])
        header += "%10d dihedral types\n"   % len(self.par_types["dih"])
        header += "%10d improper types\n\n" % len(self.par_types["oop"])
        # currently we support onyl orthorhombic cells
        if self.mol.bcond == 0:
            # in the nonperiodic case center the molecule in the origin
            self.mol.translate(-self.mol.get_com())
            xyz = self.mol.get_xyz()
            cmax = xyz.max(axis=0)+10.0
            cmin = -xyz.min(axis=0)-10.0
        else:
            cell = self.mol.get_cell()
            cmin = np.zeros([3])
            cmax = cell.diagonal()
            xyz = self.mol.get_xyz()            
        header += '%12.6f %12.6f  xlo xhi\n' % (cmin[0], cmax[0])
        header += '%12.6f %12.6f  ylo yhi\n' % (cmin[1], cmax[1])
        header += '%12.6f %12.6f  zlo zhi\n' % (cmin[2], cmax[2])
        header += '0.0 0.0 0.0  xy xz yz\n'        
        # NOTE in lammps masses are mapped on atomtypes which indicate vdw interactions (pair potentials)
        #   => we do NOT use the masses set up in the mol object because of this mapping
        #   so we need to extract the element from the vdw paramter name which is a bit clumsy (DONE IN INIT NOW)
        header += "\nMasses\n\n"        
        for i in xrange(len(self.plmps_atypes)):
            at = self.plmps_atypes[i]
            header += "%5d %10.4f # %s\n" % (i+1, self.plmps_mass[at], at)
        f.write(header)
        # write Atoms
        # NOTE ... this is MOF-FF and we silently assume that all charge params are Gaussians!!
        f.write("\nAtoms\n\n")
        for i in xrange(self.mol.get_natoms()):
            vdwt  = self.parind["vdw"][i][0]
            chat  = self.parind["cha"][i][0]
            at = vdwt+"/"+chat
            atype = self.plmps_atypes.index(at)+1
            molnumb = self.mol.molecules.whichmol[i]+1
            chrgpar    = self.par["cha"][chat]
            assert chrgpar[0] == "gaussian", "Only Gaussian type charges supported"
            chrg = chrgpar[1][0]
            x,y,z = xyz[i]
            #   ind  atype molnumb chrg x y z # comment
            f.write("%10d %5d %5d %10.5f %12.6f %12.6f %12.6f # %s\n" % (i+1, molnumb, atype, chrg, x,y,z, vdwt))
        # write bonds
        f.write("\nBonds\n\n")
        for i in xrange(len(self.rics["bnd"])):
            bndt = tuple(self.parind["bnd"][i])
            a,b  = self.rics["bnd"][i]
            f.write("%10d %5d %8d %8d  # %s\n" % (i+1, self.par_types["bnd"][bndt], a+1, b+1, bndt))
        # write angles
        f.write("\nAngles\n\n")
        for i in xrange(len(self.rics["ang"])):
            angt = tuple(self.parind["ang"][i])
            a,b,c  = self.rics["ang"][i]
            f.write("%10d %5d %8d %8d %8d  # %s\n" % (i+1, self.par_types["ang"][angt], a+1, b+1, c+1, angt))
        # write dihedrals
        f.write("\nDihedrals\n\n")
        for i in xrange(len(self.rics["dih"])):
            diht = tuple(self.parind["dih"][i])
            a,b,c,d  = self.rics["dih"][i]
            f.write("%10d %5d %8d %8d %8d %8d # %s\n" % (i+1, self.par_types["dih"][diht], a+1, b+1, c+1, d+1, diht))
        # write impropers/oops
        f.write("\nImpropers\n\n")
        for i in xrange(len(self.rics["oop"])):            
            oopt = self.parind["oop"][i]
            if oopt:
                a,b,c,d  = self.rics["oop"][i]
                f.write("%10d %5d %8d %8d %8d %8d # %s\n" % (i+1, self.par_types["oop"][tuple(oopt)], a+1, b+1, c+1, d+1, oopt))
        f.write("\n")
        f.close()
        return

    def parf(self, n):
        pf = self._settings["parformat"]+" "
        return n*pf

    def write_input(self, filename = "lmp.input", header=None, footer=None):
        """
        NOTE: add read data ... fix header with periodic info
        """
        self.input_filename = filename
        f = open(filename, "w")
        # write standard header        
        f.write("clear\n")
        f.write("units real\n")
        if self.mol.bcond == 0:
            f.write("boundary f f f\n")
        else:
            f.write("boundary p p p\n")
        f.write("atom_style full\n")
        f.write("read_data %s\n\n" % self.data_filename)
        f.write("neighbor 2.0 bin\n\n")
        # extra header
        if header:
            hf = open(header, "r")
            f.write(hf.readlines())
            hf.close()
        f.write("\n# ------------------------ MOF-FF FORCE FIELD ------------------------------\n")
        # pair style
        f.write("\npair_style buck6d/coul/gauss/dsf %10.4f\n\n" % (self._settings["cutoff"]))
        for i, ati in enumerate(self.plmps_atypes):
            for j, atj in enumerate(self.plmps_atypes[i:],i):
                r0, eps, alpha_ij = self.plmps_pair_data[(i+1,j+1)]
                A = self._settings["vdw_a"]*eps
                B = self._settings["vdw_b"]/r0
                C = eps*self._settings["vdw_c"]*r0**6
                D = 6.0*(self._settings["vdw_dampfact"]*r0)**14
                f.write(("pair_coeff %5d %5d " + self.parf(5) + "   # %s <--> %s\n") % (i+1,j+1, A, B, C, D, alpha_ij, ati, atj))            
        # bond style
        f.write("\nbond_style hybrid class2 morse\n\n")
        for bt in self.par_types["bnd"].keys():
            bt_number = self.par_types["bnd"][bt]
            for ibt in bt:
                pot_type, params = self.par["bnd"][ibt]
                if pot_type == "mm3":
                    r0 = params[1]
                    K2 = params[0]*mdyn2kcal/2.0 
                    K3 = K2*(-2.55)
                    K4 = K2*(2.55**2.)*(7.0/12.0)
                    pstring = "class2 %12.6f %12.6f %12.6f %12.6f" % (r0, K2, K3, K4)
                elif pot_type == "morse":
                    r0 = params[1]
                    E0 = params[2]
                    k  = params[0]*mdyn2kcal/2.0
                    alpha = np.sqrt(k/E0)
                    pstring = "morse %12.6f%12.6f %12.6f" % (E0, alpha, r0)
                else:
                    raise ValueError, "unknown bond potential"
                f.write("bond_coeff %5d %s    # %s\n" % (bt_number, pstring, ibt))
        # angle style
        f.write("\nangle_style hybrid class2/p6 cosine/vdwl13\n\n")                
        # f.write("\nangle_style class2/mofff\n\n")
        for at in self.par_types["ang"].keys():
            at_number = self.par_types["ang"][at]
            for iat in at:
                pot_type, params = self.par["ang"][iat]
                if pot_type == "mm3":
                    th0 = params[1]
                    K2  = params[0]*mdyn2kcal/2.0 
                    K3 = K2*(-0.014)
                    K4 = K2*5.6e-5
                    K5 = K2*-7.0e-7
                    K6 = K2*2.2e-8
                    # pstring = "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f" % (th0, K2, K3, K4, K5, K6)
                    pstring = "%12.6f %12.6f" % (th0, K2)
                    f.write("angle_coeff %5d class2/p6    %s    # %s\n" % (at_number, pstring, iat))
                    # f.write("angle_coeff %5d    %s    # %s\n" % (at_number, pstring, iat))
                    # HACk to catch angles witout strbnd
                    if len(at) == 1:
                        f.write("angle_coeff %5d class2/p6 bb 0.0 1.0 1.0\n" % (at_number))
                        f.write("angle_coeff %5d class2/p6 ba 0.0 0.0 1.0 1.0\n" % (at_number))
                elif pot_type == "strbnd":
                    ksb1, ksb2, kss = params[:3]
                    r01, r02        = params[3:5]
                    th0             = params[5]
                    f.write("angle_coeff %5d class2/p6 bb %12.6f %12.6f %12.6f\n" % (at_number, kss*mdyn2kcal, r01, r02))
                    f.write("angle_coeff %5d class2/p6 ba %12.6f %12.6f %12.6f %12.6f\n" % (at_number, ksb1*mdyn2kcal, ksb2*mdyn2kcal, r01, r02))
                    # f.write("angle_coeff %5d bb %12.6f %12.6f %12.6f\n" % (at_number, kss*mdyn2kcal, r01, r02))
                    # f.write("angle_coeff %5d ba %12.6f %12.6f %12.6f %12.6f\n" % (at_number, ksb1*mdyn2kcal, ksb2*mdyn2kcal, r01, r02))
                elif pot_type == "fourier":
                    a0 = params[1]
                    fold = params[2]
                    k = 0.5*params[0]*angleunit*rad2deg*rad2deg/fold
                    pstring = "%12.6f %5d %12.6f" % (k, fold, a0)
                    f.write("angle_coeff %5d cosine/vdwl13   %s    # %s\n" % (at_number, pstring, iat))
                else:
                    raise ValueError, "unknown angle potential"
        # dihedral style
        f.write("\ndihedral_style opls\n\n")
        for dt in self.par_types["dih"].keys():
            dt_number = self.par_types["dih"][dt]
            for idt in dt:
                pot_type, params = self.par["dih"][idt]
                if pot_type == "cos3":
                    v1, v2, v3 = params[:3]
                    pstring = "%12.6f %12.6f %12.6f %12.6f" % (v1, v2, v3, 0.0)
                elif pot_type == "cos4":
                    v1, v2, v3, v4 = params[:4]
                    pstring = "%12.6f %12.6f %12.6f %12.6f" % (v1, v2, v3, v4)
                else:
                    raise ValueError, "unknown dihedral potential"
                f.write("dihedral_coeff %5d %s    # %s\n" % (dt_number, pstring, idt))
        # improper/oop style
        f.write("\nimproper_style inversion/harmonic\n\n")
        for it in self.par_types["oop"].keys():
            it_number = self.par_types["oop"][it]
            for iit in it:
                pot_type, params = self.par["oop"][iit]
                if pot_type == "harm":
                    pstring = "%12.6f %12.6f" % (params[0]*mdyn2kcal*1.5, params[1])
                else:
                    raise ValueError, "unknown improper/oop potential"
                f.write("improper_coeff %5d %s    # %s\n" % (it_number, pstring, iit))
        f.write("\nspecial_bonds lj 0.0 0.0 1.0 coul 1.0 1.0 1.0\n\n")
        f.write("# ------------------------ MOF-FF FORCE FIELD END --------------------------\n")
        # write footer
        if footer:
            ff = open(footer, "r")
            f.write(ff.readlines())
            ff.close()
        f.close()
        return
