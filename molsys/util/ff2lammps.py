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


mdyn2kcal = 143.88

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
        for i in xrange(self.mol.get_natoms()):
            vdwt = self.parind["vdw"][i][0]
            chrt = self.parind["cha"][i][0]
            vdwt_root = vdwt.split("->")[1]
            chrt_root = chrt.split("->")[1]
            if vdwt_root != chrt_root:
                raise ValueError, "The vdw and charge types of an atom do not match: %s %s" % (vdwt, chrt)
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
        header += "%10d atom types\n"       % len(self.par_types["vdw"])
        header += "%10d bond types\n"       % len(self.par_types["bnd"]) 
        header += "%10d angle types\n"      % len(self.par_types["ang"])
        header += "%10d dihedral types\n"   % len(self.par_types["dih"])
        header += "%10d improper types\n\n" % len(self.par_types["oop"])
        # currently we support onyl orthorhombic cells
        cell = self.mol.get_cell()
        header += '0.0 %12.6f  xlo xhi\n' % cell[0,0]
        header += '0.0 %12.6f  ylo yhi\n' % cell[1,1]
        header += '0.0 %12.6f  zlo zhi\n' % cell[2,2]
        header += '0.0 0.0 0.0  xy xz yz\n'        
        # NOTE in lammps masses are mapped on atomtypes which indicate vdw interactions (pair potentials)
        #   => we do NOT use the masses set up in the mol object because of this mapping
        #   so we need to extract the element from the vdw paramter name which is a bit clumsy
        header += "\nMasses\n\n"        
        self.masses = {}
        # this is a wild hack: we want to get the sorted list of the values of the dict, sorted by the entries (numbers)
        vdw_keys = [None]*self.npar["vdw"]
        for vdwk in self.par_types["vdw"]:
            vdw_keys[self.par_types["vdw"][vdwk]-1] = vdwk
        for vdwt in vdw_keys:
            etup = vdwt[0].split("->")[1].split("|")[0]
            etup = etup[1:-2]
            e = etup.split("_")[0]
            e = filter(lambda x: x.isalpha(), e)
            header += "%5d %10.4f # %s\n" % (self.par_types["vdw"][vdwt], elements.mass[e], vdwt)
        f.write(header)
        # write Atoms
        # NOTE ... this is MOF-FF and we silently assume that all charge params are Gaussians!!
        f.write("\nAtoms\n\n")
        xyz = self.mol.get_xyz()
        for i in xrange(self.mol.get_natoms()):
            vdwt  = self.parind["vdw"][i]
            atype = self.par_types["vdw"][tuple(vdwt)]
            molnumb = self.mol.molecules.whichmol[i]+1
            chrgpar    = self.par["cha"][self.parind["cha"][i][0]]
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
            oopt = tuple(self.parind["oop"][i])
            a,b,c,d  = self.rics["oop"][i]
            f.write("%10d %5d %8d %8d %8d %8d # %s\n" % (i+1, self.par_types["oop"][oopt], a+1, b+1, c+1, d+1, oopt))
        f.close()
        return

    def write_input(self, filename = "lmp.input", header=None, footer=None):
        """
        
        
        NOTE: add read data ... fix header with periodic info
        """
        self.input_filename = filename
        f = open(filename, "w")
        if header:
            hf = open(header, "r")
            f.write(hf.readlines())
            hf.close()
        f.write("\n# ------------------------ MOF-FF FORCE FIELD ------------------------------\n")
        # pair style
        f.write("\npair_style buck6d/coul/dsf %10.4f %10.4f\n\n" % (0.05, 12))
        for vdwpair in self.mol.ff.vdwdata:
            vdwi, vdwj = vdwpair.split(":")
            i = self.par_types["vdw"][(vdwi,)]
            j = self.par_types["vdw"][(vdwj,)]
            # in order to get the corresponding combined sigma values we need to do a bit of hacking here
            # in the init it was tested that the vdw and charge type names match, so we can use the root part for a lookup.
            chai = "gaussian->" + vdwi.split("->")[1]
            chaj = "gaussian->" + vdwj.split("->")[1]
            sigma_i = self.par["cha"][chai][1][1]
            sigma_j = self.par["cha"][chaj][1][1]
            sigma_ij = np.sqrt(sigma_i*sigma_i+sigma_j*sigma_j)
            r0, eps = self.mol.ff.vdwdata[vdwpair][1]
            f.write("pair_coeff %5d %5d %12.6f %12.6f %12.6f   # %s\n" % (i,j, eps, r0, 1.0/sigma_ij, vdwpair))            
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
                    raise ValueError, "not implemented"
                else:
                    raise ValueError, "unknown bond potential"
                f.write("bond_coeff %5d %s    # %s\n" % (bt_number, pstring, ibt))
        # angle style
        f.write("\nangle_style hybrid class2/mofff cosine/mofff\n\n")
        for at in self.par_types["ang"].keys():
            at_number = self.par_types["ang"][at]
            for iat in at:
                pot_type, params = self.par["ang"][iat]
                if pot_type == "mm3":
                    th0 = params[1]
                    K2  = params[0]*mdyn2kcal/2.0 
                    K3 = K2*(-0.14)
                    K4 = K2*5.6e-5
                    K5 = K2*-7.0e-7
                    K6 = K2*2.2e-8
                    pstring = "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f" % (th0, K2, K3, K4, K5, K6)
                    f.write("angle_coeff %5d class2/mofff    %s    # %s\n" % (at_number, pstring, iat))
                elif pot_type == "strbnd":
                    ksb1, ksb2, kss = params[:3]
                    r01, r02        = params[3:5]
                    th0             = params[5]
                    f.write("angle_coeff %5d class2/mofff bb %12.6f %12.6f %12.6f\n" % (at_number, kss*mdyn2kcal, r01, r02))
                    f.write("angle_coeff %5d class2/mofff ba %12.6f %12.6f %12.6f %12.6f\n" % (at_number, ksb1*mdyn2kcal, ksb2*mdyn2kcal, r01, r02))
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
