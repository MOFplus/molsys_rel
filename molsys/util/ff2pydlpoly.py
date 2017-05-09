#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import numpy
import molsys
import string

attr_mappings = {
        "cnct": "conn",
        "types": "atypes",
        "boundarycond": "bcond",
        }

rad2deg = 180.0/numpy.pi 
class wrapper(object):

    def __init__(self,mol,is_master = True):
        self.is_master = is_master
        self.m = mol
        #
        self.virtual_atoms = []
        self.virtual_bonds = []
        self.virt_atoms_defs = []
        self.virtuals_atom_types = []
        #
        self.gaussians_used = True
        self.corechg_used = False
        self.spbasis_used = False
        #
        return

    def __getattr__(self,attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        if attr in attr_mappings.keys():
            return getattr(self.m, attr_mappings[attr])
        return getattr(self.m,attr)

    def read_tinker_xyz(self,fname):
        self.xyz_filename = fname
#        self.m.read(fname, ftype = "mfpx")
#        self.m.addon("ff")
#        self.m.ff.assign_params("MOF-FF")
        self.m.ff.setup_pair_potentials()
        self.m.addon("molecules")
        self.m.molecules()
        self.nmols = self.m.molecules.nmols
        self.molnames = self.m.molecules.molnames
        self.mols = self.m.molecules.mols
        self.moltypes = self.m.molecules.moltypes
        self.whichmol = self.m.molecules.whichmol
        return

    def empty_box(self,box):
        raise ValueError("Not implemted in wrapper")

    def get_system_from_pdlp(self,pdlpf):
        raise ValueError("Not implemented in wrapper")

    def read_pdlp_xyz(self, pdlpf, read_stage, velocities):
        raise ValueError("Not implemented in wrapper")

    def add_mol_tinker_xyz(self,file, mN, mname, offset, scale):
        raise ValueError("Not implemented in wrapper")

    def load_FF(self,keyfile):
        return

    def find_internals(self,do_smallring_check = False, do_lin_check = False):
        return

    def assign_FF(self,warning_FF):
        return

    def assign_molnames(self, *args, **kwargs):
        return

    def make_extra_mol(self, *args, **kwargs):
        raise ValueError("Not implemented in wrapper")

    def exclude_mol(self,*args, **kwargs):
        raise ValueError("Not implemented in wrapper")

    def exclude_rigid(self,*args, **kwargs):
        return

    def exclude_frozen(self,*args, **kwargs):
        return

    def report_molecules(self,*args, **kwargs):
        return

    def write2internal(self,pd):
        ric_type = {
                "bnd": self.m.ff.ric.bnd,
                "ang": self.m.ff.ric.ang,
                "dih": self.m.ff.ric.dih,
                "oop": self.m.ff.ric.oop}
        formatter = {"bnd": self.bondterm_formatter,
                "ang": self.angleterm_formatter,
                "dih": self.dihedralterm_formatter,
                "oop": self.oopterm_formatter}
        dest = {"bnd": pd.dlp_bnd.prmbnd,
                "ang": pd.dlp_ang.prmang,
                "dih": pd.dlp_dih.prmdih,
                "oop": pd.dlp_inv.prminv}
        ### bonded potentials ###
        for ict in ["bnd", "ang", "dih", "oop"]:
            counter = 0
            allprm = numpy.zeros(dest[ict].shape)
            for i, ic in enumerate(ric_type[ict]):
                ps = self.m.ff.parind[ict][i]
                for p in ps:
                    potential, params = self.m.ff.par[ict][p]
                    prm = formatter[ict](ic,potential,params,internal = True)
                    if type(prm) != type(None):
                        allprm[counter] = prm
                        counter += 1
            dest[ict] = allprm
            pd.set_atoms_moved()
        return
        

    def write_FIELD(self):
        ric_type = {
                "bnd": self.m.ff.ric.bnd,
                "ang": self.m.ff.ric.ang,
                "dih": self.m.ff.ric.dih,
                "oop": self.m.ff.ric.oop}
        formatter = {"bnd": self.bondterm_formatter,
                "ang": self.angleterm_formatter,
                "dih": self.dihedralterm_formatter,
                "oop": self.oopterm_formatter}
        header = {"bnd" :"BONDS",
                "ang": "ANGLES",
                "dih": "DIHEDRALS",
                "oop": "INVERSIONS"
                }
        f = open("FIELD", "w")
        f.write ("%s with autoassignment via MOFplus\n" % self.xyz_filename)
        f.write("UNITS kcal\n")
        f.write("\n")
        f.write("MOLECULES 1\n")
        f.write("%s\n" % self.xyz_filename)
        f.write("NUMMOLS 1\n")
        f.write("ATOMS %d\n" % self.natoms)
        ### setup and charges ###
        self.m.set_real_mass()
        chargesum = 0
        for i in xrange(self.natoms):
            p = self.m.ff.parind["vdw"][i][0]
            atype = self.m.ff.types2numbers[p]
            ### only gaussian implemented, no frozens
            p = self.m.ff.parind["cha"][i][0]
            potential, params = self.m.ff.par["cha"][p]
            if potential != "gaussian": 
                raise ValueError("Chargetype %s not implemented" % potential)
            chargesum += params[0]
            f.write("   %8s  %10.4f %10.4f %12.6f 1 %1d\n" % 
                    (atype, self.get_mass()[i], params[0], params[1], 0))
        print "Net charge of the system: %10.5f" % chargesum
        ### bonded potentials ###
        for ict in ["bnd", "ang", "dih", "oop"]:
            buffer_out = ""
            count = 0
            for i, ic in enumerate(ric_type[ict]):
                ps = self.m.ff.parind[ict][i]
                for p in ps:
                    potential, params = self.m.ff.par[ict][p]
                    term = formatter[ict](ic,potential,params)
                    if term:
                        buffer_out += term
                        count += 1
            f.write("%s %d\n" % (header[ict], count))
            f.write(buffer_out)
            print ict, count
        f.write("FINISH\n")
        ### pair potentials ###
        buffer_out = ""
        count = 0
        # RS: vdwdata contains the full matrix (both ways A:B and B:A) to work with codes
        #     that want them the other way around .. here we need to skip those which are done already
        done_pairs = []
        for pair in self.m.ff.vdwdata.keys():
            rev_pair = pair.split(":")
            rev_pair.reverse()
            if not rev_pair in done_pairs:
                potential, params = self.m.ff.vdwdata[pair]
                vdw_out = self.vdwterm_formatter(pair, potential, params)
                if vdw_out:
                    buffer_out += vdw_out
                    count += 1
                done_pairs.append(pair.split(":"))
        f.write("VDW %d\n" % count)
        f.write(buffer_out)
        f.write("CLOSE\n")
        f.close()
        return

    def bondterm_formatter(self, atoms, potential, params, unit = 71.94, iunit = 418.4, internal = False):
        assert len(atoms) == 2
        assert type(potential) == str
        assert type(params) == list
        if potential == "mm3":
            k  = 2.0*params[0]*unit
            r0 = params[1]
            if internal: return numpy.array([k*iunit, r0, 0.0, 0.0])
            return "   mm3b %5d %5d  %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, k, r0)
        elif potential == "morse":
            E0   = params[2]
            r0   = params[1]
            k    = params[0]
            alph = numpy.sqrt(unit*k/E0)
            if internal: return numpy.array([E0*iunit, r0, alph, 0.0])
            return "   mors %5d %5d  %10.5f %10.5f %10.5f\n" % (atoms[0]+1, atoms[1]+1, E0, r0, alph)
        elif potential == "quartic":
            k  = 2.0*params[0]*unit
            k2 = params[2]
            k3 = params[3]
            r0 = params[1]
            if internal: return numpy.array([k*iunit, r0, k2, k3])
            return "   aqua %5d %5d  %10.5f %10.5f %10.5f %10.5\n" % (atoms[0]+1, atoms[1]+1, k, r, k2, k3)
        else:
            raise IOError("Unknown bond potential %s" % potential)

    def angleterm_formatter(self, atoms, potential, params, unit = 0.02191418, 
            bunit = 71.94, strunit=2.51118, iunit = 418.4, internal = False):
        assert len(atoms) == 3
        assert type(potential) == str
        assert type(params) == list
        if potential == "mm3":
            k  = 2.0*params[0]*unit*rad2deg*rad2deg
            a0 = params[1]
            if internal: return numpy.array([iunit*k, a0*(1.0/rad2deg),0.0,0.0,0.0,0.0,0.0])
            return "   mm3a  %5d %5d %5d  %10.5f %10.5f\n" % (atoms[0]+1,
                    atoms[1]+1, atoms[2]+1, k, a0)
        elif potential == "quartic":
            k  = 2.0*params[0]*unit*rad2deg*rad2deg
            a0 = params[1]
            k2 = params[2]
            k3 = params[3]
            if internal: return numpy.array([iunit*k, a0*(1.0/rad2deg),k2,k3,0.0,0.0,0.0])
            return "   aqua  %5d %5d %5d  %10.5f %10.5f %10.5f %10.5f\n" % (atoms[0]+1,
                    atoms[1]+1, atoms[2]+1, k, a0, k2, k3)
        elif potential == "fourier":
            a0 = params[1]
            fold = params[2]
            k = 0.5*params[0]*unit*rad2deg*rad2deg/fold
            flag13 = " "
            vdw13s = params[3]
            chg13s = params[4]
            assert vdw13s == chg13s
            if vdw13s == 1.0: flag13 = "-"
            if internal: return numpy.array([k*iunit, a0*(1.0/rad2deg), fold, 0.0,0.0,0.0,0.0])
            return "  %1scos   %5d %5d %5d  %10.5f %10.5f %5d\n" % (flag13,atoms[0]+1,
                    atoms[1]+1, atoms[2]+1, k, a0, fold)
        elif potential == "strbnd":
            a  = 2.0*params[2]*bunit
            b  = params[0]*strunit*rad2deg
            c  = params[1]*strunit*rad2deg
            r1 = params[3]
            r2 = params[4]
            t  = params[5]
            if internal: return numpy.array([a*iunit, b*iunit, c*iunit, t*(1.0/rad2deg),r1,r2,0.0])
            return "   cmps  %5d %5d %5d  %10.5f %10.5f %10.5f % 10.5f %10.5f %10.5f\n" % (atoms[0]+1,
                    atoms[1]+1, atoms[2]+1, a,b,c,t,r1,r2)
        else:
            raise IOError("Unknown angle potential %s" % potential)

    def dihedralterm_formatter(self, atoms,potential, params, unit = 0.5, iunit =418.4,internal = False,
            settings = {"chg-14-scale":1.0, 
                "vdw-14-scale":1.0,
                "chg-13-scale": 1.0,
                "chg-12-scale": 1.0,
                }):
        assert len(atoms) == 4
        assert type(potential) == str
        assert type(params) == list
        if settings.has_key("chg-14-scale"):
            cou14=settings["chg-14-scale"]
        else:
            cou14=1.0
        if settings.has_key("vdw-14-scale"):
            vdw14=settings["vdw-14-scale"]
        else:
            vdw14=1.0
        if atoms.ring == 4:
            cou14 = 1.0
            vdw14 = 0.0
            if settings.has_key("chg-12-scale"): cou14 = settings["chg-12-scale"]
            if settings.has_key("vdw-12-scale"): cou14 = settings["vdw-12-scale"]
        elif atoms.ring == 5:
            cou14 = 1.0
            vdw14 = 0.0
            if settings.has_key("chg-13-scale"): cou14 = settings["chg-13-scale"]
            if settings.has_key("vdw-13-scale"): vdw14 = settings["vdw-13-scale"]
        if potential == "cos3":
            scale = unit*2.0
            V1 = params[0]*scale
            V2 = params[1]*scale
            V3 = params[2]*scale
            if ((V1==V2==V3==0.0) and (cou14==vdw14==1.0)):
                return None
            else:
                if internal: return numpy.array([V1*iunit, V2*iunit, V3*iunit,0.0,cou14,vdw14])
                return "   cos3  %5d %5d %5d %5d  %10.5f %10.5f %10.5f %10.5f  %10.5f  %10.5f\n" % \
                (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, V1, V2, V3, 0.0, cou14, vdw14)
        elif potential == "cos4":
            scale = unit*2.0
            V1 = params[0]*scale
            V2 = params[1]*scale
            V3 = params[2]*scale
            V4 = params[3]*scale
            if V1==V2==V3==V4==0.0:
                if cou14==vdw14==1.0: return None
            else:
                if internal: return numpy.array([V1*iunit, V2*iunit, V3*iunit,V4*iunit,cou14,vdw14])
                return "   cos4  %5d %5d %5d %5d  %10.5f %10.5f %10.5f %10.5f  %10.5f  %10.5f\n" % \
                (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, V1, V2, V3, V4, cou14, vdw14)
        else:
            raise IOError("Unknown dihedral potential %s" % potential)

    def oopterm_formatter(self, atoms, potential, params, unit=0.02191418,iunit=418.4, internal=False):
        if potential == "harm":
            k = 2.0*params[0]*3.0*unit*rad2deg*rad2deg
            a0 = params[1]*(1/rad2deg)
            if internal: return numpy.array([k*418.4, a0, 0.0, 0.0])
            return "   harm  %5d %5d %5d %5d   %10.5f  %10.5f\n" % \
                   (atoms[0]+1, atoms[1]+1, atoms[2]+1, atoms[3]+1, k, a0)
        else:
            raise IOError("Unknown oop potential %s" % potential)

    def vdwterm_formatter(self, pair, potential, params, 
            settings = {
                "a-expterm"  : 184000.0,
                "b-expterm"  : 12.0,
                "c-expterm"  : 2.25,
                "vdwdampfact": 0.25,
                }):
        t1,t2 = string.split(pair,":")
        dlp_t1 = self.m.ff.types2numbers[t1]
        dlp_t2 = self.m.ff.types2numbers[t2]
        if potential == "buck6d":
            A = settings["a-expterm"]*params[1]
            B = params[0]/settings["b-expterm"]
            C = params[1]*settings["c-expterm"]*params[0]**6
            D = 0.0
            Rcut = params[0]*settings["vdwdampfact"]
            return "   %8s %8s  ex6d %10.2f %10.5f %10.5f %10.5f\n" % (dlp_t1,dlp_t2, A, B, C, Rcut)
        else:
            raise IOError("Unknown vdw potential %s" % potential)


    def write_system_to_pdlp(self, pdlp):
        return

    def write_CONFIG(self,bcond=None, vel=False):
        f = open("CONFIG","w")
        f.write("some header\n")
        velflag = 0
        if vel: velflag = 1
        if self.cell is not None:
            if bcond: self.boundarycond = bcond
            f.write("%10d%10d%10d%20.8f\n" % (velflag, self.boundarycond, self.natoms, 0.0))
            f.write("%20.12f%20.12f%20.12f\n" % tuple(self.cell[0]))
            f.write("%20.12f%20.12f%20.12f\n" % tuple(self.cell[1]))
            f.write("%20.12f%20.12f%20.12f\n" % tuple(self.cell[2]))
        else:
            f.write("%10d%10d%10d%20.8f\n" % (velflag, 0, self.natoms, 0.0))
        for i in xrange(self.natoms):
            p = self.m.ff.parind["vdw"][i][0]
            atype = self.m.ff.types2numbers[p]
            f.write("%-8s%10d\n" % (atype, i+1))
            f.write("%20.12f%20.12f%20.12f\n" % tuple(self.xyz[i]))
            if vel:
                f.write("%20.12f%20.12f%20.12f\n" % tuple(self.vel[i]))
        f.close()
        return
        
 



