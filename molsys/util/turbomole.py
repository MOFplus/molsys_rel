"""
This module contains Turbomole related tools.
"""
import os
import numpy
import tempfile
import shutil
import molsys
from   molsys.util.units import angstrom

class DefineEndedAbnormallyError(Exception):
    def __init__(self, message=None, errors=None):
        # Calling the base class with the parameter it needs
        super().__init__(message)
        self.errors = errors
        return

class SymmetryAssignmentChangedtheStructureError(Exception):
    def __init__(self, message=None, errors=None):
        # Calling the base class with the parameter it needs
        super().__init__(message)
        self.errors = errors
        return

class dscfDidNotConverge(Exception):
    def __init__(self, message=None, errors=None):
        # Calling the base class with the parameter it needs
        super().__init__(message)
        self.errors = errors
        return



class GeneralTools:
    def __init__(self, path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        else:
            self.path = path
        return
        
    def invoke_define(self, define_in_name='define.in', coord_name='coord'):
        coord_path = os.path.join(self.path,coord_name)
        define_in_path = os.path.join(self.path,define_in_name)
        if not os.path.isfile(coord_path):
            raise FileNotFoundError("Please provide a coord file for invoking define.")
        if not os.path.isfile(define_in_path):
            raise FileNotFoundError("Please provide an input file for invoking define.")
        out_path = os.path.join(self.path, "define.out")
        err_path = os.path.join(self.path, "errormessage")
        os.system('define < %s > %s 2> %s' %(define_in_path, out_path, err_path))
        with open(err_path) as err:
            for line in err:
                if "define ended normally" not in line:
                    raise DefineEndedAbnormallyError(f"Define ended abnormally. Check the define input %s." % f_name)
                else:
                    os.remove(err_path)
        return

    def write_coord_from_mol(self, mol):
        coord_path = os.path.join(self.path,'coord')
        f=open(coord_path,'w')
        f.write('$coord\n')
        c = mol.xyz*angstrom
        for i in range(mol.natoms):
            f.write("  %19.14f %19.14f %19.14f   %-2s\n" %
                    (c[i,0],c[i,1], c[i,2], mol.elems[i]))
        f.write('$end')
        f.close()
        return   

    def coord_to_mol(self, coord_name="coord"):
        coord_path = os.path.join(self.path,coord_name)
        if not os.path.isfile(coord_path):
            raise FileNotFoundError("Please provide a coord file for invoking define.")
        coord_xyz_path = os.path.join(self.path,"coord.xyz")
        os.system("t2x %s > %s" % (coord_path, coord_xyz_path))
        mol = molsys.mol.from_file(coord_xyz_path, "xyz")
        os.remove(coord_xyz_path)
        return mol

    def read_charge_from_control(self):
        control_path = os.path.join(self.path,"control")
        if not os.path.isfile(control_path):
            raise FileNotFoundError("There is no control file in the directory %s." % self.path)
        with open(control_path) as control:
            for line in control:
                if line.startswith("$charge"):
                    charge = next(control).split()[0]
        return charge

    def get_nalpha_and_nbeta_from_ridft_output(self, ridft_out_name='ridft.out'):
        """ Read the number of alpha and beta electrons from the ridft output file. """
        ridft_path = os.path.join(self.path, ridft_out_name)
        if not os.path.isfile(ridft_path):
            raise FileNotFoundError("There is no ridft output file named %s in the directory %s." %(ridft_out_name, path))
        nalpha = 0
        nbeta = 0
        with open(ridft_path) as ridft_out:
            for line in ridft_out:
                ### get the number of alpha and beta shell occupations ###
                if line.startswith('   sum'):
                    nalpha = float(line.split()[1])
                    nbeta = float(line.split()[2])
        return nalpha, nbeta

    def calculate_spin_multiplicity_from(nalpha, nbeta):
        """ Calculate the spin multiplicity from the number of alpha and beta electrons:
            M = 2S+1 = 2*(nalpha-nbeta)*(1/2)+1 = (naplha-nbeta)+1
        """
        M = abs(round(nalpha-nbeta))+1
        return M

    
    def round_fractional_occupation(self, control_name='control'):
        """ Rounds the occupations and writes a new control file. """
        control_path = os.path.join(self.path,"control")
        if not os.path.isfile(control_path):
            raise FileNotFoundError("There is no control file in the directory %s." % self.path)
        THRESHOLD = 1.0E-7
        lines = []
        with open(control_path,'r') as control:
            for line in control:
                lines.append(line)
        new_lines = lines
        for i,line in enumerate(lines):
            if '$alpha shells' in line:
                j = 1
                while '$' not in lines[i+j]:
                    split_line = lines[i+j].strip().split()
                    occ = float(split_line[3])
                    if abs(occ-round(occ)) > THRESHOLD:
                        new_lines[i+j] = ' %s       %s                                     ( %s )\n' %(split_line[0],split_line[1],str(round(occ)))
                    j += 1
            if '$beta shells' in line:
                j = 1
                while '$' not in lines[i+j]:
                    split_line = lines[i+j].strip().split()
                    occ = float(split_line[3])
                    if abs(occ-round(occ)) > THRESHOLD:
                        new_lines[i+j] = ' %s       %s                                     ( %s )\n' %(split_line[0],split_line[1],str(round(occ)))
                    j += 1
        os.remove(control_path)
        with open(control_path,'w') as new_control:
            for line in new_lines:
                new_control.write(line)
        return

    def make_tmole_dft_input_for_fermi_smearing(self, elems, xyz, M_start, max_mem, title, lot_DFT):
        """Creates a tmole input called 'turbo.in'.

        Parameters
        ----------
        elems   : the list of elements
        xyz     : numpy array of shape (len(elems),3)
        M_start : the initial multiplicity before Fermi smearing
                   = 3 for systems with even number of electrons
                   = 2 for systems with odd number of electrons
        max_mem : Maximum memory per core to use in the calculations.
        title   : title of the job
        lot_DFT : DFT level of theory, must be string
        """
        turbo_in_path = os.path.join(self.path,"turbo.in")
        c = xyz*angstrom
        f = open(turbo_in_path,"w")
        f.write("%title\n")
        f.write("%s\n" % title)
        f.write("%method\n")
        f.write("ENRGY :: ri-u%s [gen_mult = %d, scf_dsta = 1, gen_symm = c1, scf_msil = 1000, scf_rico = %d, for_maxc = %d]\n" %
                (lot_DFT, M_start, max_mem, max_mem))
        f.write("%coord\n")
        for i in range(len(elems)):
            f.write("  %19.14f %19.14f %19.14f   %-2s\n" %
                    (c[i,0],c[i,1], c[i,2], elems[i]))
        f.write("%add_control_commands\n")
        f.write("$disp3\n")
        f.write("$fermi tmstrt=300.00 tmend= 50.00 tmfac=0.900 hlcrt=1.0E-01 stop=1.0E-03\n")
        f.write("ADD END\n")
        f.write("%end\n")
        f.close()
        return

    def run_tmole(self):
        maindir = os.getcwd()
        os.chdir(self.path)
        os.system("tmole")
        self.check_dscf_converged()
        os.chdir(maindir)
        return

    def check_dscf_converged(self):
        if os.path.isfile(os.path.join(self.path,"dscf_problem")):
            raise dscfDidNotConverge("SCF did not converge! Check the directory %s." % self.path)
        return

    def ridft(self):
        maindir = os.getcwd()
        os.chdir(self.path)
        os.system("ridft > ridft.out")
        os.chdir(maindir)
        return


    def kdg(self, dg_name=""):
        """Removes the data group named <dg_name>."""
        maindir = os.getcwd()
        os.chdir(self.path)
        os.system("kdg %s" % dg_name)
        os.chdir(maindir)
        return

class GeometryTools:

    def add_noise(mol, active_atoms = []):
        """ Returns xyz coordinates with a noise using a uniform distribution.

        First shifts the noise to the origin by subtracting 0.5,
        then divides it by 10 to get a noise up to 0.1 Angstrom.

        No noise is added to the active atoms.
        """
        if active_atoms != []:
            noise_ini = (numpy.random.rand(mol.natoms-len(active_atoms),3)-0.5)/10.0
            noise = GeometryTools._make_noise_without_active_atoms(noise_ini, active_atoms) 
        else:
            noise = (numpy.random.rand(mol.natoms,3)-0.5)/10.0
        # add the noise array to the coordinates array
        new_xyz = mol.xyz + noise
        return new_xyz

    def _make_noise_without_active_atoms(noise_ini, active_atoms):
        """ For the case of active atoms (in transition states), makes 
        sure that there is no noise added to the coordinates of these atoms.
        """
        noise_list = noise_ini.tolist()
        noise_to_add = []
        j = 0
        for i in range(mol.natoms):
            if (i+1) in active_atoms:
                noise_to_add.append([0.0,0.0,0.0])
            else:
                noise_to_add.append(noise_list[j])
                j += 1
        noise = numpy.array(noise_to_add)
        return noise

    def get_point_group_from_mol(mol):
        """
        invokes define in a temporary folder to get the assigned symmetry.
        """
        ### create a tmp directory ### 
        tmp_path = tempfile.mkdtemp()
        GeneralTools(tmp_path).write_coord_from_mol(mol)
        ### write a define file ###
        GeometryTools._write_define_check_symmetry(tmp_path)
        ### invoke define ###
        GeneralTools(tmp_path).invoke_define()
        ### read the assigned point group ###
        point_group = GeometryTools.get_point_group_from_control()
        ### remove the tmp directory ###
        shutil.rmtree(tmp_path)
        return point_group

    def _write_define_check_symmetry(path=os.getcwd()):
        """ Writes a 'define.in' file to get the detected point group. """
        define_in_path = os.path.join(path,'define.in')
        f=open(define_in_path,'w')
        f.write('\n\na coord\ndesy\n*\nno\nqq')
        f.close()
        return
    
    def get_point_group_from_control(path=os.getcwd()):
        """ Reads the point group from control file. """
        with open(os.path.join(path,'control')) as control:
             for lines in control:
                 if 'symmetry' in lines:
                     point_group = lines.strip().split()[1]
        return point_group


    def change_to_abelian_symmetry_group(path=os.getcwd()):
        """
        Changes the symmetry group to an abelian subgroup. 
        nel = number of singly occupied electrons
        """
        abelian_point_groups = { 'oh':'d2h',    'td':'c2v',    'th':'s6',      't':'d2',
                                'd3d':'s6',    'd2d':'c2v',   'd4h':'d2h',   'c4v':'c2v',
                                'c3v':'cs',    'c6v':'c2v',    'd4':'d2',     'd3':'c3'}
        point_group = GeometryTools.get_point_group_from_control(path)
        if point_group not in ['c1','cs']:
            ### determine the point group to be assigned ###
            try:
                new_point_group = abelian_point_groups[point_group]
            except:
                new_point_group = 'c1'
                print("The point group %s is not in the dictionary pointgroupdict." % point_group)
            print("It will try to assign the subgroup %s." % new_point_group)
            natoms = GeometryTools._get_natoms_from_control(path)
            ### write the define.in file ###
            GeometryTools._write_define_new_point_group(new_point_group, path)
            ### invoke define ###
            GeneralTools(path).invoke_define()
            ### check if the number of atoms is still the same ###
            newnatoms = GeometryTools._get_natoms_from_control(path)
            if natoms != newnatoms:
                raise SymmetryAssignmentChangedtheStructureError("The structure is far from the symmetry group %s. Therefore, while trying to change the symmetry group, new atoms are added to enforce it." % new_point_group)
                print('The number of atoms is not the same!')
        else: 
            print('The detected point group is already c1 or cs!')
        return


    def _write_define_new_point_group(new_point_group='c1', path=os.getcwd()):
        ### get the number of alpha and beta electrons ###
        nalpha, nbeta = GeneralTools(path).get_nalpha_and_nbeta_from_ridft_output()
        ### get the charge ###
        charge = GeneralTools(path).read_charge_from_control()
        ### write a define input ###
        define_in_path = os.path.join(path,"define.in")
        f = open(define_in_path,"w")
        f.write("\n")
        f.write("y\n")
        f.write("desy \n")
        f.write("sy %s\n" % new_point_group)
        f.write("ired\n")
        f.write("*\n")
        f.write("\n")
        f.write("eht\n")
        f.write("\n")
        f.write("%s\n" % charge)
        f.write("n\n")
        f.write("u %d\n" % abs(nalpha-nbeta))
        f.write("*\n\n\nq")
        f.close()
        return

    def _get_natoms_from_control(path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        with open(os.path.join(path,'control'),'r') as control:
            for lines in control:
                if 'natoms' in lines:
                    natoms = int(lines.strip().split('=')[-1])
        return natoms

    def assign_symmetry(path=os.getcwd()):
        changed = False
        ### get the number of alpha and beta electrons ###
        nalpha, nbeta = GeneralTools(path).get_nalpha_and_nbeta_from_ridft_output()
        ### get the charge ###
        charge = GeneralTools(path).read_charge_from_control()
        ### change the directory ###
        curdir = os.getcwd()
        os.chdir(path)
        os.system("cpc symmetry")
        tmp_path = os.path.join(path,"symmetry")
        define_in_file = os.path.join(tmp_path, "define.in")
        f = open(define_in_file,"w")
        f.write("\n")
        f.write("y\n")
        f.write("desy 1d-5\n")
        f.write("ired\n")
        f.write("*\n")
        f.write("\n")
        f.write("eht\n")
        f.write("\n")
        f.write("%s\n" % charge)
        f.write("n\n")
        f.write("u %d\n" % abs(nalpha-nbeta)) 
        f.write("*\n\n\nq")
        f.close()
        tmp_path = os.path.join(path,"symmetry")
        try:
            GeneralTools(tmp_path).invoke_define()
        except:
            GeometryTools.change_to_abelian_symmetry_group(tmp_path)
        point_group = GeometryTools.get_point_group_from_control(tmp_path)
        if point_group != "c1":
            changed = True
            for f in os.listdir(path):
                if os.isfile(f):
                    os.remove(f)
            for f in os.listdir(tmp_path):
                shutil.move(f,path)
        shutil.rmtree(tmp_path)
        return changed, point_group

        




class Mol:

    def __init__(self, mol):
        # The element symbols
        self.symbol = ['H', 'He', 'C', 'N', 'O', 'F', 'Ne', 'S', 'Cl', 'Ar']
        # The dictionary of number of electrons
        self.nelectrons = {'H':1,'He':2, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'S':16, 'Cl':17, 'Ar':18}
        # The dictionary of atomic m<F2>asses in a.u. 
        self.mass = {'H':1.0079,'He':4.0026, 'C':12.0107, 'N':14.0067, 'O':15.9994, 'F':18.9984, 'Ne':20.1797, 'S':32.065, 'Cl':35.453, 'Ar':39.948}
        # Molecular Mass
        self.MM = 0
        # Number of electrons
        self.nel = 0
        # The formula of the molecule. e.g. for methane CH4
        self.formula = ''
        self.types = []
        self.elems = []
        for el in mol.elems:
            self.elems.append(el.capitalize())
            if el.capitalize() not in self.types:
                self.types.append(el.capitalize())
        return

    def make_formula(self):
        """ Sorts the types alphabetically and counts for the number of each element type and write a formula. """
        for t in self.types:
            amount = self.elems.count(t)
            if amount != 1: 
                self.formula += t+str(amount)
            else: 
                self.formula += t
        return self.formula

    def count_number_of_electrons(self, charge=0):
        """ Counts the number of electrons. """
        ### count the number of electrons in the system ###
        for t in self.types:
           amount = self.elems.count(t)
           self.nel += self.nelectrons[t]*amount
        ### account for the charge of the molecule ###
        self.nel -= charge
        return self.nel

    def molecular_mass(self, mol):
        """ Calculates the molecular mass in a.u. """
        self.MM = 0
        for t in self.types:
           amount = self.elems.count(t)
           self.MM += self.mass[t]*amount
        return self.MM


class Submit:
    """Performs submission of the jobs."""
    def __init__(self, path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        else:
            self.path = path
        return

    def _freeze_atoms(self, active_atoms):
        """ writes a new define file to fix active atoms in internal coordinates. """
        define_in_path = os.path.join(self.path,'define.in')
        f = open(define_in_path, 'w')
        f.write(' \n')
        f.write('y\n')
        f.write('fix %s\n' %(active_atoms))
        f.write('ired\n')
        f.write('fix none\n')
        f.write('*\n')
        f.write('qq\n')
        f.close()
        GeneralTools(self.path).invoke_define()
        return

    def _ired_and_itvc_1(self):
        """writes a new define file to define internal redundant coordinates,
           changes the itvc to 1 (for TS optimization), 
           and changes the coordinates to the redundant internal coordinates.
        """
        define_in_path = os.path.join(self.path,'define.in')
        f = open(define_in_path, 'w')
        f.write(' \n')
        f.write('y\n')
        f.write('ired\n')
        f.write('*\n')
        f.write(' \n')
        f.write(' \n')
        f.write(' \n')
        f.write(' \n')
        f.write('stp\n')
        f.write('itvc\n')
        f.write('1\n')
        f.write('rmax 3e-2\n')
        f.write('*\n')
        f.write('q\n')
        f.close()
        GeneralTools(self.path).invoke_define()
        return

    def _check_geo_opt_converged(self):
        converged = False
        f_path = os.path.join(self.path,"GEO_OPT_CONVERGED")
        if os.path.isfile(f_path):
                converged = True
                os.remove(f_path)
        return converged

    def check_imaginary_frequency(self):
        f_vibspec = os.path.join(self.path,"vibspectrum")
        with open(f_vibspec) as vibspec:
           linenum = 0
           inum = 0 # number of imaginary frequencies
           imfreq = []
           for lines in vibspec:
               linenum += 1
               if (linenum > 3 and 'end' not in lines):
                   s = lines.rstrip('\n').split(None,3)
                   v = float(s[2])
                   if (v < 0):
                       inum += 1
                       imfreq.append(float(v))
        vibspec.close()
        #print('The number of imaginary frequencies is ', self.inum,'and they are/it is', self.imfreq)
        return inum, imfreq

    def jobex(self, ts=False):
        maindir = os.getcwd()
        os.chdir(self.path)
        if ts:
            os.system("jobex -ri -trans > jobex.out")
        else:
            os.system("jobex -ri > jobex.out")
        os.chdir(maindir)
        return

    def aoforce(self):
        maindir = os.getcwd()
        os.chdir(self.path)
        os.system("aoforce > aoforce.out")
        os.remove("dh")
        os.chdir(maindir)
        return


    def transition_state_workflow(self, active_atoms):
        found = False
        self._freeze_atoms(active_atoms)
        self.jobex()
        if not self._check_geo_opt_converged():
            print("Transition state pre-optimization did not converge.")
        return 


    def minima_workflow(self):
        found = False
        self.jobex()
        if not self._check_geo_opt_converged():
            print("Geometry optimization did not converge.")
        else:
            changed, point_group = GeometryTools.assign_symmetry(self.path)
            if changed:
                print("The point group is changed to %s." % point_group)
                self.jobex()
            self.aoforce()
            inum, imfreq = self.check_imaginary_frequency()
            if inum == 0:
                found = True
                print("The minima found succesfully.")
            else:
                print("There are imaginary frequencies. This is not a minima.")
        return found

class Harvest:
    
    def __init__(self, path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        else:
            self.path = path
        return

    def get_energy(self, f_ZPE = 1.0):
        aoforce_path = os.path.join(self.path, "aoforce.out")
        with open(aoforce_path) as aoforce:
            if '  zero point VIBRATIONAL energy  ' in line:
                self.ZPE = float(line.split()[6]) # The zero point vibrational energy in Hartree
            if 'SCF-energy' in line:
                self.SPE = float(line.split()[3]) # Energy in Hartree
        return (self.SPE + f_ZPE*self.ZPE)



