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

class ChargeIsNotWritteninControlFileError(Exception):
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

    def invoke_define(define_in_name='define.in', coord_name='coord', path=os.getcwd()):
        coord_path = os.path.join(path,coord_name)
        define_in_path = os.path.join(path,define_in_name)
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        if not os.path.isfile(coord_path):
            raise FileNotFoundError("Please provide a coord file for invoking define.")
        if not os.path.isfile(define_in_path):
            raise FileNotFoundError("Please provide an input file for invoking define.")
        os.system('define < define.in > define.out 2> errormessage')
        with open(os.path.join(path,'errormessage')) as err:
            for line in err:
                if "define ended normally" not in line:
                    raise DefineEndedAbnormallyError(f"Define ended abnormally. Check the define input %s." % f_name)
                else:
                    os.remove(os.path.join(path,"errormessage"))
        return

    def write_coord_from_mol(mol, path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        coord_path = os.path.join(path,'coord')
        f=open(coord_path,'w')
        f.write('$coord\n')
        c = mol.xyz*angstrom
        for i in range(mol.natoms):
            f.write("  %19.14f %19.14f %19.14f   %-2s\n" %
                    (c[i,0],c[i,1], c[i,2], mol.elems[i]))
        f.write('$end')
        f.close()
        return   

    def coord_to_mol(coord_name="coord", path=os.getcwd()):
        coord_path = os.path.join(path,coord_name)
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        if not os.path.isfile(coord_path):
            raise FileNotFoundError("Please provide a coord file for invoking define.")
        coord_xyz_path = os.path.join(path,"coord.xyz")
        os.system("t2x %s > %s" % (coord_path, coord_xyz_path))
        mol = molsys.mol.from_file(coord_xyz_path, "xyz")
        os.remove(coord_xyz_path)
        return mol

    def read_charge_from_control(path=os.getcwd()):
        control_path = os.path.join(path,"control")
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        if not os.path.isfile(control_path):
            raise FileNotFoundError("There is no control file in the directory %s." % path)
        with open(control_path) as control:
            for line in control:
                if line.startswith("$charge"):
                    charge = next(control).split()[0]
                else:
                    raise ChargeIsNotWritteninControlFileError("Please make sure that $charge is written in the control file.")
        return charge

    def get_nalpha_and_nbeta_from_ridft_output(ridft_out_name='ridft.out', path=os.getcwd()):
        """ Read the number of alpha and beta electrons from the ridft output file. """
        ridft_path = os.path.join(path, ridft_out_name)
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
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

    
    def round_fractional_occupation(control_name='control', path=os.getcwd()):
        """ Rounds the occupations and writes a new control file. """
        control_path = os.path.join(path,"control")
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        if not os.path.isfile(control_path):
            raise FileNotFoundError("There is no control file in the directory %s." % path)
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
                    print(lines[i+j])
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

    def make_tmole_dft_input_for_fermi_smearing(elems, xyz, M_start, path, max_mem, title, lot_DFT):
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
        turbo_in_path = os.path.join(path,"turbo.in")
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

    def run_tmole(path=os.getcwd()):
        maindir = os.getcwd()
        os.chdir(path)
        os.system("tmole")
        if os.path.isfile("dscf_problem"):
            raise dscfDidNotConverge("dscf did not converge! Check the directory %s." % path)
        os.chdir(maindir)
        return

    def kdg(path=os.getcwd(), dg_name=""):
        """Removes the data group named <dg_name>."""
        maindir = os.getcwd()
        os.chdir(path)
        os.system("kdg $%s" % dg_name)
        os.chdir(maindir)
        return





class GeometryTools:

    def __init__(self):
        self.abelian_point_groups = { 'oh':'d2h',    'td':'c2v',    'th':'s6',      't':'d2',
                                     'd3d':'s6',    'd2d':'c2v',   'd4h':'d2h',   'c4v':'c2v',
                                     'c3v':'cs',    'c6v':'c2v',    'd4':'d2',     'd3':'c3'}
        self.symbol = {1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 16: 'S', 17: 'Cl', 18: 'Ar'}
        return

    def add_noise(self, mol, active_atoms = []):
        """ Returns xyz coordinates with a noise using a uniform distribution.

        First shifts the noise to the origin by subtracting 0.5,
        then divides it by 10 to get a noise up to 0.1 Angstrom.

        No noise is added to the active atoms.
        """
        if active_atoms != []:
            noise_ini = (numpy.random.rand(mol.natoms-len(active_atoms),3)-0.5)/10.0
            noise = self._make_noise_without_active_atoms(noise_ini, active_atoms) 
        else:
            noise = (numpy.random.rand(mol.natoms,3)-0.5)/10.0
        # add the noise array to the coordinates array
        new_xyz = mol.xyz + noise
        return new_xyz

    def _make_noise_without_active_atoms(self, noise_ini, active_atoms):
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

    def get_point_group_from_mol(self, mol):
        """
        invokes define in a temporary folder to get the assigned symmetry.
        """
        ### create a tmp directory ### 
        tmp_path = tempfile.mkdtemp()
        GeneralTools.write_coord_from_mol(mol=mol, path=tmp_path)
        ### write a define file ###
        self._write_define_check_symmetry(path=tmp_path)
        ### invoke define ###
        GeneralTools.invoke_define(path=tmp_path)
        ### read the assigned point group ###
        point_group = self.get_assigned_point_group()
        ### remove the tmp directory ###
        shutil.rmtree(tmp_path)
        return point_group

    def _write_define_check_symmetry(self, path=os.getcwd()):
        """ Writes a 'define.in' file to get the detected point group. """
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        define_in_path = os.path.join(path,'define_in')
        f=open(define_in_path,'w')
        f.write('\n\na coord\ndesy\n*\nno\nqq')
        f.close()
        return
    
    def get_point_group_from_control(self, path=os.getcwd()):
        """ Reads the point group from control file. """
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        with open(os.path.join(path,'control')) as control:
             for lines in control:
                 if 'symmetry' in lines:
                     point_group = lines.strip().split()[1]
        return point_group


    def change_to_abelian_symmetry_group(self, path=os.getcwd()):
        """
        Changes the symmetry group to an abelian subgroup. 
        nel = number of singly occupied electrons
        """
        point_group = self.get_point_group_from_control(path)
        if point_group not in ['c1','cs']:
            ### determine the point group to be assigned ###
            try:
                new_point_group = self.abelian_point_groups[point_group]
            except:
                new_point_group = 'c1'
                print("The point group %s is not in the dictionary pointgroupdict." % point_group)
            print("It will try to assign the subgroup %s." % new_point_group)
            natoms = _get_natoms_from_control(path)
            ### write the define.in file ###
            self._write_define_new_point_group(new_point_group=new_point_group, path=path)
            ### invoke define ###
            GeneralTools.invoke_define(path=path)
            ### check if the number of atoms is still the same ###
            newnatoms = _get_natoms_from_control(path)
            if natoms != newnatoms:
                raise SymmetryAssignmentChangedtheStructureError("The structure is far from the symmetry group %s. Therefore, while trying to change the symmetry group, new atoms are added to enforce it." % new_point_group)
                print('The number of atoms is not the same!')
        else: 
            print('The detected point group is already c1 or cs!')
        return


    def _write_define_new_point_group(self, new_point_group='c1', path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        ### get the number of alpha and beta electrons ###
        nalpha, nbeta = GeneralTools.get_nalpha_and_nbeta_from_ridft_output(path=path)
        ### get the charge ###
        charge = GeneralTools.read_charge_from_control(path=path)
        ### write a define input ###
        f = open(os.path.join(path,'define.in'),'w')
        f.write('\ny\ndesy \nsy %s\nired\n*\n\neht\n\n%s\nn\nu %d\n*\n\n\nq' %(new_point_group,charge,abs(nalpha-nbeta)))
        f.close()
        return

    def _get_natoms_from_control(self, path):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        with open(os.path.join(path,'control'),'r') as control:
            for lines in control:
                if 'natoms' in lines:
                    natoms = int(lines.strip().split('=')[-1])
        return natoms





class Mol:

    def __init__(self, mol):
        # The element symbols
        self.symbol = ['H', 'He', 'C', 'N', 'O', 'F', 'Ne', 'S', 'Cl', 'Ar']
        # The dictionary of number of electrons
        self.nelectrons = {'H':1,'He':2, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'S':16, 'Cl':17, 'Ar':18}
        # The dictionary of atomic masses in a.u. 
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
