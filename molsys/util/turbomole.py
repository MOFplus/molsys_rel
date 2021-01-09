"""
This module contains Turbomole related tools.
"""
import sys
import os
import numpy
import tempfile
import time
import shutil
import molsys
from   molsys.util.units import angstrom
import matplotlib.pyplot as plt
from molsys.addon import zmat
import graph_tool
from graph_tool import Graph, GraphView

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


class GeneralTools:
    def __init__(self, path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        else:
            self.path    = os.path.abspath(path)
            self.maindir = os.getcwd()
        return
        
    def invoke_define(self, define_in_name='define.in', define_out_name = 'define.out', coord_name='coord'):
        os.chdir(self.path)
        coord_path = os.path.join(self.path,coord_name)
        define_in_path = os.path.join(self.path,define_in_name)
        if not os.path.isfile(coord_path):
            raise FileNotFoundError("Please provide a coord file for invoking define.")
        if not os.path.isfile(define_in_path):
            raise FileNotFoundError("Please provide an input file for invoking define.")
        out_path = os.path.join(self.path, define_out_name)
        err_path = os.path.join(self.path, "errormessage")

        # 1. invoke define
        os.system('define < %s > %s 2> %s' %(define_in_path, out_path, err_path))

        # 2. check if define was succesful
        with open(err_path) as err:
            for line in err:
                if "define ended normally" not in line:
                    raise DefineEndedAbnormallyError("Define ended abnormally. Check the define input under %s." % define_in_path)
                else:
                    os.remove(err_path)
        os.chdir(self.maindir)
        return

    def write_coord_from_mol(self, mol, coord_name = 'coord'):
        coord_path = os.path.join(self.path,coord_name)
        f=open(coord_path,'w')
        f.write('$coord\n')
        c = mol.xyz*angstrom
        for i in range(mol.natoms):
            f.write("  %19.14f %19.14f %19.14f   %-2s\n" %
                    (c[i,0],c[i,1], c[i,2], mol.elems[i]))
        f.write('$end')
        f.close()
        return   

    def mol_to_coord(self, mol):
        coord_str = '$coord\n'
        c = mol.xyz*angstrom
        for i in range(mol.natoms):
            coord_str += '  %19.14f %19.14f %19.14f   %-2s\n' %(c[i,0],c[i,1], c[i,2], mol.elems[i])
        coord_str += '$end'
        return coord_str

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
        c = None
        control_path = os.path.join(self.path,"control")
        if not os.path.isfile(control_path):
            raise FileNotFoundError("There is no control file in the directory %s." % self.path)
        with open(control_path) as control:
            for line in control:
                if line.startswith("$charge"):
                    charge = float(next(control).split()[0])
        c = round(charge)
        return c

    def read_ssquare_from_control(self):
        ssquare = None
        control_path = os.path.join(self.path,"control")
        if not os.path.isfile(control_path):
            raise FileNotFoundError("There is no control file in the directory %s." % self.path)
        with open(control_path) as control:
            for line in control:
                if line.startswith("$ssquare"):
                    ssquare = float(next(control).split()[0])
        return ssquare
  

    def get_nalpha_and_nbeta_from_ridft_output(self, ridft_out_name='ridft.out'):
        """ Read the number of alpha and beta electrons from the ridft output file. """
        ridft_path = os.path.join(self.path, ridft_out_name)
        if not os.path.isfile(ridft_path):
            raise FileNotFoundError("There is no ridft output file named %s in the directory %s." %(ridft_out_name, self.path))
        nalpha = 0
        nbeta = 0
        with open(ridft_path) as ridft_out:
            for line in ridft_out:
                # Get the number of alpha and beta shell occupations
                if line.startswith('   sum'):
                    nalpha = float(line.split()[1])
                    nbeta = float(line.split()[2])
        return nalpha, nbeta


    def calculate_spin_multiplicity_from(self, nalpha, nbeta):
        """ Calculate the spin multiplicity from the number of alpha and beta electrons:
            M = 2S+1 = 2*(nalpha-nbeta)*(1/2)+1 = (nalpha-nbeta)+1
        """
        M = abs(nalpha-nbeta)+1
        return M

    
    def round_fractional_occupation(self, THRESHOLD = 1.0E-7):
        """ Rounds the occupations and writes a new control file. """
        control_path = os.path.join(self.path,'control')
        if not os.path.isfile(control_path):
            raise FileNotFoundError("There is no control file in the directory %s." % self.path)

        # 1. Read the lines of the control file
        lines = []
        with open(control_path,'r') as control:
            for line in control:
                lines.append(line)
        new_lines = lines

        # 2. If the difference between the occupation and the rounded occupation is less than the threshold, round the occupation number
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

        # 3. Rewrite the control file with the round occupation numbers
        os.remove(control_path)
        with open(control_path,'w') as new_control:
            for line in new_lines:
                new_control.write(line)
        return


    def for_c1_sym_change_multiplicity_in_control_by(self, N, nalpha, nbeta):
        ''' Changes the multiplicity iby N, by changing the occupation of the alpha and beta electrons in the control file. 
        '''
        if N % 2 != 0:
            print('Provide an even number.')
        else:
            control_path = os.path.join(self.path,'control')
            if not os.path.isfile(control_path):
                raise FileNotFoundError("There is no control file in the directory %s." % self.path)

            # 1. Get the line number where alpha shells are written
            with open(control_path,'r') as control:
                for i,line in enumerate(control):
                    if '$alpha shells' in line:
                       line_alpha = i
            
            # 2. Determine the new occupations
            if nalpha >= nbeta:
                new_alpha = nalpha+N/2
                new_beta  = nbeta-N/2
            elif nalpha < nbeta:
                new_alpha = nalpha-N/2
                new_beta  = nbeta+N/2
            
            # 3. Remove the old occupations from the control file
            self.kdg('alpha shells')
            self.kdg('beta shells')

            # 4. Add the new occupations to the control file
            newlines = ''
            with open(control_path) as control:
               for i,line in enumerate(control):
                   if i == line_alpha:
                       newlines += '$alpha shells\n'
                       newlines += ' a       1-%d   ( 1 )\n' %(new_alpha)
                       newlines += '$beta shells\n' 
                       newlines += ' a       1-%d   ( 1 )\n' %(new_beta)
                   newlines += line
            f = open(control_path,'w')
            f.write(newlines)
            f.close()
        return


    def make_tmole_dft_input(self, elems, xyz, M, max_mem, title, lot, genprep = 0, scf_dsta = 1.0, fermi = True, nue = False):
        """Creates a tmole input called 'turbo.in' with c1 symmetry.
        Parameters
        ----------
        elems   : the list of elements
        xyz     : numpy array of shape (len(elems),3)
        M       : the (initial) spin multiplicity
        max_mem : Maximum memory per core to use in the calculations.
        title   : title of the job
        scf_dsta: start value for the SCF damping.
        lot     : The QM level of theory, must be string, and according to Tmole manual
        genprep : 0 -> perform calculation, 1 -> only prepare the input
        fermi   : Boolean for Fermi smearing
        """
        turbo_in_path = os.path.join(self.path,"turbo.in")
        c = xyz*angstrom
        f = open(turbo_in_path,"w")
        f.write("%title\n")
        f.write("%s\n" % title)
        f.write("%method\n")
        f.write("ENRGY :: %s [gen_mult = %d, gen_symm = c1, scf_dsta = %f, scf_msil = 1000, scf_rico = %d, for_maxc = %d, genprep = %d]\n" %
                (lot, M, scf_dsta, 0.3*max_mem, 0.7*max_mem, genprep))
        f.write("%coord\n")
        for i in range(len(elems)):
            f.write("  %19.14f %19.14f %19.14f   %-2s\n" %
                    (c[i,0],c[i,1], c[i,2], elems[i]))
        f.write("%add_control_commands\n")
        f.write("$disp3\n")
        if fermi:
            if nue:
                f.write("$fermi tmstrt=300.00 tmend= 50.00 tmfac=0.900 hlcrt=1.0E-01 stop=1.0E-03 nue=%d\n" %M)
            else:
                f.write("$fermi tmstrt=300.00 tmend= 50.00 tmfac=0.900 hlcrt=1.0E-01 stop=1.0E-03\n")
        f.write("ADD END\n")
        f.write("%end\n")
        f.close()
        return

    def run_tmole(self):
        os.chdir(self.path)
        os.system("tmole &>/dev/null")
        os.chdir(self.maindir)
        return

    def check_ridft_converged(self, ridft_out_name = 'ridft.out'):
        converged = True
        ridftout = open(os.path.join(self.path, ridft_out_name),'r')
        for line in ridftout:
            if 'ATTENTION: ridft did not converge!' in line:
                converged = False
        return converged

    def get_energy_from_ridft_out(self, ridft_out_name = 'ridft.out'):
        ridftout = open(os.path.join(self.path, ridft_out_name),'r')
        for line in ridftout:
           l = line.split()
           if  'total' in l and 'energy' in l and '=' in l:
               SPE = float(l[4])
        return SPE

    def get_energy_from_aoforce_out(self, aoforce_out_name = 'aoforce.out'):
        aoforce_path = os.path.join(self.path, aoforce_out_name)
        with open(aoforce_path) as aoforce:
            for line in aoforce:
                if '  zero point VIBRATIONAL energy  ' in line:
                    ZPE = float(line.split()[6]) # The zero point vibrational energy in Hartree
                if 'SCF-energy' in line:
                    SPE = float(line.split()[3]) # Energy in Hartree
        return SPE, ZPE

    def ridft(self):
        os.chdir(self.path)
        os.system("ridft > ridft.out")
        SPE =  self.get_energy_from_ridft_out()
        os.chdir(self.maindir)
        return SPE

    def change_dsta_to(self, dsta = 1.0):
        ''' Changes the start value for SCF damping. 
        '''
        control_path = os.path.join(self.path,'control')
        if not os.path.isfile(control_path):
            raise FileNotFoundError("There is no control file in the directory %s." % self.path)
    
        newlines = []
        with open(control_path,'r') as control:
            for i,line in enumerate(control):
                if '$scfdamp' in line:
                    newline = line.split()[0]+' start=%f ' %dsta + " ".join(line.split()[2:])+'\n'
                    newlines.append(newline)
                else:
                    newlines.append(line)
    
        os.remove(control_path)
        with open(control_path,'w') as new_control:
            for line in newlines:
                new_control.write(line)
        return


    def kdg(self, dg_name=""):
        """Removes the data group named <dg_name>."""
        os.chdir(self.path)
        os.system("kdg %s" % dg_name)
        print("The data group %s is removed from the control file." %dg_name)
        os.chdir(self.maindir)
        return


    def similaritycheck_from_mol(self, mol_1, mol_2):
         xyz_1 = os.path.join(self.path,'mol_1.xyz')
         mol_1.write(xyz_1)
         xyz_2 = os.path.join(self.path,'mol_2.xyz')
         mol_2.write(xyz_2)
         is_similar = self.similaritycheck_from_xyz(xyz_1, xyz_2)
         os.remove(xyz_1)
         os.remove(xyz_2)
         return is_similar


    def similaritycheck_from_xyz(self, xyz_1, xyz_2):
        out = os.popen('''echo "'%s' '%s'" | similaritycheck''' %(xyz_1, xyz_2)).read()
        is_similar = False
        if out.split()[-1] == 'T':
           is_similar = True
        if os.path.isfile('ISSIMILAR'):
           os.remove('ISSIMILAR')
        return is_similar




class GeometryTools:

    def add_noise(mol, active_atoms = [], upto = 0.1):
        """ Returns xyz coordinates with a noise using a uniform distribution.

        First shifts the noise to the origin by subtracting 0.5,
        then divides it by 10 to get a noise up to 0.1 Angstrom by default.

        No noise is added to the active atoms.
        """
        if active_atoms != []:
            noise_ini = (numpy.random.rand(mol.natoms-len(active_atoms),3)-0.5)*upto
            noise = GeometryTools._make_noise_without_active_atoms(mol, noise_ini, active_atoms) 
        else:
            noise = (numpy.random.rand(mol.natoms,3)-0.5)*upto
        # add the noise array to the coordinates array
        new_xyz = mol.xyz + noise
        return new_xyz

    def _make_noise_without_active_atoms(mol, noise_ini, active_atoms):
        """ For the case of active atoms (in transition states), makes 
        sure that there is no noise added to the coordinates of these atoms.
        """
        noise_list = noise_ini.tolist()
        noise_to_add = []
        j = 0
        for i in range(mol.natoms):
            if i+1 in active_atoms:
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
        # 1. Create a tmp directory
        curdir = os.getcwd()
        tmp_path = os.path.join(curdir,"tmp")
        os.mkdir(tmp_path)
        GeneralTools(tmp_path).write_coord_from_mol(mol)
        # 2. Write a define file 
        GeometryTools._write_define_check_symmetry(tmp_path)
        # 3. Invoke define 
        os.chdir(tmp_path)
        os.system("define < define.in > define.out")
        os.chdir(curdir)
        # 4. Read the assigned point group
        point_group = GeometryTools.get_point_group_from_control(tmp_path)
        # 5. Remove the tmp directory
        shutil.rmtree(tmp_path)
        return point_group

    def get_point_group_from_coord(path=os.getcwd(), coord_name='coord'):
        mol = GeneralTools(path).coord_to_mol(coord_name)
        point_group = GeometryTools.get_point_group_from_mol(mol)
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


    def _write_define_new_point_group(new_point_group='c1', path=os.getcwd()):
        # 1. Get the number of alpha and beta electrons
        nalpha, nbeta = GeneralTools(path).get_nalpha_and_nbeta_from_ridft_output()
        # 2. Get the charge
        charge = GeneralTools(path).read_charge_from_control()
        # 3. Write a define input
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
        f.write("%d\n" % charge)
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

    def change_point_group(path=os.getcwd(), point_group='c1'):
        non_abelian_point_groups = {  'o':'d2',  'oh':'d2h',  'td':'c2v',  'th':'s6',    't':'d2',
                                'd2d':'c2v','d3d':'s6',  'd4d':'s8',  'd5d':'s10', 'd6d':'s12','d7d':'s14','d8d':'s16',
                                'd3h':'c2v','d4h':'c4h', 'd5h':'c2v', 'd6h':'c2h', 'd7h':'c2v', 'd8h':'c2h',
                                'c3v':'cs', 'c4v':'c2v', 'c5v':'c5',  'c6v':'c2v', 'c7v':'c7', 'c8v':'c8',
                                 'd3':'c3',  'd4':'c4',   'd5':'c5',   'd6':'c6',   'd7':'c7',  'd8':'c8'}
        if point_group in non_abelian_point_groups:
            point_group_to_assign = non_abelian_point_groups[point_group] 
        else:
            point_group_to_assign = point_group
        # 1. Get the number of atoms
        natoms = GeometryTools._get_natoms_from_control(path)
        # 2. Write the define.in file
        GeometryTools._write_define_new_point_group(point_group_to_assign, path)
        # 3. Invoke define
        GeneralTools(path).invoke_define(define_out_name = "define-sym.out")
        # 4. Check if the number of atoms is still the same
        newnatoms = GeometryTools._get_natoms_from_control(path)
        if natoms != newnatoms:
            raise SymmetryAssignmentChangedtheStructureError("The structure is does not follow the  point group %s symmetry operations. Therefore, while trying to change the symmetry group, new atoms are added to enforce it." % new_point_group)
        return point_group_to_assign

        

class Mol:

    def __init__(self, mol):
        self.mol = mol
        return

    def count_number_of_electrons(self, charge=0):
        """ Counts the number of electrons. """
        # The dictionary of number of electrons
        nelectrons = {'H':1,'He':2, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'S':16, 'Cl':17, 'Ar':18} 
        nel = 0
        # 1) Count the number of electrons in the system
        for t in set(self.mol.elems):
           amount = self.mol.elems.count(t)
           nel += nelectrons[t.capitalize()]*amount
        # 2) Account for the charge of the molecule
        nel -= charge
        return nel

    def make_molecular_graph(self):
        # 1) Detect the connectivity information by bond order
        self.mol.detect_conn_by_bo()
        # 2) Make the molecular graph
        self.mol.addon("graph")
        self.mol.graph.make_graph()
        mg = self.mol.graph.molg
        return mg

    def separate_molecules(self):
       """Returns a dictionary of mol objects."""
       self.make_molecular_graph()
       mg = self.mol.graph.molg
       # 1) Label the components of the molecular graph to which each vertex in the the graph belongs
       from graph_tool.topology import label_components
       labels = label_components(mg)[0].a.tolist()

       # 2) Number of molecules
       n_mols = len(set(labels))

       # 3) Now create mol objects with the separated molecules and append them into a list
       mols = []
       for i in set(labels):
           n_atoms = labels.count(i)
           mol_str = '%d\n\n' %n_atoms
           counter = 0
           for j,label in enumerate(labels):
               if i == label:
                   mol_str += '%s %5.10f %5.10f %5.10f' %(self.mol.elems[j], self.mol.xyz[j,0], self.mol.xyz[j,1], self.mol.xyz[j,2])
                   counter += 1
                   if counter != n_atoms:
                       mol_str += '\n'
           mol_tmp = molsys.mol.from_string(mol_str, 'xyz')
           mol_tmp.detect_conn_by_bo()
           mols.append(mol_tmp)
       return mols

    def get_max_M(self):
        """ Calculates the number of core + valence MOs (nMOs)
        and from there calculates the maximum multiplicity that molecule possibly have.
        """
        # The dictionary of the minimal number of atomic orbitals
        # i.e. core + valence shells
        nAOs = {'H':1,'He':1, 'C':5, 'N':5, 'O':5, 'F':5, 'Ne':5, 'S':9, 'Cl':9, 'Ar':9}
        nMOs = 0       
        for t in set(self.mol.elems):
           amount = self.mol.elems.count(t)
           nMOs += nAOs[t.capitalize()]*amount
        nel = self.count_number_of_electrons()
        alphashells = nMOs
        betashells = nel - nMOs
        NumberofMaxUnpairedElectrons = alphashells - betashells
        Max_M = NumberofMaxUnpairedElectrons + 1
        return Max_M


 
class OptimizationTools:
    """Methods for the optimization of QM species with ridft."""
    def __init__(self, path=os.getcwd(), lot = 'ri-utpss/SVP', max_mem = 500):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        else:
            self.path = os.path.abspath(path)
            self.maindir = os.getcwd()
            self.lot = lot # e.g. ri-utpss/SVP
            self.max_mem = max_mem
            self.QM_paths = {}
        return

    def get_mol_from_coord(self, path = ''):
        if path == '': path = self.path
        coord = os.path.join(path, 'coord')
        xyz   = os.path.join(path, 'coord.xyz')
        os.system("t2x %s > %s" %(coord, xyz))
        mol =  molsys.mol.from_file(xyz,'xyz')
        mol.detect_conn_by_bo()
        return mol

    def freeze_atoms(self, active_atoms, path = ''):
        """ writes a new define file to fix active atoms in internal coordinates. """
        if path == '': path = self.path
        a_a = ''
        for i in active_atoms:
            a_a += str(i)
            if i != active_atoms[-1]:
               a_a += ',' 
        define_in_path = os.path.join(path,'define.in')
        f = open(define_in_path, 'w')
        f.write(' \n')
        f.write('y\n')
        # adds the letter f to the coord file next to the active atoms
        f.write('fix %s\n' %(a_a))
        # defines internal coordinates, also taking into account of the active atoms
        f.write('ired\n')
        # removes the letter f in the coord file, so that only internal coordinates are frozen
        f.write('fix none\n')
        # leaves the geometry menu
        f.write('*\n')
        # exits define
        f.write('qq\n')
        f.close()
        os.chdir(self.path)
        os.system('define < %s > %s' %(define_in_path, os.path.join(path,'define.out')))
        os.chdir(self.maindir)
        return

    def ired_and_itvc_1(self, rmax = 3e-2, path = ''):
        """writes a new define file to define internal redundant coordinates,
           changes the itvc to 1 (for TS optimization), 
           and changes the coordinates to the redundant internal coordinates.
        """
        if path == '': path = self.path
        define_in_path = os.path.join(path,'define.in')
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
        f.write('rmax %e\n' %rmax)
        f.write('*\n')
        f.write('q\n')
        f.close()
        GeneralTools(path).invoke_define()
        return

    def change_rmax(self, rmax = 3e-2, path = ''):
        " to change the maximum thrust radius for the geometry optimizations. "
        if path == '': path = self.path
        define_in_path = os.path.join(path,'define.in')
        f = open(define_in_path, 'w')
        f.write(' \n')
        f.write(' \n')
        f.write(' \n')
        f.write(' \n')
        f.write(' \n')
        f.write(' \n')
        f.write('stp\n')
        f.write('rmax %e\n' %rmax)
        f.write('*\n')
        f.write('q\n')
        f.close()
        GeneralTools(path).invoke_define()
        return


    def ired(self, path = ''):
        """writes a new define file to define internal redundant coordinates,
           and changes the coordinates to the redundant internal coordinates.
        """
        if path == '': path = self.path
        define_in_path = os.path.join(path,'define.in')
        f = open(define_in_path, 'w')
        f.write(' \n')
        f.write('y\n')
        f.write('ired\n')
        f.write('*\n')
        f.write(' \n')
        f.write(' \n')
        f.write(' \n')
        f.write(' \n')
        f.write('q\n')
        f.close()
        GeneralTools(path).invoke_define()
        return

    def check_imaginary_frequency(self, verbose = False, path = ''):
        if path == '': path = self.path
        f_vibspec = os.path.join(path,"vibspectrum")
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
        if verbose: print('The number of imaginary frequencies is ', inum,'and they are/it is', imfreq)
        return inum, imfreq

    def jobex(self, ts=False, cycles=150, gcart=4):
        converged = False
        os.chdir(self.path)
        if ts:
            os.system("jobex -ri -c %d -trans -gcart %d > jobex.out" %(cycles,gcart))
        else:
            os.system("jobex -ri -c %d -gcart %d > jobex.out" %(cycles,gcart))
        f_path = os.path.join(self.path,"GEO_OPT_CONVERGED")
        if os.path.isfile(f_path):
                converged = True
                os.remove(f_path)
        os.chdir(self.maindir)
        return converged

    def aoforce(self):
        os.chdir(self.path)
        os.system("aoforce > aoforce.out")
        if os.path.isfile("dh"):
            os.remove("dh")
        os.chdir(self.maindir)
        return

    def IRC(self):
        os.chdir(self.path)
        os.system("DRC -i -c 150 > IRC.out")
        os.chdir(self.maindir)
        return

    def common_workflow(self, mol, title = '', add_noise = True, TS = False, M = None, active_atoms = [], fermi = True):
        ''' Adds a noise within 0.1 Angstrom to the reference structure and performs single point calculation
        M              : the spin multiplicity
        fermi          : apply Fermi smearing or not?
        '''
        # mol = molsys.mol.from_file(path_ref)
        GT = GeneralTools(self.path)
        atom = False
        if add_noise:
            xyz = GeometryTools.add_noise(mol, active_atoms = active_atoms, upto = 0.05)
        else:
            xyz = mol.xyz

        if mol.natoms == 1:
            atom = True

        #############################
        # CASES WITH FERMI SMEARING #
        #############################
        if fermi:

            # 1. Determine the start multiplicity
            nel = Mol(mol).count_number_of_electrons(charge = 0)
            if (nel % 2) == 0:
                M_start = 3
            else:
                M_start = 2

            # 2. Perform single point calculation with Fermi smeatring
            GT.make_tmole_dft_input(
                    elems = mol.elems, 
                    xyz = xyz, 
                    M = M_start, 
                    max_mem = self.max_mem, 
                    title = title, 
                    lot = self.lot, 
                    scf_dsta = 1.0, # SCF start damping
                    fermi = True, # True for Fermi smearing 
                    nue = False) # True to enforce a multiplicity in the Fermi smearing
            GT.run_tmole()
            converged = GT.check_ridft_converged()

            # 3. If SCF did not converge, increase the SCF start damping to 2.0 instead of 1.0.
            if not converged:
                print('The SCF calculation for Fermi smearing with start damping 1.0 did not converge. Increasing it to 2.0.' )
                for f in os.listdir(self.path):
                    os.remove(os.path.join(self.path,f))
                GT.make_tmole_dft_input(
                        elems = mol.elems,  
                        xyz = xyz,
                        M = M_start,
                        max_mem = self.max_mem,
                        title = title,
                        lot = self.lot,
                        scf_dsta = 2.0, # SCF start damping
                        fermi = True, # True for Fermi smearing 
                        nue = False) # True to enforce a multiplicity in the Fermi smearing
                GT.run_tmole()
                converged = GT.check_ridft_converged()
                if not converged:
                    print('The SCF calculation with Fermi smearing did not converge also with start damping 2.0.')
                    sys.exit()
                else:
                    print('The SCF calculation converged with start damping 2.0. Decreasing it back to 1.0 and re-performing the SCF calculation.')
                    GT.change_dsta_to(1.0)
                    energy = GT.ridft()
                    converged = GT.check_ridft_converged()
                    if not converged:
                        print('The SCF calculation with start damping 1.0 did not converge.')
                        sys.exit()

            # 4. Remove the data group $fermi from the control file
            GT.kdg("fermi")

            # 5. If there are partial occupations round them to integers
            GT.round_fractional_occupation()
            energy = GT.ridft()
            print("The energy of the structure is %f Hartree." %energy)

            # 6. Now get the spin multiplicity from Fermi
            nalpha, nbeta = GT.get_nalpha_and_nbeta_from_ridft_output()
            M_fermi = GT.calculate_spin_multiplicity_from(nalpha, nbeta)

            ### TS ###
            if TS:
                assert M != None, "You must provide a multiplicity for a TS!"
 
                # 7. If the multiplicity changes in the Fermi smearing, change the multiplicity such that the desired multiplicity is used.
                if M_fermi != M:
                    print("The spin multiplicity after Fermi smearing is %d, it will be changed to %d." %(M_fermi, M))
                    GT.for_c1_sym_change_multiplicity_in_control_by(M-M_fermi, nalpha, nbeta)
                    energy = GT.ridft()
                    converged = GT.check_ridft_converged()
                    if not converged:
                        print('The SCF calculation did not converge starting with orbitals from Fermi and multiplicity %d.' %M)
                        sys.exit()

            ### Equilibrium Structure ###
            else:
                # 7. Now calculate two lower spin multiplicities
                dict_energies = {energy:M_fermi}
                if nel > 1:
                    new_dirs = []
                    if M_fermi-2 > 0:
                        m2_path = os.path.join(self.path,'M_%d' %(M_fermi-2))
                        new_dirs.append(m2_path)
                        os.chdir(self.path)
                        os.system('cpc %s' %m2_path)
                        os.chdir(m2_path)
                        GT_m2 = GeneralTools(m2_path)
                        GT_m2.for_c1_sym_change_multiplicity_in_control_by(-2, nalpha, nbeta)
                        energy_m2 = GT_m2.ridft()
                        dict_energies[energy_m2] = M_fermi-2
                        os.chdir(self.maindir)
                    # The maximum possible multiplicity that the molecule can have
                    M_max = Mol(mol).get_max_M()
                    if M_fermi < M_max:
                        p2_path = os.path.join(self.path,'M_%d' %(M_fermi+2))
                        new_dirs.append(p2_path)
                        os.chdir(self.path)
                        os.system('cpc %s' %p2_path)
                        os.chdir(p2_path)
                        GT_p2 = GeneralTools(p2_path)
                        GT_p2.for_c1_sym_change_multiplicity_in_control_by(+2, nalpha, nbeta)
                        energy_p2 = GT_p2.ridft()
                        dict_energies[energy_p2] = M_fermi+2
                        os.chdir(self.maindir)
                    M_final = dict_energies[min(dict_energies)]
                    print('The dictionary of energies and multiplicities:')
                    print(dict_energies)
                    if M_final != M_fermi:
                        print('The multiplicity %d results in lower energy than %d.' %(M_final, M_fermi))
                        path_min = os.path.join(self.path,'M_%d' %(M_final))
                        path_fermi = os.path.join(self.path,'M_%d' %(M_fermi))
                        os.mkdir(path_fermi)
                        # move the files the path where the original Fermi smearing was done to a new directory.
                        for fname in os.listdir(self.path):
                            f = os.path.join(self.path,fname)
                            if os.path.isfile(f) and fname not in ['submit.py','submit.out']:
                                shutil.move(f, path_fermi)
                        # move the files of the multiplicity which gives the lowest energy to a higher directory.
                        for fname in os.listdir(path_min):
                            f = os.path.join(path_min,fname)
                            shutil.move(f, self.path)
                        shutil.rmtree(path_min)
                        print('The calculations will proceed using multiplicity %d.' % M_final )
                    else:
                        print('The multiplicity %d will be used.' %M_final)

        #################################
        # CASES WITHOUT FERMI SMEARING  #
        #################################
        # For cases without Fermi smearing starting from extented Hueckel Theory guess for a given multiplicity
        else:
            assert M != None, "You must provide a multiplicity without Fermi smearing!"
            GT.make_tmole_dft_input(
                    elems   = mol.elems,
                    xyz     = xyz,
                    M       = M,
                    max_mem = self.max_mem,
                    title   = title,
                    lot     = self.lot,
                    fermi   = False, # True for Fermi smearing
                    nue     = False) # True to enforce a multiplicity in the Fermi smearing
            GT.run_tmole()
            converged = GT.check_ridft_converged()
            # If SCF did not converge, increase the SCF start damping to 2.0 instead of 1.0.
            if not converged:
                print('The SCF calculation with start damping 1.0 did not converge. Increasing it to 2.0.' )
                for f in self.path:
                    os.remove(f)
                GT.make_tmole_dft_input(
                        elems = mol.elems,
                        xyz = xyz,
                        M = M,
                        max_mem = self.max_mem,
                        title = title,
                        lot = self.lot,
                        scf_dsta = 2.0, # SCF start damping
                        fermi = False, # True for Fermi smearing 
                        nue = False) # True to enforce a multiplicity in the Fermi smearing
                GT.run_tmole()
                converged = GT.check_ridft_converged()
                if not converged:
                    print('The SCF calculation did not converge also with start damping 2.0.')
                    sys.exit()
                else:
                    print('The SCF calculation converged with start damping 2.0. Decreasing it back to 1.0 and re-performing the SCF calculation.')
                    GT.change_dsta_to(1.0)
                    energy = GT.ridft()
                    converged = GT.check_ridft_converged()
                    if not converged:
                        print('The SCF calculation with start damping 1.0 did not converge.')
                        sys.exit()
            energy = GT.get_energy_from_ridft_out()

        converged = GT.check_ridft_converged()
        if not converged:
            print('The ridft calculation did not converge.')
            sys.exit()
        elif atom:
            os.system('touch %s' %os.path.join(self.path,'FOUND'))
        return atom, energy

    def list_of_list_rbonds(self, rbonds):
        bonds = []
        bond_tmp = []
        for i, index in enumerate(rbonds):
            if i%2 == 0:
               bond_tmp.append(index)
            else:
               bond_tmp.append(index)
               bonds.append(sorted(bond_tmp))
               bond_tmp = []
        return bonds

    def match_wrt_ts(self, mol_spec, mol_ts, label, n, rbonds, atom_ids_dict):
        ''' This method matches the indices of the mol_spec (an equilibrium species) with that of the mol_ts. 
            atoms_ids_dict   : The dictionary which holds the list of atom_ids from the ReaxFF trajectory
                               The keys of the dictionary should be as "label_n"; 
                               e.g. label = "educt", n = 1 => The corresponding key is "educt_1"
            mol_spec, mol_ts : The mol objects, with detected connectivities
            rbonds           : 
            The indices cannot be directly used to re-order the species, because they are not always written to
            the opt_species table in the database, to avoid the redundant species. Therefore, the indices of the
            species should be mapped on an the indices of an "extracted species" from the transition state
            structure.
        '''
        # 1. Convert the rbonds to a list of lists
        rbonds = self.list_of_list_rbonds(rbonds)

        # 2. Extract the "equilibrium" species from the TS structure and create a mol object
        counter = 0
        natoms = len(mol_spec.elems)
        mol_str = '%d\n\n' %natoms
        ReaxFF2match = {} # dictionary to match the indices of the ReaxFF optimized TS with the DFT optimized equilibrium species.
        match2ReaxFF = {}           
        # Loop over the indices of the ReaxFF optimized TS and its atom indices from the MD simulation
        for i_ReaxFF, i_fromMD in enumerate(atom_ids_dict['ts']):
            # Loop over the indices of the "equilibrium" species from the MD simulation
            for j in atom_ids_dict['%s_%d' %(label,n)]:

                # Match the indices of TS from MD simulation to those of the "equilibrium" species from MD simulation
                if i_fromMD == j:
                    # Mapping of the ReaxFF optimized TS indices to the extracted equilibrium structure
                    ReaxFF2match[i_ReaxFF] = counter
                    match2ReaxFF[counter] = i_ReaxFF
                    mol_str += '%s %5.10f %5.10f %5.10f\n' %(mol_ts.elems[i_ReaxFF], mol_ts.xyz[i_ReaxFF,0], mol_ts.xyz[i_ReaxFF,1], mol_ts.xyz[i_ReaxFF,2])
                    counter += 1

        # 3. Create a mol object for the extracted equilibrium structure from the TS structure.
        mol_match = molsys.mol.from_string(mol_str, 'xyz')
        mol_match.detect_conn_by_bo()

        # 4. Check if a reactive bond is within the extracted species
        # As the "equilibrium" species is extracted from the TS, the bond might not be connected.
        # This will cause a non-isomorphism; and therefore, the matching of the indices will fail.
        # To avoid this, assure that the reactive bonds are artificially connected.
        ratom = None
        rbond_btw = []
        for rb in rbonds:
            if rb[0] not in mol_ts.conn[rb[1]]:
                mol_ts.conn[rb[0]].append(rb[1])
                mol_ts.conn[rb[1]].append(rb[0])
            if set(ReaxFF2match).intersection(rb) == set(rb):
                print('The reactive bond %d---%d belongs to the %s_%d.' %(rb[0]+1,rb[1]+1,label,n))
                i0 = ReaxFF2match[rb[0]]
                i1 = ReaxFF2match[rb[1]]
                if i0 not in mol_match.conn[i1]:
                    mol_match.conn[i0].append(i1)
                    mol_match.conn[i1].append(i0)
            elif set(ReaxFF2match).intersection(rb) != set():
                 ratom = list(set(ReaxFF2match).intersection(rb))[0]
                 rbond_btw = rb
                 print('The reactive bond %d---%d is between the two %ss and %d belongs to the %s_%d.' %(rb[0]+1,rb[1]+1,label,ratom+1,label,n))
        if n==2 and ratom == None:
            print('There is no reactive bond between the two molecules! Exiting...')
            sys.exit()

        # 5. Make molecular graphs of the DFT opt species and extracted one
        mol_spec.addon("graph")
        mol_spec.graph.make_graph()
        mg_spec = mol_spec.graph.molg
        mol_spec.graph.plot_graph('mg_%s_%d'  %(label, n))
        mol_spec.write("mol_%s_%d.xyz" %(label, n))

        mol_match.addon("graph")
        mol_match.graph.make_graph()
        mg_match = mol_match.graph.molg
        mol_match.graph.plot_graph('mg_match')
        mol_match.write("mol_match.xyz")

        # 6. Now compare the molecular graphs of the DFT optimized species and the "extracted" equilibrium structure to get the matching indices.
        # a) compare the graph
        # isomorphism is buggy ---> first indices are mapped ---> then vertex types are compared => prone to fail
        #is_equal, isomap = graph_tool.topology.isomorphism(mg_spec,mg_match, vertex_inv1=mg_spec.vp.type, vertex_inv2=mg_match.vp.type, isomap=True)
        # Therefore, as a workaround use subgraph_isomorphism...
        masterg = Graph(mg_match)
        masterg.add_vertex()
        vertex_maps = graph_tool.topology.subgraph_isomorphism(mg_spec, masterg, max_n=0, vertex_label=(mg_spec.vp.type,masterg.vp.type), edge_label=None, induced=False, subgraph=True, generator=False) 
        is_equal = len(vertex_maps) > 0
        if is_equal == False:
            print('The graph of the %s_%d is not isomorphic to the graph of the fragment of %s_%d in the transition state. Exiting...' %(label,n,label,n))
            sys.exit()
        isomap = vertex_maps[0]
        # b) get matching indices
        vts2vspec = {}
        print("------------------")
        print(" Matching indices")
        print("------------------")
        print("  spec   |   ts   ")
        print("------------------")
        for vspec, vmatch in zip(mg_spec.vertices(), isomap):
            vts = match2ReaxFF[vmatch]
            print("  %2s%3d  |  %2s%3d " %(mol_spec.elems[int(vspec)].capitalize(), int(vspec)+1, mol_ts.elems[vts].capitalize(), vts+1))
            vts2vspec[vts]  = int(vspec)
        return ratom, rbond_btw, vts2vspec

    def reorder_wrt_ts(self, QM_path, ts_path, label, n, rbonds, atom_ids_dict):
        # 1. Make a mol object for the species
        path_spec = os.path.join(os.path.join(QM_path,"%s_%d" %(label,n)),'coord.xyz')
        mol_spec = molsys.mol.from_file(path_spec)
        mol_spec.detect_conn_by_bo()

        # 2. Make a mol object for the TS
        mol_ts = molsys.mol.from_file(ts_path)
        mol_ts.detect_conn_by_bo()

        # 3. Get the matching indices
        ratom, rbond_btw, vts2vspec = self.match_wrt_ts(mol_spec, mol_ts, label, n, rbonds, atom_ids_dict)

        # 4. Order the species wrt TS
        mol_str = '%d\n\n' %mol_ts.natoms
        for vts in vts2vspec:
            vspec = vts2vspec[vts]
            atom = mol_spec.elems[vspec]
            x = mol_spec.xyz[vspec][0]
            y = mol_spec.xyz[vspec][1]
            z = mol_spec.xyz[vspec][2]
            mol_str += '%s %5.6f %5.6f %5.6f\n' %(atom,x,y,z)
        mol_ordered = molsys.mol.from_string(mol_str,'xyz')

        return mol_ordered


    def make_rxn_complex(self, rbonds, atom_ids_dict, label, n_eq, QM_path, ts_path, distance = 3.0):
        ''' This method adds the 2nd species based on the internal coordinates of the reference TS; e.g. ReaxFF optimized TS,
            at a distance (3.0 Angstrom by default) to the 1st species.
        '''
        if n_eq > 3:
            print('Cannot handle more than two species. Exiting...')
            sys.exit()

        # 1. Make a mol object for the TS
        mol_ts = molsys.mol.from_file(ts_path)
        mol_ts.detect_conn_by_bo()

        # 2. Loop over the equilibrium species
        for n in range(1, n_eq+1):
           print('\n================')
           print('    %s_%d' %(label,n))
           print('================')
           # a) Create a mol object of the DFT optimized species
           path_opt = os.path.join(os.path.join(QM_path,"%s_%d" %(label,n)),'coord.xyz')
           mol_opt = molsys.mol.from_file(path_opt)
           mol_opt.detect_conn_by_bo()
           natoms = len(mol_opt.elems)

           # 3. Get the matching indices
           ratom, rbond_btw, vts2vopt = self.match_wrt_ts(mol_opt, mol_ts, label, n, rbonds, atom_ids_dict)
           iopt2its = {}
           vopt2vts = {}
           for vts in vts2vopt:
               vopt = vts2vopt[vts]
               vopt2vts[vopt] = vts
               iopt2its[vopt+1] = vts+1

           if n == 1: vts2vopt_1 = vts2vopt

           # 4. Make the Z-matrix of the optimized species
           # 1st species
           if n == 1:
               # 5.a) Build construction table for the 1st species
               xyz_1 = zmat(mol_opt).to_Cartesian()
               # b) Form the Z-matrix of the 1st species (to replace the internal coordinates later)
               zmat_opt_1 = xyz_1.get_zmat()
               # c) Change indices to that of the TS using iopt2its dictionary
               const_table_1 = xyz_1.get_construction_table()
               const_table_1 = const_table_1.replace(iopt2its)
               new_index = [iopt2its[iopt] for iopt in const_table_1.index]
               const_table_1.index = new_index

           #  2nd species
           elif n == 2:
               # 5.a) Build the construction table for the 2nd species
               xyz_2 = zmat(mol_opt).to_Cartesian()
               const_table_2 = xyz_2.get_construction_table()
               # b) Make sure that the Z-matrix of the 2nd species starts with the reacting atom on the 2nd species
               first_atom_idx = const_table_2.index[0]
               ratom_idx = vts2vopt[ratom]+1
               if first_atom_idx != ratom_idx:
                   print('The Z-matrix will be modified to have the reacting atom %d on %s_%d as the first.' %(ratom+1,label,n))
                   for i, idx in enumerate(const_table_2.index):
                       if idx == ratom_idx:
                          ratom_pos = i
                   new_index = list(const_table_2.index)
                   new_index[0] = ratom_idx
                   new_index[ratom_pos] = first_atom_idx
                   idx2newidx = {}
                   for i, idx in enumerate(new_index):
                       idx2newidx[const_table_2.index[i]] = idx
                   const_table_2.index = new_index
                   const_table_2 = const_table_2.replace(idx2newidx)
               # c) Form the Z-matrix of the 2nd species (to replace the internal coordinates later)
               zmat_opt_2 = xyz_2.get_zmat(const_table_2)
               # d) Change indices to that of the TS
               const_table_2 = const_table_2.replace(iopt2its)
               new_index = [iopt2its[iopt] for iopt in const_table_2.index]
               const_table_2.index = new_index

               # 6. Build the construction table of the complex/transition state

               # a) Append the construction table of the 2nd species to that of the 1st species.
               const_table = const_table_1.append(const_table_2)

               # b) Replace the empty references based on the connectivity of the TS structure.

               # 1st atom of the 2nd molecule
               # ----------------------------

               # BOND:
               ratom_on_mol1 = list(set(rbond_btw)-{ratom})[0]
               const_table.loc[ratom+1,'b'] = ratom_on_mol1+1 
               
               # ANGLE:
               for i in mol_ts.conn[ratom_on_mol1]:
                    if i != ratom_on_mol1  and i != ratom and mol_ts.elems[i] == 'c': 
                        a = i
               try:
                   const_table.loc[ratom+1,'a'] = a + 1
               except:
                    print('The connected atom to the reacting atom %d different than the other reacting atom %d could not be assigned.' %(ratom_on_mol1, ratom))
                    sys.exit()
              
               # DIHEDRAL:
               for i in mol_ts.conn[a]:
                    if i != ratom_on_mol1 and i != ratom and i != a: 
                        d = i
               try:
                    const_table.loc[ratom+1,'d'] = d + 1
               except:
                    print('The connected atom to the atom %d, different than the reactive atom %d could not be assigned.' %(a,ratom))
                    sys.exit()

               # 2nd atom of the 2nd molecule
               # ----------------------------
               if natoms >= 2:
                   idx = const_table_2.index[1]
                   const_table.loc[idx,'a'] = a + 1
                   const_table.loc[idx,'d'] = d + 1

               # 3rd atom of the 2nd molecule
               # ----------------------------
               if natoms >= 3:
                   idx = const_table_2.index[2]
                   const_table.loc[idx,'d'] = d + 1

               # 7. Construct the Z-matrix of the TS using the construction table created
               xyz_ts = zmat(mol_ts).to_Cartesian()
               zmat_complex = xyz_ts.get_zmat(const_table)

               # 8. Replace all of the coordinates for the 1st molecule
               for i in const_table_1.index:
                   if zmat_complex.loc[i,'atom'] != zmat_opt_1.loc[vts2vopt_1[i-1]+1,'atom']:
                      print('Something went wrong with matching the atom indices. Exiting...')
                      sys.exit()
                   #zmat_complex.safe_loc[i,'atom'] = zmat_opt_1.loc[vts2vopt_1[i-1]+1,'atom']
                   zmat_complex.safe_loc[i,'bond'] = zmat_opt_1.loc[vts2vopt_1[i-1]+1,'bond']
                   zmat_complex.safe_loc[i,'angle'] = zmat_opt_1.loc[vts2vopt_1[i-1]+1,'angle']
                   zmat_complex.safe_loc[i,'dihedral'] = zmat_opt_1.loc[vts2vopt_1[i-1]+1,'dihedral']
               # 9. Replace the coordinates independent coordinates of the 2nd molecule
               for i in const_table_2.index:
                   if zmat_complex.safe_loc[i,'atom'] != zmat_opt_2.loc[vts2vopt[i-1]+1,'atom']:
                      print('Something went wrong with matching the atom indices. Exiting...')
                      sys.exit()
                   # first atom
                   if i == const_table_2.index[0]:
                       zmat_complex.safe_loc[i,'bond'] = distance
                   # second atom
                   elif i == const_table_2.index[1]:
                       zmat_complex.safe_loc[i,'bond'] = zmat_opt_2.loc[vts2vopt[i-1]+1,'bond']
                   # third atom
                   elif i == const_table_2.index[2]:
                       zmat_complex.safe_loc[i,'bond'] = zmat_opt_2.loc[vts2vopt[i-1]+1,'bond']
                       zmat_complex.safe_loc[i,'angle'] = zmat_opt_2.loc[vts2vopt[i-1]+1,'angle']
                   # rest of the 2nd molecule
                   else:
                       zmat_complex.safe_loc[i,'bond'] = zmat_opt_2.loc[vts2vopt[i-1]+1,'bond']
                       zmat_complex.safe_loc[i,'angle'] = zmat_opt_2.loc[vts2vopt[i-1]+1,'angle']
                       zmat_complex.safe_loc[i,'dihedral'] = zmat_opt_2.loc[vts2vopt[i-1]+1,'dihedral']

               # 10. Convert the Z-matrix of the complex back to the carte
               xyz_complex = zmat_complex.get_cartesian()
               print('The complex is succesfully created.')
               mol_str = '%d\n\n' %mol_ts.natoms
               for i in range(mol_ts.natoms):
                   atom = xyz_complex.loc[i+1, 'atom']
                   x = xyz_complex.loc[i+1, 'x']
                   y = xyz_complex.loc[i+1, 'y']
                   z = xyz_complex.loc[i+1, 'z']
                   mol_str += '%s %5.6f %5.6f %5.6f\n' %(atom,x,y,z)
               mol_complex = molsys.mol.from_string(mol_str,'xyz')
        return mol_complex 


    def find_end_points_from_IRC(self, IRC_path = '', displaced = 'minus', gcart = 4):
        """ Optimizes the end points of IRC as found in the directories 'displaced_minus' and
            'displaced_plus'. Separates the molecules and returns the corresponding mol objects.
        """
        if IRC_path == '': IRC_path = self.path

        # 1. Make a sub directory for the optimization
        path = os.path.abspath(os.path.join(IRC_path, 'displaced_%s' %displaced))
        os.chdir(path)
        os.system('cpc %s' %displaced)
        QM_path = os.path.join(path, displaced)

        # 2. Remove the gradient
        os.remove(os.path.join(QM_path,'gradient'))

        # 3. Define the internal coordinates
        OT = OptimizationTools(QM_path, lot = self.lot, max_mem = self.max_mem)
        OT.ired()

        # 4. Optimize the geometry
        converged = OT.jobex(gcart = gcart)
        if not converged:
            print('The geometry optimization of the end point of IRC has failed.')
            sys.exit()

        # 5. Get the molecules at the end points
        mol = OT.get_mol_from_coord()
        mols = Mol(mol).separate_molecules()

        # 6. Go to the main directory
        os.chdir(self.maindir)

        return mols, QM_path


    def compare_mols(self, mols_1, mols_2, mode = 'mg'):
        ''' Compares the list of reference molecules (mols_1) with the list of molecules to compare (mols_2).
        '''
        mols_similar = True # similarity of the list of molecules
        index_dict = {}
        n_mols_1 = len(mols_1)
        n_mols_2  = len(mols_2)        
        if n_mols_1 != n_mols_2:
            mols_similar = False
        else:
            # Loop over the reference molecules
            for i,mol_1 in enumerate(mols_1):

                mol_similar = False # similarity of individual molecules

                if mode == 'mg':
                    Mol(mol_1).make_molecular_graph()
                    mg_1 = mol_1.graph.molg

                # Loop over the molecules to compare
                for j,mol_2 in enumerate(mols_2):

                    if mode == 'mg':
                        Mol(mol_2).make_molecular_graph()
                        mg_2 = mol_2.graph.molg
                        is_equal = molsys.addon.graph.is_equal(mg_1, mg_2, use_fast_check=False)[0]
                    elif mode == 'similaritycheck':
                        is_equal = GeneralTools().similaritycheck_from_mol(mol_1, mol_2)
                    else:
                        print('Please specify a valid comparison method!')
                        sys.exit()

                    if is_equal:
                        index_dict[i] = j

                    mol_similar = mol_similar or is_equal

                # if all of the reference molecules were similar to one of the molecules to compare, it is True
                mols_similar = mols_similar and mol_similar

        return mols_similar, index_dict


    def check_end_points(self, mols_minus, mols_plus, mols_ed, mols_prod, mode = 'mg'):
        """
        Compares the molecular graphs of the output of the IRC calculation to those of reference structures.
        mols_minus    : List of mol objects created by separating the molecules from IRC output, displaced_minus
        mols_plus     : List of mol objects created by separating the molecules from IRC output, displaced_plus
        mols_ed       : List of mol objects of the reference educts   (e.g. from ReaxFF optimized structures)
        mols_prod     : List of mol objects of the reference products (e.g. from ReaxFF optimized structures)
        mode          : The comparison method: 
                        -> 'mg' molecular graph, 
                        -> 'similaritycheck' mapping the structures in 3D
                           as described in https://doi.org/10.1002/jcc.21925
        """
        is_similar = False
        match = {}
        index_dict = {}
        reason = ''
        n_mol_minus  = len(mols_minus)
        n_mol_plus   = len(mols_plus)
        n_mol_educts = len(mols_ed)
        n_mol_products = len(mols_prod)
        if (n_mol_minus == n_mol_educts and n_mol_plus == n_mol_products) or (n_mol_minus == n_mol_products and n_mol_plus == n_mol_educts):

            # 1. Compare the educts with minus and plus from IRC
            educt_minus_is_similar, index_dict_ed_minus  = self.compare_mols(mols_ed, mols_minus, mode = mode)
            educt_plus_is_similar , index_dict_ed_plus   = self.compare_mols(mols_ed, mols_plus, mode = mode)

            # 2. Compare the products with minus and plus from IRC
            prod_minus_is_similar, index_dict_prod_minus = self.compare_mols(mols_prod, mols_minus, mode = mode)
            prod_plus_is_similar , index_dict_prod_plus  = self.compare_mols(mols_prod, mols_plus, mode = mode)

            if (educt_minus_is_similar and prod_plus_is_similar):
                is_similar = True
                match['educt']        = 'minus'
                index_dict['educt']   =  index_dict_ed_minus
                match['product']      = 'plus'
                index_dict['product'] =  index_dict_prod_plus
            elif (educt_plus_is_similar and prod_minus_is_similar):
                is_similar = True
                match['educt']        = 'plus'
                index_dict['educt']   =  index_dict_ed_plus
                match['product']      = 'minus'
                index_dict['product'] =  index_dict_prod_minus
            else:
                reason = 'This transition state do not connect the reference educts and products.'
                print(reason)
        return is_similar, match, index_dict, reason


    def get_max_energy_struc(self, path, plot = False):
        '''get the max energy structure
        path : string  : The path to where the woelfling calculation have been performed.
        plot : boolean : True if you want to plot the energy profile
        '''
        barrierless = False

        f_woelfling_out = os.path.join(path,"woelfling_current.out")
        if not os.path.isfile(f_woelfling_out):
            print('No woelfling calculations have been performed under this directory.')
            sys.exit()
        else:
            x = [] # structure number
            y = [] # energy profile
            with open(f_woelfling_out) as woelfling_out:
               energy_profile = {}
               for lines in woelfling_out:
                   if 'structure ' in lines:
                       line = lines.strip().split()
                       energy   = float(line[5])
                       struc_no = int(line[1])
                       y.append(energy)
                       x.append(struc_no)
                       energy_profile[struc_no] = energy
                       
            if plot:
                plt.plot(x,y)
                plt.ylabel('Energy Profile (Hartree)')
                plt.savefig(os.path.join(path,'energy_profile.pdf'), format='pdf')

            # This function takes a 1-D array and finds all local maxima by simple comparison of neighboring values. 
            from scipy.signal import find_peaks
            peaks, _ = find_peaks(y)
            print(y,x)
            if len(peaks) == 0:
                barrierless = True
            elif len(peaks) >  1: 
                print('WARNING!!! There are multiple local maxima! Only the highest one will be considered.') 
                print('The local maxima:')
                en_peaks = []
                en2peak = {}
                for peak in peaks:
                    print(peak+1)
                    en = energy_profile[peak+1]
                    en_peaks.append(en)
                    en2peak[en] = peak
                max_en = max(en_peaks)
                max_struc = en2peak[max_en] + 1
            else:
                max_struc = peaks[0] + 1

        return max_struc, barrierless

    def add_woelfling_to_control(self, control_path = 'control', ninter = 24, ncoord = 2, maxit = 40):
        newlines = ''
        with open(control_path) as control:
           for line in control:
               print(line)
               if '$end' in line:
                   newlines += '$woelfling\n ninter  %d\n riter   0\n ncoord  %d\n align   0\n maxit   %d\n dlst    3.0\n thr     1.0E-4\n method  q\n$end' %(ninter, ncoord, maxit)
               else:
                   newlines += line
        f = open(control_path,'w')
        f.write(newlines)
        f.close()
        return

    def woelfling_workflow(self, woelfling_path = '', ninter = 24, ncoord = 2, maxit = 40):
        '''workflow to perform a TS search using woelfling
        '''
        if woelfling_path == '': woelfling_path = self.path
        assert os.path.isfile(os.path.join(woelfling_path,'coords')), 'Please provide coords file!'
        control_path = os.path.join(woelfling_path,'control')
        assert os.path.isfile(control_path), 'Please provide a control file!'

        # 1. Add woelfling group to the control file
        self.add_woelfling_to_control(control_path = control_path, ninter = ninter, ncoord = ncoord, maxit = maxit)

        # 2. Perform the woelfling calculation
        os.chdir(woelfling_path)
        os.system('woelfling-job -ri > woelfling.out')

        # 3. If exists, get the highest energy local maxima as a TS start guess
        try:
            max_struc, barrierless = self.get_max_energy_struc(path = woelfling_path, plot = True)
        except:
            max_struc, barrierless = self.get_max_energy_struc(path = woelfling_path, plot = False)
        print('The maximum energy structure is %d, it will perform the calculations in the directory rechnung-%d.' %(max_struc, max_struc))

        ts_path = ''
        if not barrierless:
            # Get into the corresponding directory and calculate its Hessian
            max_struc_path = os.path.join(woelfling_path, 'rechnung-%d' %(max_struc))
            ts_path = os.path.join(max_struc_path, 'ts-test')
            os.chdir(max_struc_path)
            os.system('cpc %s' %ts_path)
            self.ired(path = ts_path)
        os.chdir(self.maindir)

        return barrierless, ts_path 
#converged, ts_path


    def ts_pre_optimization(self, mol, M, rbonds, gcart = 3, add_noise = True):
        # we only fix the internal coordinates between the atoms involved in bond-order change. So convert the bonds into atoms list.
        # indexing of active_atoms goes as 1,2,3,... whereas that of rbonds goes as 0,1,2,...
        active_atoms = []
        for i in set(rbonds):
            active_atoms.append(i+1)

        # create a QM path for the pre-optimization
        QM_path = os.path.join(self.path, 'ts')
        os.mkdir(QM_path)
        OT = OptimizationTools(QM_path, lot = self.lot, max_mem = self.max_mem)
        GT = GeneralTools(QM_path)

        # Add noise to the structure and perform the single point energy calculation
        atom, energy = OT.common_workflow(mol = mol, TS = True, fermi = True, M = M, active_atoms = active_atoms, add_noise = add_noise)

        # freeze the active atoms
        OT.freeze_atoms(active_atoms)

        # pre-optimization
        converged = OT.jobex(gcart = gcart)

        if converged:
            # remove the gradient left from pre-optimization
            os.remove(os.path.join(QM_path, 'gradient'))

        # define internal coordinates without constraints and set itvc to 1.
        GT.kdg('intdef')
        GT.kdg('redundant')
        OT.ired_and_itvc_1()

        return converged, QM_path

    def optimize_irc_end_points(self, M, label = '', mols_QM_ref = [], match = {}, index_dict = {}, irc_mols = {}, irc_path = {}, gcart = 4, is_similar = False):
        displaced   = match['%s' %label]   # 'minus' path or 'plus' path
        mols        = irc_mols[displaced]  # the separated educt mol objects from the irc end points, those we want to optimize now
        path       = irc_path[displaced] # the pathway to the .../displaced_minus or .../displaced_plus directories 
        QM_ref2irc  = index_dict['%s' %label] # dictionary which matches of e.g., QM_reference educts to the corresponding IRC species
        for i, mol_QM_ref   in enumerate(mols_QM_ref):
            # 1. Get the corresponding mol and directory of the IRC calculation
            irc_mol  = mols[QM_ref2irc[i]]

            M_QM_ref = mol_QM_ref.multiplicity # Multiplicity of the QM reference structure

            # 2. If there is only one molecule and the multiplicity matches with the ref QM species, 
            if len(mols) == 1 and M_QM_ref == M:
                OT = OptimizationTools(path, lot = self.lot, max_mem = self.max_mem)
                OT.aoforce()
                inum, imfreq = OT.check_imaginary_frequency()
                if inum != 0:
                    print("This is not a minimum!")
                    sys.exit()
                else:
                    self.QM_paths['%s_%d' %(label, i)] = irc_path
            else:
                print("The multiplicity of the educt and the irc end point does not match or the number of molecules is not 1...")
                path_i = os.path.join(path, '%s_%d'  %(label, i))
                os.mkdir(path_i)
                OT = OptimizationTools(path_i, lot = self.lot, max_mem = self.max_mem)
                # Optimize the end points
                multiplicities, QM_paths, mols_opt = OT.educts_and_products_workflow(mols = irc_mol, add_noise = False, label = label, gcart = gcart)
                M_spec = multiplicities[0]
                OT_opt = OptimizationTools(QM_paths[0], lot = self.lot, max_mem = self.max_mem)
                OT_opt.aoforce()
                inum, imfreq = OT_opt.check_imaginary_frequency()
                if inum != 0:
                    print("This is not a minimum!")
                    sys.exit()
                self.QM_paths['%s_%d' %(label, i)] = QM_paths[0]
                if M_QM_ref != M_spec:
                    print('The multiplicity of the reference QM %s %d does not match with the re-optimized IRC end point.' %(label, i))
        return


    def ts_workflow(self, QM_path_ts, mols_ed_QM_ref, mols_prod_QM_ref, gcart = 4, M = None, mode = 'mg'):
        converged = False
        is_similar = False
        reason = ''
        OT = OptimizationTools(QM_path_ts, lot = self.lot, max_mem = self.max_mem)

        # 1. Calculate the Hessian
        OT.aoforce()
        
        # 2. Check the number of imaginary frequencies
        inum, imfreq = OT.check_imaginary_frequency()

        if inum == 0:        
            reason += 'No imaginary frequency at the start structure.'
        else:
            if inum > 1:
                reason += 'There are more than one imaginary frequencies at the start structure. But it will try to optimize.'
                print(reason)

            # 3. Optimize the TS with eigenvector following
            converged = OT.jobex(ts=True, gcart = gcart)
 
            if not converged:
               reason += 'The transition state optimization did not converge.'

            else:
               # 4. Calculate Hessian and check if it is a saddle point
               OT.aoforce()
               inum, imfreq = OT.check_imaginary_frequency()

               if inum == 1: # if a saddle point
                  print('There is only one imaginary frequency. The intrinsic reaction coordinate is being calculated.')

                  # 5. Perform intrinsic reaction coordinate calculation
                  OT.IRC()
                  irc_mols = {}
                  irc_path = {}

                  # 6. Optimize the end points
                  irc_mols['minus'], irc_path['minus'] = OT.find_end_points_from_IRC(displaced = 'minus')
                  irc_mols['plus' ], irc_path['plus' ] = OT.find_end_points_from_IRC(displaced = 'plus') 

                  # 7. Compare the end points with the educts and products
                  is_similar, match, index_dict, reason_comparison = OT.check_end_points(irc_mols['minus'], irc_mols['plus'], mols_ed_QM_ref, mols_prod_QM_ref, mode = mode)
                  reason += reason_comparison

                  # 8. If the molecular graphs are similar, then optimize the end points.
                  self.optimize_irc_end_points(M = M, label = 'educt',   mols_QM_ref = mols_ed_QM_ref,   match = match, index_dict = index_dict, irc_mols = irc_mols, irc_path = irc_path, gcart = gcart, is_similar = is_similar)
                  self.optimize_irc_end_points(M = M, label = 'product', mols_QM_ref = mols_prod_QM_ref, match = match, index_dict = index_dict, irc_mols = irc_mols, irc_path = irc_path, gcart = gcart, is_similar = is_similar)

        return converged, is_similar, reason


    def educts_and_products_workflow(self, mols, add_noise = True, up_to = 0.1, label = 'eq_spec', gcart = 4):
        ''' This is meant to be called from the reaction_workflow method. But can also probably called separately.
            mols     : The mol objects of the structures to be optimised.
            add_noise: Should a noise be added to the structures?
            up_to    : How much maximum noise should be added to the structures? (in Angstrom)
            label    : Label to name the directories; e.g., 'product', 'educt', etc...
            gcart    : Threshold for geometry optimization, converge maximum norm of cartesian gradient up to 10^(-gcart) atomic units.
        '''
        multiplicities = []
        QM_paths       = []
        mols_opt        = []
        for i, mol_ini in enumerate(mols):
            QM_path = os.path.join(self.path, '%s_%d' %(label,i+1))
            os.mkdir(QM_path)
            OT = OptimizationTools(QM_path, lot = self.lot, max_mem = self.max_mem)
            GT = GeneralTools(QM_path)

            # 1. Make molecular graph of the reference structure
            Mol(mol_ini).make_molecular_graph()
            mg_ini = mol_ini.graph.molg            

            # 2. Add noise to the structure and perform the single point energy calculation
            atom, energy = OT.common_workflow(mol = mol_ini, add_noise = add_noise)

            # 3. Get the multiplicity
            nalpha, nbeta = GeneralTools(QM_path).get_nalpha_and_nbeta_from_ridft_output()
            M = abs(nalpha-nbeta)+1

            if not atom:
                # 4. Perform the geometry optimization
                converged = OT.jobex(gcart = gcart)
                if not converged:
                    print('The geometry optimization did not converge for %s_%d.' %(label,i+1))
                    exit()
 
                # 5. Make mol object after the geometry optimization
                mol_opt = OT.get_mol_from_coord()

                # 6. Check if the molecule stays as a single molecule, e.g., the awkward H-O-H-O-H-O complexes of ReaxFF...
                mols =  Mol(mol_opt).separate_molecules()
                if len(mols) != 1:
                    print('The optimised structure has more than a single molecule. This cannot be handled automatically.')
                    exit()

                # 7. Compare the graph of the initial and optimized structure 
                Mol(mol_opt).make_molecular_graph()
                mg_opt = mol_opt.graph.molg
                equal = molsys.addon.graph.is_equal(mg_ini, mg_opt)[0]

                # 8. If the graph is different, try to increase multiplicity by two
                if not equal:
                    print('The graph changes. The program will try multiplicity %d.' %(M+2))
                    path_tmp = os.path.join(QM_path, 'M_%d' %(M+2))
                    if os.path.isdir(path_tmp):
                        OT_tmp = OptimizationTools(path_tmp, lot = self.lot, max_mem = self.max_mem)
                        converged = OT_tmp.jobex(gcart = gcart)
                        if not converged:
                            print('The geometry optimization with multiplicity %d has failed.' %(M+2))
                            exit()
                        mol_tmp = OT_tmp.get_mol_from_coord()
                        Mol(mol_tmp).make_molecular_graph()
                        mg_tmp = mol_tmp.graph.molg
                        equal = molsys.addon.graph.is_equal(mg_ini, mg_tmp)[0]
                        if not equal:
                             print('The graph still changes with the multiplicity %d.' %(M+2))
                             exit()
                        else:
                             print('The graph does not change with the multiplicity %d.' %(M+2))
                             mol_opt = mol_tmp
                             files = os.listdir(QM_path)
                             os.mkdir(os.path.join(QM_path,'M_%d' %M))
                             for f in files:
                                 f_to_move = os.path.join(QM_path,f)
                                 if os.path.isfile(f_to_move) and f not in ['submit.py','submit.out']:
                                     shutil.move(f_to_move, os.path.join(QM_path,'M_%d' %M))
                             for f in os.listdir(path_tmp):
                                 shutil.move(os.path.join(path_tmp, f), QM_path)
                             shutil.rmtree(path_tmp)
                             M = M+2

           # 9. Return the final multiplicities and the QM paths of each structure as a list
            mol_opt.multiplicity = M
            multiplicities.append(M)
            QM_paths.append(QM_path)
            mols_opt.append(mol_opt)

        return multiplicities, QM_paths, mols_opt

    def write_coords(self, coords_path, mol_ini, mol_fin):
        coords = open(coords_path,'w')
        GT = GeneralTools()
        coords.write(GT.mol_to_coord(mol_ini))
        coords.write('\n')
        coords.write(GT.mol_to_coord(mol_fin))
        coords.close()
        return


    def reaction_workflow(self, rbonds = [], path_ref_educts = [], path_ref_products = [], path_ref_ts = '', atom_ids_dict = {}, gcart = 4):
        ''' This method considers the reaction event and optimizes the species accordingly.
            All of the input variables are retrieved from the RDB database, but should also work if one would like to just provide some reference structures...
            rbonds           : list of integers : The indices of the reactive bonds in the TS structure. The atom indexing is like in python, 0,1,2,...
            path_ref_educts  : list of strings  : The list of paths to the reference educts cartesian coordinate files.
            path_ref_products: list of strings  : The list of paths to the reference products cartesian coordinate files.
            path_ref_ts      : string           : The path to the TS cartesian coordinate files.
        '''
        QM_path = self.path

        GT = GeneralTools(QM_path)

        n_ed   = len(path_ref_educts)
        n_prod = len(path_ref_products)

        # 1. Optimize the educts and the products
        mols_ini_ed   = []
        mols_ini_prod = []
        for path_ref in path_ref_educts:
            mol_ini = molsys.mol.from_file(path_ref) 
            mol_ini.detect_conn_by_bo()
            mols_ini_ed.append(mol_ini)
        for path_ref in path_ref_products:
            mol_ini = molsys.mol.from_file(path_ref)
            mol_ini.detect_conn_by_bo()
            mols_ini_prod.append(mol_ini)
        multiplicities_ed,   QM_paths_ed  , mols_opt_ed    = self.educts_and_products_workflow(mols = mols_ini_ed,   label = 'educt', add_noise = True, gcart = gcart)
        multiplicities_prod, QM_paths_prod, mols_opt_prod  = self.educts_and_products_workflow(mols = mols_ini_prod, label = 'product',  add_noise = True, gcart = gcart)

        unimolecular = False
        if n_ed == 1 and n_prod == 1: unimolecular = True

        # 2. Determine the multiplicity of the reaction
        M_ed   = sum(multiplicities_ed)   - n_ed   + 1 
        M_prod = sum(multiplicities_prod) - n_prod + 1

        if unimolecular:
            if M_ed != M_prod:
                print('The multiplicities of the educt and product is not the same. Not implemented yet!')
                exit()
            else:
                M = M_ed
        else:
            M = min(M_ed, M_prod)
        print("The multiplicity of the reaction is assigned as %d." %M)


        print("=========== TS BY ONLY EIGENVECTOR FOLLOWING ============")
        # 3. Pre-optimize the TS contraining the atoms of the reactive bonds
        mol_ini_ts = molsys.mol.from_file(path_ref_ts)
        mol_ini_ts.detect_conn_by_bo()
        converged, QM_path_ts = self.ts_pre_optimization(mol = mol_ini_ts, M = M, rbonds = rbonds, gcart = 3, add_noise = True)

        # ( initiate the woelfling calculation using the orbitals from the pre-optimized TS
        try_woelfling = False
        woelfling_path = os.path.join(self.path,'woelfling')
        os.chdir(QM_path_ts)
        os.system("cpc %s" %woelfling_path)
        os.chdir(self.maindir)
        # )

        if not converged:
            # -> If not converged set up a woelfling calculation.
            print('The TS pre-optimization did not converge.') 
            try_woelfling = True
        else:
            # 4. Optimize the TS
            # * Calculate the Hessian of the pre-optimized structure.
            # * If the number of imaginary frequencies are not zero, try to optimize the geometry of the transition state with eigenvector following.
            # * If the geometry optimization is converged, perform intrinsic reaction coordinate (IRC) calculation and optimize the end points.
            # * Compare the end points with from the IRC with the educts end products under QM_paths_ed and QM_paths_prod
            #    NOTE: This could also be done as described in  https://doi.org/10.1002/jcc.21925
            #          but at the moment done based on the comparison of the molecular graph
            converged, is_similar, reason = self.ts_workflow(QM_path_ts = QM_path_ts, mols_ed_QM_ref = mols_opt_ed, mols_prod_QM_ref = mols_opt_prod, M = M, gcart = gcart, mode = 'mg')
            print(reason)
            # * If not converged then set up a woelfling calculation.
            if not converged:
                try_woelfling = True
            else:
                if is_similar:
                    print('The TS connecting the minima is found.')
                    print('The TS is under:', QM_path_ts)
                    print('The educts are under:', self.QM_path_ed)
                    print('The products are under:', self.QM_path_prod)
                    # TODO: Add to the database.
                else: 
                    try_woelfling = True
                    print('A TS is found but does not connect the original minima.')

        # 5. Perform the woelfling calculation
        if try_woelfling:
            print("============ WOELFLING CALCULATION =============")
            # the coords file
            coords_path = os.path.join(woelfling_path,'coords')
            if unimolecular:
                mol_ed = self.reorder_wrt_ts(QM_path, path_ref_ts, 'educt', 1, rbonds, atom_ids_dict)
                mol_prod = self.reorder_wrt_ts(QM_path, path_ref_ts, 'product', 1, rbonds, atom_ids_dict)
                is_similar = GT.similaritycheck_from_mol(mol_ed, mol_prod)
                if is_similar:
                    print('The educt and the product is the same. A woelfling calculation is not possible.')
                    sys.exit()
                self.write_coords(coords_path, mol_ed, mol_prod)
            elif n_ed > 1 and n_prod == 1:
                mol_ed_complex = self.make_rxn_complex(rbonds, atom_ids_dict, 'educt', n_ed, self.path, path_ref_ts)
                mol_prod = self.reorder_wrt_ts(QM_path, path_ref_ts, 'product', 1, rbonds, atom_ids_dict)
                self.write_coords(coords_path, mol_ed_complex, mol_prod)
            elif n_ed > 1 and n_prod > 1:
                mol_ed_complex = self.make_rxn_complex(rbonds, atom_ids_dict, 'educt', n_ed, self.path, path_ref_ts)
                mol_prod_complex = self.make_rxn_complex(rbonds, atom_ids_dict, 'product', n_prod, self.path, path_ref_ts)
                self.write_coords(coords_path, mol_ed_complex, mol_prod_complex)
            elif n_ed == 1 and n_prod > 1:
                mol_ed = self.reorder_wrt_ts(QM_path, path_ref_ts, 'educt', 1, rbonds, atom_ids_dict)
                mol_prod_complex = self.make_rxn_complex(rbonds, atom_ids_dict, 'product', n_prod, self.path, path_ref_ts)
                self.write_coords(coords_path, mol_ed, mol_prod_complex)
            barrierless, ts_path = self.woelfling_workflow(woelfling_path)
            print(barrierless, ts_path)
            if barrierless:
                print('The reaction is barrierless.')
                # TODO: Add to the database
            else:
                converged, is_similar, reason = self.ts_workflow(QM_path_ts = ts_path, mols_ed_QM_ref = mols_opt_ed, mols_prod_QM_ref = mols_opt_prod, gcart = gcart, M = M, mode = 'mg')
                print(reason)
                # * If not converged then set up a woelfling calculation.
                if not converged:
                    print('The TS guess from woelfling calculation did not converge.')
                else:
                    if is_similar:
                        print('The TS connecting the minima is found.')
                        print('The TS is under:', QM_path_ts)
                        print('The equilibrium species:')
                        for eq_spec in self.QM_paths:
                            print(self.QM_paths[eq_spec])
                        # TODO: Add to the database.
                    else:
                        print('A TS is found but does not connect the original minima.')
                
        else:
            shutil.rmtree(woelfling_path)
        return

 

    def write_submit_py(self, rbonds = [], path_ref_educts = [], path_ref_products = [], path_ref_ts = '', atom_ids_dict = {}):
        ''' In order to write a script which will run the desired routine, in this case, the reaction_workflow. Therefore, it needs to be modified if some other task is needed.
            This can later be used to submit jobs to the queuing system by writing a job submission script. See in the Slurm class the write_submission_script.
        '''
        f_path = os.path.join(self.path,"submit.py")
        f = open(f_path,"a")
        f.write("import os\n")
        f.write("import molsys\n")
        f.write("from molsys.util import turbomole\n\n")
        f.write("lot               = '%s'\n" %self.lot)
        f.write("max_mem           = %d\n" %self.max_mem)
        f.write("rbonds            = %s\n"   %str(rbonds))
        f.write("path_ref_educts   = %s\n"   %str(path_ref_educts))
        f.write("path_ref_products = %s\n"   %str(path_ref_products))
        f.write("path_ref_ts       = '%s'\n" %path_ref_ts)
        f.write("atom_ids_dict     = %s\n" %str(atom_ids_dict))
        f.write("OT = turbomole.OptimizationTools(lot = lot, max_mem = max_mem)\n")
        f.write("OT.reaction_workflow(rbonds, path_ref_educts, path_ref_products, path_ref_ts, atom_ids_dict)\n")
        f.close()
        return

class Slurm:

    def get_avail_nodes():
        ''' Returns a list of available nodes.
        '''
        sinfo = os.popen('sinfo --format="%n %t %c %m"').read()
        n_t_c_m = sinfo.split("\n")
        avail_nodes = []
        for lines in n_t_c_m:
            line = lines.split()
            if len(line)==4 and line[1] == "idle":
                avail_nodes.append((line[0],int(line[2]),int(line[3])))
        return avail_nodes


    def get_avail_nodes_of_(partition="normal"):
        ''' Returns a list of available nodes under that partition.
        '''
        sinfo = os.popen('sinfo --format="%n %t %P"').read()
        n_t_P = sinfo.split("\n")
        avail_nodes = []
        for lines in n_t_P:
            line = lines.split()
            if len(line)==3 and line[1] == "idle" and line[2].startswith(partition):
                avail_nodes.append(line[0])
        return avail_nodes

    def get_partition_info(partition="normal"):
        ''' Returns the number of CPUs and memory of a node which belongs to that partition.
        '''
        sinfo = os.popen('sinfo --format="%P %c %m"').read()
        P_c_m = sinfo.split("\n")
        for lines in P_c_m:
            line = lines.split()
            if len(line)==3 and line[0].startswith(partition):
                CPUS = int(line[1])
                MEM  = int(line[2])
        return CPUS, MEM

    def write_submission_script(path=os.getcwd(), TURBODIR="", ntasks=8, partition="normal"):
        ''' Writes a SLURM job submission script which will run whatever is written in submit.py. See write_submit_py under OptimizationTools.
        '''
        s_script_path = os.path.join(path, "submit.sh")
        s_script = open(s_script_path,"w")
        s_script.write("#!/bin/bash\n")
        s_script.write("#SBATCH --ntasks=%d\n" %ntasks)
        s_script.write("#SBATCH --nodes=1\n")
        s_script.write("#SBATCH --error=job.%J.err\n")
        s_script.write("#SBATCH --output=job.%J.out\n")
        s_script.write("#SBATCH --time=999:00:00\n")
        s_script.write("#SBATCH --partition=%s\n" %partition)
        s_script.write("#=====================================\n")
        s_script.write("# Setup the environment for Turbomole \n")
        s_script.write("#=====================================\n")
        s_script.write("export PARA_ARCH=SMP\n")
        s_script.write("export OMP_NUM_THREADS=%d\n" %ntasks)
        s_script.write("export PARNODES=%d\n" %ntasks)
        s_script.write("export TURBODIR=%s\n" %TURBODIR)
        s_script.write("export PATH=$TURBODIR/bin/`sysname`:$PATH\n")
        s_script.write("export PATH=$TURBODIR/scripts:$PATH\n")
        s_script.write("export TURBOTMPDIR=$TMP\n")
        s_script.write("#=====================================\n")
        s_script.write("#  Copy every file and run the job    \n")
        s_script.write("#=====================================\n")
        s_script.write("sbcast submit.py $TMPDIR/submit.py\n")
        s_script.write("cd $TMPDIR\n")
        s_script.write("python3 submit.py > submit.out\n")
        s_script.write("#=====================================\n")
        s_script.write("# Copy everything back                \n")
        s_script.write("#=====================================\n")
        s_script.write("cp -r $TMPDIR/* $SLURM_SUBMIT_DIR/\n")
        s_script.write("exit\n")
        s_script.close()
        return



class Harvest:
    
    def __init__(self, path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        else:
            self.path = path
        return





