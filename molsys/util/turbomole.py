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
        coord_path = os.path.join(self.path,coord_name)
        define_in_path = os.path.join(self.path,define_in_name)
        if not os.path.isfile(coord_path):
            raise FileNotFoundError("Please provide a coord file for invoking define.")
        if not os.path.isfile(define_in_path):
            raise FileNotFoundError("Please provide an input file for invoking define.")
        out_path = os.path.join(self.path, define_out_name)
        err_path = os.path.join(self.path, "errormessage")
        os.system('define < %s > %s 2> %s' %(define_in_path, out_path, err_path))
        with open(err_path) as err:
            for line in err:
                if "define ended normally" not in line:
                    raise DefineEndedAbnormallyError("Define ended abnormally. Check the define input under %s." % define_in_path)
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
                ### get the number of alpha and beta shell occupations ###
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

    
    def round_fractional_occupation(self):
        """ Rounds the occupations and writes a new control file. """
        control_path = os.path.join(self.path,'control')
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
                    #print(split_line, occ)
                    if abs(occ-round(occ)) > THRESHOLD:
                        new_lines[i+j] = ' %s       %s                                     ( %s )\n' %(split_line[0],split_line[1],str(round(occ)))
                    j += 1
            if '$beta shells' in line:
                j = 1
                while '$' not in lines[i+j]:
                    split_line = lines[i+j].strip().split()
                    occ = float(split_line[3])
                    #print(split_line, occ)
                    if abs(occ-round(occ)) > THRESHOLD:
                        new_lines[i+j] = ' %s       %s                                     ( %s )\n' %(split_line[0],split_line[1],str(round(occ)))
                    j += 1
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
            lines = []
            with open(control_path,'r') as control:
                for i,line in enumerate(control):
                    lines.append(line)
                    if '$alpha shells' in line:
                       line_alpha = i
                    if '$beta shells' in line:
                       line_beta  = i

            if nalpha >= nbeta:
                new_alpha = nalpha+N/2
                new_beta  = nbeta-N/2
            elif nalpha < nbeta:
                new_alpha = nalpha-N/2
                new_beta  = nbeta+N/2

            new_lines = []
            for i,line in enumerate(lines):
                if '$alpha shells' in line:
                    if new_alpha != 0:
                        new_lines.append('$alpha shells\n')
                        new_lines.append(' a       1-%d   ( 1 )\n' %(new_alpha))
                if '$beta shells' in line:
                    if new_beta != 0:
                        new_lines.append('$beta shells\n')
                        new_lines.append(' a       1-%d   ( 1 )\n' %(new_beta))
                elif i not in [line_beta, line_alpha, line_beta+1, line_alpha+1]:
                    new_lines.append(line)
            os.remove(control_path)
            with open(control_path,'w') as new_control:
                for line in new_lines:
                    new_control.write(line)
        return


    def make_tmole_dft_input(self, elems, xyz, M, max_mem, title, lot, scf_dsta = 1.0, fermi = True, nue = False):
        """Creates a tmole input called 'turbo.in' with c1 symmetry.
        Parameters
        ----------
        elems   : the list of elements
        xyz     : numpy array of shape (len(elems),3)
        M       : the (initial) spin multiplicity
        max_mem : Maximum memory per core to use in the calculations.
        title   : title of the job
        scf_dsta: start value for the SCF damping.
        lot     : The QM level of theory, must be string
        fermi   : Boolean for Fermi smearing
        """
        turbo_in_path = os.path.join(self.path,"turbo.in")
        c = xyz*angstrom
        f = open(turbo_in_path,"w")
        f.write("%title\n")
        f.write("%s\n" % title)
        f.write("%method\n")
        f.write("ENRGY :: %s [gen_mult = %d, gen_symm = c1, scf_dsta = %f, scf_msil = 1000, scf_rico = %d, for_maxc = %d]\n" %
                (lot, M, scf_dsta, 0.3*max_mem, 0.7*max_mem))
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

class GeometryTools:

    def add_noise(mol, active_atoms = [], upto = 0.1):
        """ Returns xyz coordinates with a noise using a uniform distribution.

        First shifts the noise to the origin by subtracting 0.5,
        then divides it by 10 to get a noise up to 0.1 Angstrom ny default.

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
        print(active_atoms)
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
        ### create a tmp directory ### 
        curdir = os.getcwd()
        tmp_path = os.path.join(curdir,"tmp")
        os.mkdir(tmp_path)
        #tempfile.mkdtemp()
        GeneralTools(tmp_path).write_coord_from_mol(mol)
        ### write a define file ###
        GeometryTools._write_define_check_symmetry(tmp_path)
        ### invoke define ###
        os.chdir(tmp_path)
        os.system("define < define.in > define.out")
        os.chdir(curdir)
        ### read the assigned point group ###
        point_group = GeometryTools.get_point_group_from_control(tmp_path)
        ### remove the tmp directory ###
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
        ### Get the number of atoms ###
        natoms = GeometryTools._get_natoms_from_control(path)
        ### write the define.in file ###
        GeometryTools._write_define_new_point_group(point_group_to_assign, path)
        ### invoke define ###
        GeneralTools(path).invoke_define(define_out_name = "define-sym.out")
        ### check if the number of atoms is still the same ###
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
        ### count the number of electrons in the system ###
        for t in set(self.mol.elems):
           amount = self.mol.elems.count(t)
           nel += nelectrons[t.capitalize()]*amount
        ### account for the charge of the molecule ###
        nel -= charge
        return nel

    def make_molecular_graph(self, thresh = 0.2):
        # if the connectivity information not defined before, detect it.
        #if not any(self.mol.conn):
        self.mol.detect_conn(thresh = thresh)
        self.mol.addon("graph")
        self.mol.graph.make_graph()
        return

    def separate_molecules(self, thresh = 0.2):
       """Returns a dictionary of mol objects."""
       self.make_molecular_graph()
       mg = self.mol.graph.molg
       # label the components of the molecular graph to which each vertex in the the graph belongs
       from graph_tool.topology import label_components
       labels = label_components(mg)[0].a.tolist()
       # number of molecules
       n_mols = len(set(labels))
       # now create mol objects with the separated molecules
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
           mol_tmp.detect_conn(thresh)
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
    def __init__(self, path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        else:
            self.path = os.path.abspath(path)
            self.maindir = os.getcwd()
        return

    def freeze_atoms(self, active_atoms):
        """ writes a new define file to fix active atoms in internal coordinates. """
        a_a = ''
        for i in active_atoms:
            a_a += str(i)
            if i != active_atoms[-1]:
               a_a += ',' 
        define_in_path = os.path.join(self.path,'define.in')
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
        os.system('define < %s > %s' %(define_in_path, os.path.join(self.path,'define.out')))
        os.chdir(self.maindir)
        return

    def ired_and_itvc_1(self):
        """writes a new define file to define internal redundant coordinates,
           changes the itvc to 1 (for TS optimization), 
           and changes the coordinates to the redundant internal coordinates.
        """
        os.chdir(self.path)
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
        os.chdir(self.maindir)
        return

    def ired(self):
        """writes a new define file to define internal redundant coordinates,
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
        f.write('q\n')
        f.close()
        GeneralTools(self.path).invoke_define()
        return

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
        os.remove("dh")
        os.chdir(self.maindir)
        return

    def IRC(self):
        os.system("DRC -i -c 150 > IRC.out")
        return


    def common_workflow(self, path_ref, lot, TS = False, M = None, active_atoms = [], max_mem = 500, fermi = True):
        ''' Performs single point calculation, adds a noise within 0.1 Angstrom to the reference structure. 
        TS             : for transition state or not?
        M              : the spin multiplicity
        fermi          : apply Fermi smearing or not?
        max_mem        : maximum memory per core
        path_ref       : the path to the reference structure; e.g. from ReaxFF
        lot            : string, e.g. ri-utpss/SVP 
        '''
        mol = molsys.mol.from_file(path_ref)
        GT = GeneralTools(self.path)
        atom = False
        new_xyz = GeometryTools.add_noise(mol, active_atoms = active_atoms, upto = 0.05)

        if mol.natoms == 1:
            atom = True

        # For cases with Fermi smearing
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
                    xyz = new_xyz, 
                    M = M_start, 
                    max_mem = max_mem, 
                    title = path_ref, 
                    lot = lot, 
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
                        xyz = new_xyz,
                        M = M_start,
                        max_mem = max_mem,
                        title = path_ref,
                        lot = lot,
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

        # For cases without Fermi smearing starting from extented Hueckel Theory guess for a given multiplicity
        else:
            assert M != None, "You must provide a multiplicity without Fermi smearing!"
            GT.make_tmole_dft_input(
                    mol.elems,
                    new_xyz,
                    M,
                    max_mem,
                    path_ref,
                    lot,
                    False, # True for Fermi smearing
                    False) # True to enforce a multiplicity in the Fermi smearing
            GT.run_tmole()
            converged = GT.check_ridft_converged()
            # If SCF did not converge, increase the SCF start damping to 2.0 instead of 1.0.
            if not converged:
                print('The SCF calculation with start damping 1.0 did not converge. Increasing it to 2.0.' )
                for f in self.path:
                    os.remove(f)
                GT.make_tmole_dft_input(
                        elems = mol.elems,
                        xyz = new_xyz,
                        M = M,
                        max_mem = max_mem,
                        title = path_ref,
                        lot = lot,
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


    def find_end_points_from_IRC(self):
        ### MINUS
        path_minus = os.path.join(self.path, 'displaced_minus')
        os.chdir(path_minus)
        os.system('cpc minus')
        sub_path_minus = os.path.join(path_minus, 'minus')
        os.chdir(sub_path_minus)
        # remove the gradient
        os.remove('gradient')
        # define internal coordinates
        self.ired()
        # optimize the geometry
        converged = self.jobex()
        # TODO I should check if these are converged... 
        if not converged:
            print('The geometry optimization of the end point of IRC has failed.')
        # get the molecules at the end points
        os.system("t2x %s/coord > %s/coord.xyz" %(sub_path_minus, sub_path_minus))
        opt_mol_minus =  molsys.mol.from_file('%s/coord.xyz' %sub_path_minus,'xyz')
        opt_mol_minus.detect_conn()

        mols_minus = Mol(opt_mol_minus).separate_molecules()

        ### PLUS
        path_plus  = os.path.join(self.path, 'displaced_plus')
        os.chdir(path_plus)
        os.system('cpc plus')
        sub_path_plus = os.path.join(path_plus, 'plus')
        os.chdir(sub_path_plus)
        # remove the gradient
        os.remove('gradient')
        # define internal coordinates
        self.ired()
        # optimize the geometry
        converged = self.jobex()
        # TODO I should check if these are converged... 
        if not converged:
            print('The geometry optimization of the end point of IRC has failed.')
        # get the molecules at the end points
        os.system("t2x %s/coord > %s/coord.xyz" %(sub_path_plus, sub_path_plus))
        opt_mol_plus =  molsys.mol.from_file('%s/coord.xyz' %sub_path_plus,'xyz')
        opt_mol_plus.detect_conn()
        mols_plus = Mol(opt_mol_plus).separate_molecules()

        # go to the main directory
        os.chdir(self.maindir)
        return mols_minus, mols_plus


    def check_end_points(self, mols_minus, mols_plus, path_ref_educts, path_ref_products):
        """
        Compares the molecular graphs of the output of the IRC calculation to those of reference structures.
        mols_minus        : List of mol objects created by separating the molecules from IRC output, displaced_minus
        mols_plus         : List of mol objects created by separating the molecules from IRC output, displaced_plus
        path_ref_educts   : List of path of the reference educts   (e.g. from ReaxFF optimized structures)
        path_ref_products : List of path of the reference products (e.g. from ReaxFF optimized structures)
        """
        is_similar = False
        reason = ''
        n_mol_minus  = len(mols_minus)
        n_mol_plus   = len(mols_plus)
        n_mol_educts = len(path_ref_educts)
        n_mol_products = len(path_ref_products)
        if (n_mol_minus == n_mol_educts and n_mol_plus == n_mol_products) or (n_mol_minus == n_mol_products and n_mol_plus == n_mol_educts):
            educt_minus_is_similar = True
            educt_plus_is_similar = True
            for ed in path_ref_educts:
                mol_ed   = molsys.mol.from_file(ed)
                Mol(mol_ed).make_molecular_graph()
                mg_ed = mol_ed.graph.molg

                educt_minus_tmp = False
                for mol_minus in mols_minus:
                    Mol(mol_minus).make_molecular_graph()
                    mg_minus = mol_minus.graph.molg
                    is_equal = molsys.addon.graph.is_equal(mg_ed, mg_minus, use_fast_check=False)[0]
                    if is_equal: 
                        educt_minus_tmp = educt_minus_tmp or True
                    else:
                        educt_minus_tmp = educt_minus_tmp or False
                educt_minus_is_similar = educt_minus_is_similar and educt_minus_tmp

                educt_plus_tmp  = False
                for mol_plus in mols_plus:
                    Mol(mol_plus).make_molecular_graph()
                    mg_plus = mol_plus.graph.molg
                    is_equal = molsys.addon.graph.is_equal(mg_ed, mg_plus, use_fast_check=False)[0]
                    if is_equal:
                        educt_plus_tmp = educt_plus_tmp or True
                    else:
                        educt_plus_tmp = educt_plus_tmp or False
                educt_plus_is_similar = educt_plus_is_similar and educt_plus_tmp

            prod_minus_is_similar = True
            prod_plus_is_similar = True                       
            for prod in path_ref_products:
                mol_prod = molsys.mol.from_file(prod)
                Mol(mol_prod).make_molecular_graph()
                mg_prod = mol_prod.graph.molg

                product_minus_tmp = False
                for mol_minus in mols_minus:
                    mg_minus = mol_minus.graph.molg
                    is_equal = molsys.addon.graph.is_equal(mg_prod, mg_minus, use_fast_check=False)[0]
                    if is_equal: 
                        product_minus_tmp = product_minus_tmp or True
                    else:
                        product_minus_tmp = product_minus_tmp or False
                prod_minus_is_similar = prod_minus_is_similar and product_minus_tmp

                product_plus_tmp  = False
                for mol_plus in mols_plus:
                    mg_plus = mol_plus.graph.molg
                    is_equal = molsys.addon.graph.is_equal(mg_prod, mg_plus, use_fast_check=False)[0]
                    if is_equal: 
                        product_plus_tmp = product_plus_tmp or True
                    else:
                        product_plus_tmp = product_plus_tmp or False
                prod_plus_is_similar = prod_plus_is_similar and product_plus_tmp
            if (educt_minus_is_similar and prod_plus_is_similar) or (educt_plus_is_similar and prod_minus_is_similar):
                is_similar = True
            else:
                reason = 'This transition state do not connect the reference educts and products.'
                print(reason)
        return is_similar, reason


    def getthemaxenergystruc(self, path, plot):
        '''get the max energy structure
        path : string  : The path to where the woelfling calculation have been performed.
        plot : boolean : True if you want to plot the energy profile
        '''
        # TODO: I need to exclude probably the first and last one or two points
        max_struc = None
        f_woelfling_out = os.path.join(path,"woelfling_current.out")
        if os.path.isfile(f_woelfling_out):
            with open(f_woelfling_out) as woelfling_out:
               energy_profile = {}
               for lines in woelfling_out:
                   if 'structure ' in lines:
                       line = lines.strip().split()
                       energy_profile[float(line[5])] = int(line[1]) 
            max_energy = max(energy_profile)
            max_struc = energy_profile[max_energy]
            print(max_struc)
            if plot:
                x = []
                y = []
                for en in energy_profile:
                    y.append(en)
                    x.append(energy_profile[en])
                plt.plot(x,y)
                plt.ylabel('Energy Profile (Hartree)')
                plt.savefig(os.path.join(path,'energy_profile.pdf'), format='pdf')
        else:
            print('No woelfling calculations have been performed.')
        return max_struc

    def woelfling_workflow(self, unimolecular, M, active_atoms, max_mem, lot, path_ref_ed = '', path_ref_prod = '', path_ed_complex = '', path_prod_complex = ''):
        '''workflow to perform a TS search using woelfling
        '''
        reason = ''
        maindir = os.getcwd()
        woelfling_dir = os.path.join(maindir,'woelfling')
        os.mkdir(woelfling_dir)
        os.chdir(woelfling_dir)
        if not unimolecular:
            # optimize the ReaxFF complex structures
            path_ref_ed   = os.path.join(woelfling_dir, 'educt_complex')
            path_ref_prod = os.path.join(woelfling_dir, 'product_complex')
            converged_ed,   reason_ed   = self.optimize_rxn_complex(path_ref_ed,   path_ed_complex,   M, max_mem, lot, active_atoms)
            converged_prod, reason_prod = self.optimize_rxn_complex(path_ref_prod, path_prod_complex, M, max_mem, lot, active_atoms)
        else:
            out = os.popen('similaritycheck %s %s' %(os.path.join(path_ref_ed,'coord.xyz'),os.path.join(path_ref_prod,'coord.xyz'))).read()
        os.chdir(path_ref_prod)
        os.system('cpc %s' % woelfling_dir)
        os.chdir(woelfling_dir) 
        os.system('cp %s %s' %(os.path.join(path_ref_ed,   'coord'), os.path.join(woelfling_dir, 'ini')))
        os.system('cp %s %s' %(os.path.join(path_ref_prod, 'coord'), os.path.join(woelfling_dir, 'fin')))
        os.system('cat ini fin > coords')
        # add to the control file $woelfling key
        newlines = ''
        with open('control') as control:
           for line in control:
               if 'end' in line:
                   newlines += '$woelfling\n ninter  24 \n riter   0 \n ncoord  2 \n align   0 \n maxit   40 \n dlst    3.0 \n thr     1.0E-4 \n method  q \n$end'
               else:
                   newlines += line
        f = open('control','w')
        f.write(newlines)
        f.close()
        print(newlines)
        # perform the woelfling calculation
        os.system('woelfling-job > woelfling.out')
        # determine which of the points is a TS guess structure
        max_struc = self.getthemaxenergystruc(worlfling_dir, True)
        print('The maximum energy structure is %d, it will perform the calculations in the directory rechnung-%d.' %(max_struc, max_struc))
        # get into the corresponding directory and calculate its Hessian
        os.chdir('rechnung-%d' %(max_struc))
        os.system('cpc ts-test')
        os.chdir('ts-test')
        defineinput = open('define.in', 'w')
        defineinput.write(' \n')
        defineinput.write('y\n')
        defineinput.write('ired\n')
        defineinput.write('*\n')
        defineinput.write(' \n')
        defineinput.write(' \n')
        defineinput.write(' \n')
        defineinput.write(' \n')
        defineinput.write('stp\n')
        defineinput.write('itvc\n')
        defineinput.write('1\n')
        defineinput.write(' \n')
        defineinput.write('q\n')
        defineinput.close()
        os.system('define < define.in > define.out')
        os.system('ridft > ridft.out')
        self.aoforce()
        os.chdir(self.maindir)
        return


    def make_rxn_complex(self, rbonds, atom_ids_dict, label, n_eq, QM_path, ts_path, distance = 3.0):
        ''' This method adds the 2nd species based on the internal coordinates of the reference TS; e.g. ReaxFF optimized TS,
            at a distance (3.0 Angstrom by default) to the 1st species.
        '''
        if n_eq > 3:
            print('Cannot handle more than two species. Exiting...')
            sys.exit()

        # 1. Convert the rbonds to a list of lists
        bonds = []
        bond_tmp = []
        for i, index in enumerate(rbonds):
            if i%2 == 0:
               bond_tmp.append(index)
            else:
               bond_tmp.append(index)
               bonds.append(sorted(bond_tmp))
               bond_tmp = []
        rbonds = bonds
        ratom = None

        # 2. Make a mol object for the TS
        mol_ts = molsys.mol.from_file(ts_path)
        mol_ts.detect_conn()

        # 3. Loop over the equilibrium species
        for n in range(n_eq):
           print('%s_%d' %(label,n+1))
        
           # a) Create a mol object of the DFT optimized species
           path_opt = os.path.join(os.path.join(QM_path,"%s_%d" %(label,n+1)),'coord.xyz')
           mol_opt = molsys.mol.from_file(path_opt)
           mol_opt.detect_conn()
        
           # b) Extract the "equilibrium" species from the TS structure and create a mol object
           counter = 0
           natoms = len(atom_ids_dict['%s_%d' %(label,n+1)])
           if n == 0: natoms_1 = len(atom_ids_dict['%s_%d' %(label,n+1)])
           mol_str = '%d\n\n' %natoms
           ReaxFF_match = {} # dictionary to match the indices of the ReaxFF optimized TS with the DFT optimized equilibrium species.
           # Loop over the indices of the ReaxFF optimized TS and its atom indices from the MD simulation
           for i_ReaxFF, i_fromMD in enumerate(atom_ids_dict['ts']):
               # Loop over the indices of the "equilibrium" species from the MD simulation
               for j in atom_ids_dict['%s_%d' %(label,n+1)]:

                   # c) Match the indices of TS from MD simulation to those of the "equilibrium" species from MD simulation
                   if i_fromMD == j:
                       # Mapping of the ReaxFF optimized TS indices to the extracted equilibrium structure
                       ReaxFF_match[i_ReaxFF] = counter
                       mol_str += '%s %5.10f %5.10f %5.10f\n' %(mol_ts.elems[i_ReaxFF], mol_ts.xyz[i_ReaxFF,0], mol_ts.xyz[i_ReaxFF,1], mol_ts.xyz[i_ReaxFF,2])
                       counter += 1

           # d) Create a mol object for the extracted equilibrium structure from the TS structure.
           mol_match = molsys.mol.from_string(mol_str, 'xyz')
           mol_match.detect_conn()

           # 4. Check if a reactive bond is within the extracted species
           # As the "equilibrium" species is extracted from the TS, the bond might not be connected.
           # This will cause a non-isomorphism; and therefore, the matching of the indices will fail.
           # To avoid this, assure that the reactive bonds are artificially connected.
           for rb in rbonds:
               if set(ReaxFF_match).intersection(rb) == set(rb):
                   print('The reactive bond %s belongs to the %s_%d.' %(str(rb),label,n+1))
                   i0 = ReaxFF_match[rb[0]]
                   i1 = ReaxFF_match[rb[1]]
                   if i0 not in mol_match.conn[i1]:
                       mol_match.conn[i0].append(i1)
                       mol_match.conn[i1].append(i0)
               elif set(ReaxFF_match).intersection(rb) != set():
                    ratom = list(set(ReaxFF_match).intersection(rb))[0]
                    rbond_btw = rb
                    print('The reactive bond %s is between the two %ss and %d belongs to the %s_%d.' %(str(rb),label,ratom,label,n+1))
        
           if n==1 and ratom == None:
               print('There is no reactive bond between the two molecules! Exiting...')
               sys.exit()

           # 5. Make molecular graphs of the DFT opt species and extracted one
           mol_opt.addon("graph")
           mol_opt.graph.make_graph()
           mg_opt = mol_opt.graph.molg
           mol_match.addon("graph")
           mol_match.graph.make_graph()
           mg_match = mol_match.graph.molg
        
           # 6. Now compare the molecular graphs of the DFT optimized species and the "extracted" equilibrium structure to get the matching indices.
           # a) compare the graph
           is_equal, isomap = graph_tool.topology.isomorphism(mg_opt,mg_match,isomap=True)
           if is_equal == False:
               print('The graph of the %s_%d is not isomorphic to the graph of the fragment of %s_%d in the transition state. Exiting...' %(label,n+1,label,n+1))
               sys.exit()
           # b) get matching indices
           vts2vopt = {}
           iopt2its = {}
           vopt2vts = {}
           for vopt, vmatch, vts in zip(isomap, mg_opt.vertices(), ReaxFF_match):
               vts2vopt[vts]  = vopt
               iopt2its[vopt+1] = vts+1
               vopt2vts[vopt] = vts
               if n == 0: vts2vopt_1 = vts2vopt

           # 7. Make sure that the ordering of the DFT optimized species are proper for forming a Z-matrix for the TS
   
           #    2nd species
           if n == 1 and vts2vopt[ratom] != 0:
               print('The second %s should start with the atom which belongs to the breaking/forming bond.' %label)
               print('A new re-ordered mol object will be created.')
               mol_str = '%d\n\n' %natoms
               # a) first add the reacting atom
               mol_str += '%s %5.10f %5.10f %5.10f\n' %(mol_opt.elems[vts2vopt[ratom]], mol_opt.xyz[vts2vopt[ratom],0], mol_opt.xyz[vts2vopt[ratom],1], mol_opt.xyz[vts2vopt[ratom],2])
               added = [vts2vopt[ratom]]
               print('The reactive atom %d is connected to the atoms %s.' %(vts2vopt[ratom], str(mol_opt.conn[vts2vopt[ratom]])))
               # b) then add the connected atoms
               for i in mol_opt.conn[vts2vopt[ratom]]:
                    added.append(i)
                    mol_str += '%s %5.10f %5.10f %5.10f\n' %(mol_opt.elems[i], mol_opt.xyz[i,0], mol_opt.xyz[i,1], mol_opt.xyz[i,2])
               # c) then add the rest
               for i in set(range(natoms))-set(added):
                    added.append(i)
                    mol_str += '%s %5.10f %5.10f %5.10f\n' %(mol_opt.elems[i], mol_opt.xyz[i,0], mol_opt.xyz[i,1], mol_opt.xyz[i,2])
               # d) create the new mol object
               mol_opt_new = molsys.mol.from_string(mol_str, 'xyz')
               mol_opt_new.detect_conn()
               mol_opt_new.addon("graph")
               mol_opt_new.graph.make_graph()
               mol_opt = mol_opt_new
               mg_opt_new = mol_opt.graph.molg
               # e) once more, compare the molecular graphs of the DFT optimized species and the "extracted" equilibrium structure to get the matching indices.
               is_equal, isomap = graph_tool.topology.isomorphism(mg_opt_new,mg_match,isomap=True)
               vts2vopt = {}
               vopt2vts = {}
               iopt2its = {}
               for vopt, vmatch, vts in zip(isomap, mg_opt_new.vertices(), ReaxFF_match):
                   vts2vopt[vts]  = vopt
                   vopt2vts[vopt] = vts
                   iopt2its[vopt+1] = vts + 1
               # f) if with the new index matching the ratom is not the first one (which might happen for, e.g., O2, H2O), then replace the two
               if vts2vopt[ratom] != 0:
                   print('Symmetric molecule... Graph tool could not differentiate the two equivalent vertices.')
                   vts2vopt[vopt2vts[0]] = vts2vopt[ratom]
                   iopt2its[vts2vopt[ratom]+1] = vopt2vts[0]+1
                   vopt2vts[vts2vopt[ratom]] = vopt2vts[0]
                   vts2vopt[ratom] = 0
                   iopt2its[1] = ratom+1
                   vopt2vts[0] = ratom

           # 8. Make the Z-matrix of the optimized species
           # 1st species
           if n == 0:
               # 9. a) Build construction table for the 1st species and change indices to that of the TS using iopt2its dictionary
               xyz_1 = zmat(mol_opt).to_Cartesian()
               zmat_opt_1 = xyz_1.get_zmat()
               const_table_1 = xyz_1.get_construction_table()
               const_table_1 = const_table_1.replace(iopt2its)
               new_index = [iopt2its[iopt] for iopt in const_table_1.index]
               const_table_1.index = new_index

           #  2nd species
           elif n == 1:
               # 9. b) Build construction table for the 2nd species and change indices to that of the TS using iopt2its dictionary
               xyz_2 = zmat(mol_opt).to_Cartesian()
               zmat_opt_2 = xyz_2.get_zmat()
               const_table_2 = xyz_2.get_construction_table()
               const_table_2 = const_table_2.replace(iopt2its)
               new_index = [iopt2its[iopt] for iopt in const_table_2.index]
               const_table_2.index = new_index

               # 10. Build the construction table of the transition state

               # a) Append the construction table of the 2nd species to that of the 1st species.
               const_table = const_table_1.append(const_table_2)

               # b) Replace the empty references based on the connectivity of the TS structure.

               # 1st atom of the 2nd molecule
               for i in mol_ts.conn[rbond_btw[0]]:
                    if i != rbond_btw[0] and i != ratom: a = i
               for i in mol_ts.conn[a]:
                    if i != rbond_btw[0] and i != ratom and i != a: d = i
               const_table.loc[ratom+1,'b'] = rbond_btw[0]+1
               const_table.loc[ratom+1,'a'] = a + 1
               const_table.loc[ratom+1,'d'] = d + 1

               # 2nd atom of the 2nd molecule
               if natoms >= 2:
                   const_table.loc[vopt2vts[1]+1,'a'] = a + 1
                   const_table.loc[vopt2vts[1]+1,'d'] = d + 1

               # 3rd atom of the 2nd molecule
               elif natoms >= 3:
                   const_table.loc[vopt2vts[2]+1,'d'] = d + 1

               # 11. Construct the Z-matrix of the TS using the construction table created
               xyz_ts = zmat(mol_ts).to_Cartesian()
               zmat_complex = xyz_ts.get_zmat(const_table)

               # 12. Replace all of the coordinates for the 1st molecule
               for i in const_table_1.index:
                   zmat_complex.safe_loc[i,'bond'] = zmat_opt_1.loc[vts2vopt_1[i-1]+1,'bond']
                   zmat_complex.safe_loc[i,'angle'] = zmat_opt_1.loc[vts2vopt_1[i-1]+1,'angle']
                   zmat_complex.safe_loc[i,'dihedral'] = zmat_opt_1.safe_loc[vts2vopt_1[i-1]+1,'dihedral']
               # 13. Replace the coordinates independent coordinates of the 2nd molecule
               for i in const_table_2.index:
                   if i == ratom+1:
                       zmat_complex.safe_loc[i,'bond'] = distance
                   elif i == vopt2vts[1]+1:
                       zmat_complex.safe_loc[i,'bond'] = zmat_opt_2.safe_loc[vts2vopt[i-1]+1,'bond']
                   elif i == vopt2vts[2]+1:
                       zmat_complex.safe_loc[i,'bond'] = zmat_opt_2.safe_loc[vts2vopt[i-1]+1,'bond']
                       zmat_complex.safe_loc[i,'angle'] = zmat_opt_2.safe_loc[vts2vopt[i-1]+1,'angle']
                   else:
                       zmat_complex.safe_loc[i,'bond'] = zmat_opt_2.safe_loc[vts2vopt[i-1]+1,'bond']
                       zmat_complex.safe_loc[i,'angle'] = zmat_opt_2.safe_loc[vts2vopt[i-1]+1,'angle']
                       zmat_complex.safe_loc[i,'dihedral'] = zmat_opt_2.safe_loc[vts2vopt[i-1]+1,'dihedral']

               # 14. Convert the Z-matrix of the complex back to the carte
               xyz_complex = zmat_complex.get_cartesian()

        return xyz_complex






    def optimize_rxn_complex(self, path, mfpx_complex, M, max_mem, lot, active_atoms):
        ''' optimizes the reaction complex by constraining internal coordinates of the active atoms.
        path         : string  : the path to where the optimization should be made
        mfpx_complex : string  : the path to the mfpx structure of the complex
        M            : integer : the multiplicity as determined based on that of educts and products
        '''
        converged = False
        reason = ''
        os.mkdir(path)
        os.chdir(path)
        mol = molsys.mol.from_file(mfpx_complex)
        GT = GeneralTools(path)
        GT.make_tmole_dft_input(elems = mol.elems, xyz = mol.xyz, M = M, max_mem = max_mem, title = mfpx_complex, lot = lot, fermi = True, nue = True)
        GT.run_tmole()
        GT.kdg("fermi")
        GT.ridft()
        converged = GT.check_ridft_converged()
        if not converged:
            reason = 'The self consistent field calculation did not converge for the reactive complex.'
        else:
            OptimizationTools(path).freeze_atoms(active_atoms)
            converged = OptimizationTools(path).jobex()
            if not converged:
                reason = 'The geometry optimization did not converge for the reactive complex.' 
            else:
                os.system('kdg intdef')
                os.system('kdg redundant')
                self.ired()
                converged = True
        os.chdir(self.maindir)
        return converged, reason


    def transition_state_workflow(self, active_atoms, path_ref_educts, path_ref_products):
        found = False
        reason = ""    
        n_ed = len(path_ref_educts)
        n_prod = len(path_ref_products)
        # freeze the active atoms
        self.freeze_atoms(active_atoms)
        # pre-optimization
        converged = self.jobex()
        if not converged:
            reason = 'Transition state pre-optimization did not converge.'
            print(reason)
        else:
            # remove the gradient left from the previous calculation.
            os.remove('gradient')

            # define internal coordinates without constraints and set itvc to 1.
            os.system('kdg intdef')
            os.system('kdg redundant')
            self.ired_and_itvc_1()

            # assign symmetry
            point_group_final = GeometryTools.get_point_group_from_coord(self.path,'coord')
            print("point_group_final", point_group_final)
            if point_group_final != "c1":
                point_group_assigned = GeometryTools.change_point_group(self.path, point_group_final)
                if point_group_assigned == point_group_final:
                    print("The point group is changed to %s." % point_group_final)
                else:
                    print("The molecule has point group of %s. However, the abelian point group %s is assigned." % ( point_group_final, point_group_assigned))

            # perform aoforce calculation     
            self.aoforce()
            
            # check the number of imaginary frequencies
            inum, imfreq = self.check_imaginary_frequency()
            if inum == 0:        
                reason += 'No imaginary frequency at the start structure.'
                print(reason)
            elif inum == 1:
                converged = self.jobex(ts=True)
                if not converged:
                   reason += 'The transition state optimzation did not converge.'
                   print(reason)
                else:
                   self.aoforce()
                   inum, imfreq = self.check_imaginary_frequency()
                   if inum == 1: 
                      print('There is only one imaginary frequency. The intrinsic reaction coordinate is being calculated.')
                      self.IRC()
                      mols_minus, mols_plus = self.find_end_points_from_IRC()
                      found, reason = self.check_end_points(mols_minus, mols_plus, path_ref_educts, path_ref_products)
            elif inum > 1:
                print('There are more than one imaginary frequencies. But it will try to optimize.')
                converged = self.jobex(ts=True)
                if not converged:
                   reason += 'The transition state optimzation did not converge.'
                   print(reason)
                else:
                   self.aoforce()
                   inum, imfreq = self.check_imaginary_frequency()
                   if inum == 1:
                      print('There is only one imaginary frequency. The intrinsic reaction coordinate is being calculated.')
                      self.IRC()
                      mols_minus, mols_plus = self.find_end_points_from_IRC()
                      found, reason = self.check_end_points(mols_minus, mols_plus, path_ref_educts, path_ref_products)
                   else:
                      reason += 'The final number of imaginary frequency is not 1.'
        return found, reason

    def minima_workflow(self):
        found = False
        reason = ""
        point_group_initial = GeometryTools.get_point_group_from_coord(self.path,'coord')
        print("point_group_initial", point_group_initial)
        converged = self.jobex()
        if converged:
            point_group_final = GeometryTools.get_point_group_from_coord(self.path,'coord')
            print("point_group_final", point_group_final)
            if point_group_final != "c1":
                point_group_assigned = GeometryTools.change_point_group(self.path, point_group_final)
                if point_group_assigned == point_group_final:
                    print("The point group is changed to %s." % point_group_final)
                else:
                    print("The molecule has point group of %s. However, the abelian point group %s is assigned." % (point_group_final, point_group_assigned))
                GeneralTools(self.path).ridft()
                converged = self.jobex()
                if not converged:
                    reason += "The geometry optimization is failed."
            self.aoforce()
            inum, imfreq = self.check_imaginary_frequency()
            if inum == 0:
                found = True
                print("The equilibrium structure is found succesfully!")
            else:
                reason += "There are imaginary frequencies. This is not an equilibrium structure."
                print(reason)
        else:
            reason += "The geometry optimization is failed."
        return found, reason

    def submit(self, foffset, active_atoms):
        if foffset == 0:
            found = self.transition_state_workflow(active_atoms)
        else:
            found = self.minima_workflow()
        return found

    def get_mol(self):
        os.system("t2x %s/coord > %s/coord.xyz" %(self.path, self.path))
        mol = molsys.mol.from_file(os.path.join(self.path,"coord.xyz"))
        return mol


    def ts_pre_optimization(self, path_ref, lot, M, rbonds):
        # we only fix the internal coordinates between the atoms involved in bond-order change. So convert the bonds into atoms list.
        # indexing of active_atoms goes as 1,2,3,... whereas that of rbonds goes as 0,1,2,...
        active_atoms = []
        for i in set(rbonds):
            active_atoms.append(i+1)

        # create a QM path for the pre-optimization
        QM_path = os.path.join(self.path, 'ts')
        os.mkdir(QM_path)
        OT = OptimizationTools(QM_path)
        GT = GeneralTools(QM_path)

        # Add noise to the structure and perform the single point energy calculation
        atom, energy = OT.common_workflow(path_ref = path_ref, lot = lot, TS = True, fermi = True, M = M, active_atoms = active_atoms)

        # freeze the active atoms
        OT.freeze_atoms(active_atoms)

        # pre-optimization
        converged = OT.jobex()

        if converged:
            # remove the gradient left from pre-optimization
            os.remove(os.path.join(QM_path, 'gradient'))

        # define internal coordinates without constraints and set itvc to 1.
        GT.kdg('intdef')
        GT.kdg('redundant')
        OT.ired_and_itvc_1()

        return converged, QM_path

    def ts_eigenvector_following(self, QM_path, path_ref_educts, path_ref_products):
        found = False
        reason = ''
        OT = OptimizationTools(QM_path)

        # 1) perform aoforce calculation     
        OT.aoforce()
        
        # 2) check the number of imaginary frequencies
        inum, imfreq = OT.check_imaginary_frequency()

        # case 1: 
        if inum == 0:        
            reason += 'No imaginary frequency at the start structure.'

        # case 2:
        elif inum == 1:
            converged = OT.jobex(ts=True)
            if not converged:
               reason += 'The transition state optimzation did not converge.'
            else:
               OT.aoforce()
               inum, imfreq = OT.check_imaginary_frequency()
               if inum == 1: 
                  print('There is only one imaginary frequency. The intrinsic reaction coordinate is being calculated.')
                  OT.IRC()
                  mols_minus, mols_plus = OT.find_end_points_from_IRC()
                  found, reason_comparison = OT.check_end_points(mols_minus, mols_plus, path_ref_educts, path_ref_products)
                  reason += reason_comparison

        # case 3:
        elif inum > 1:
            reason += 'There are more than one imaginary frequencies. But it will try to optimize.'
            print(reason)
            converged = OT.jobex(ts=True)
            if not converged:
               reason += 'The transition state optimization did not converge.'
            else:
               OT.aoforce()
               inum, imfreq = OT.check_imaginary_frequency()
               if inum == 1:
                  print('There is only one imaginary frequency. The intrinsic reaction coordinate is being calculated.')
                  OT.IRC()
                  mols_minus, mols_plus = OT.find_end_points_from_IRC()
                  found, reason_comparison = OT.check_end_points(mols_minus, mols_plus, path_ref_educts, path_ref_products)
                  reason += reason_comparison
               else:
                  reason += 'The final number of imaginary frequency is not 1.'
        return found, reason


    def educts_and_products_workflow(self, paths_ref, lot, label = 'eq_spec'):
        ''' This is meant to be called from the reaction_workflow method. But can also probably called separately.
            paths_ref: list of strings : The path to the coordinate files of the reference structures.
            lot      : string          : It should be given in accordance with the tmole manual.
        '''
        multiplicities = []
        QM_paths       = []
        for i, path_ref in enumerate(paths_ref):
            QM_path = os.path.join(self.path, '%s_%d' %(label,i+1))
            os.mkdir(QM_path)
            OT = OptimizationTools(QM_path)
            GT = GeneralTools(QM_path)

            # 1. Make mol object and molecular graph of the reference structure
            mol_ini = molsys.mol.from_file(path_ref)
            Mol(mol_ini).make_molecular_graph()
            mg_ini = mol_ini.graph.molg            

            # 2. Add noise to the structure and perform the single point energy calculation
            atom, energy = OT.common_workflow(path_ref = path_ref, lot = lot)

            # 3. Get the multiplicity
            nalpha, nbeta = GeneralTools(QM_path).get_nalpha_and_nbeta_from_ridft_output()
            M = abs(nalpha-nbeta)+1

            if not atom:
                # 4. Perform the geometry optimization
                converged = OT.jobex(gcart = 3)
                if not converged:
                    print('The geometry optimization did not converge for %s_%d.' %(label,i+1))
                    exit()
 
                # 5. Make mol object after the geometry optimization
                mol_fin = OT.get_mol()

                # 6. Check if the molecule stays as a single molecule, e.g., the awkward H-O-H-O-H-O complexes of ReaxFF...
                mols =  Mol(mol_fin).separate_molecules()
                if len(mols) != 1:
                    print('The optimised structure has more than a single molecule. This cannot be handled automatically.')
                    exit()

                # 7. Compare the graph of the initial and optimized structure 
                Mol(mol_fin).make_molecular_graph()
                mg_fin = mol_fin.graph.molg
                equal = molsys.addon.graph.is_equal(mg_ini, mg_fin)[0]

                # 8. If the graph is different, try to increase multiplicity by two
                if not equal:
                    print('The graph changes. The program will try multiplicity %d.' %(M+2))
                    path_tmp = os.path.join(QM_path, 'M_%d' %(M+2))
                    if os.path.isdir(path_tmp):
                        OT_tmp = OptimizationTools(path_tmp)
                        converged = OT_tmp.jobex(gcart = 3)
                        if not converged:
                            print('The geometry optimization with multiplicity %d has failed.' %(M+2))
                            exit()
                        mol_tmp = OT_tmp.get_mol()
                        Mol(mol_tmp).make_molecular_graph()
                        mg_tmp = mol_tmp.graph.molg
                        equal = molsys.addon.graph.is_equal(mg_ini, mg_tmp)[0]
                        if not equal:
                             print('The graph still changes with the multiplicity %d.' %(M+2))
                             exit()
                        else:
                             print('The graph does not change with the multiplicity %d.' %(M+2))
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
            multiplicities.append(M)
            QM_paths.append(QM_path)

        return multiplicities, QM_paths


    def reaction_workflow(self, lot = '', rbonds = [], path_ref_educts = [], path_ref_products = [], path_ref_ts = '', atom_ids_dict = {}):
        ''' This method considers the reaction event and optimizes the species accordingly.
            All of the input variables are retrieved from the RDB database, but should also work if one would like to just provide some reference structures...
            rbonds           : list of integers : The indices of the reactive bonds in the TS structure. The atom indexing is like in python, 0,1,2,...
            path_ref_educts  : list of strings  : The list of paths to the reference educts cartesian coordinate files.
            path_ref_products: list of strings  : The list of paths to the reference products cartesian coordinate files.
            path_ref_ts      : string           : The path to the TS cartesian coordinate files.
        '''
        n_ed   = len(path_ref_educts)
        n_prod = len(path_ref_products)

        multiplicities_ed,   QM_paths_ed   = self.educts_and_products_workflow(path_ref_educts,   lot, 'educt')
        multiplicities_prod, QM_paths_prod = self.educts_and_products_workflow(path_ref_products, lot, 'product')

        unimolecular = False
        if n_ed == 1 and n_prod == 1: unimolecular = True

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

        try_woelfling = False
        converged, QM_path_ts = self.ts_pre_optimization(path_ref_ts, lot, M, rbonds)
        if not converged:
            print('The TS pre-optimization did not converge.')     
            try_woelfling = True
        else:
            found, reason = self.ts_eigenvector_following(QM_path_ts, QM_paths_ed, QM_paths_prod)
            if not found:
                print(reason)
                try_woelfling = True

        if try_woelfling:
            if not unimolecular:
                if n_ed > 1 and n_prod == 1:
                    xyz_ed_complex = self.make_rxn_complex(rbonds, atom_ids_dict, 'educt', n_ed, self.path, path_ref_ts)
                    print(xyz_ed_complex)
                elif n_ed > 1 and n_prod > 1:
                    xyz_ed_complex = self.make_rxn_complex(rbonds, atom_ids_dict, 'educt', n_ed, self.path, path_ref_ts)
                    xyz_prod_complex = self.make_rxn_complex(rbonds, atom_ids_dict, 'product', n_prod, self.path, path_ref_ts)
                    print(xyz_ed_complex)
                    print(xyz_prod_complex)
                elif n_ed == 1 and n_prod > 1:
                    xyz_prod_complex = self.make_rxn_complex(rbonds, atom_ids_dict, 'product', n_prod, self.path, path_ref_ts)
                    print(xyz_prod_complex)
        return

 

    def write_submit_py(self, lot = '', rbonds = [], path_ref_educts = [], path_ref_products = [], path_ref_ts = '', atom_ids_dict = {}):
        ''' In order to write a script which will run the desired routine, in this case, the reaction_workflow. Therefore, it needs to be modified if some other task is needed.
            This can later be used to submit jobs to the queuing system by writing a job submission script. See in the Slurm class the write_submission_script.
        '''
        f_path = os.path.join(self.path,"submit.py")
        f = open(f_path,"a")
        f.write("import os\n")
        f.write("import molsys\n")
        f.write("from molsys.util import turbomole\n\n")
        f.write("lot               = '%s'\n" %lot)
        f.write("rbonds            = %s\n"   %str(rbonds))
        f.write("path_ref_educts   = %s\n"   %str(path_ref_educts))
        f.write("path_ref_products = %s\n"   %str(path_ref_products))
        f.write("path_ref_ts       = '%s'\n" %path_ref_ts)
        f.write("atom_ids_dict     = %s\n" %str(atom_ids_dict))
        f.write("OT = turbomole.OptimizationTools()\n")
        f.write("OT.reaction_workflow(lot, rbonds, path_ref_educts, path_ref_products, path_ref_ts, atom_ids_dict)\n")
        f.close()
        return


#    def write_submit_py(self, TS, M, active_atoms, max_mem, ref_struc_path, lot, unimolecular = False, path_ref_educts = [], path_ref_products = [],  path_ed_complex = '',  path_prod_complex = ''):
#        f_path = os.path.join(self.path,"submit.py")
#        f = open(f_path,"a")
#        f.write("import os\n")
#        f.write("import molsys\n")
#        f.write("from molsys.util import turbomole\n\n")
#        f.write("max_mem           = %d\n" % max_mem)
#        f.write("ref_struc_path    = '%s' \n" %ref_struc_path)
#        f.write("lot           = '%s'\n" % lot)
#        if TS:
#            f.write("M                 = %d\n" %M)
#            f.write("TS             = True\n")
#            f.write("active_atoms      = %s\n" %str(active_atoms))
#            f.write("unimolecular      = %s\n" %str(unimolecular))
#            f.write("path_ref_educts   = %s\n" %str(path_ref_educts))
#            f.write("path_ref_products = %s\n" %str(path_ref_products))
#            f.write("path_ed_complex   = '%s'\n" %path_ed_complex)
#            f.write("path_prod_complex = '%s'\n\n" %path_prod_complex)
#        else:
#            f.write("M                 = None\n")
#            f.write("TS             = False\n")
#            f.write("active_atoms      = []\n\n")
#        f.write("# Optimization Tools\n")
#        f.write("ST = turbomole.OptimizationTools()\n\n")
#        f.write("atom, converged, init_energy = ST.common_workflow(")
#        f.write("TS = TS,\n")
#        f.write("                   M = M,\n")     
#        f.write("        active_atoms = active_atoms,\n")
#        f.write("             max_mem = max_mem,\n")
#        f.write("      ref_struc_path = ref_struc_path,\n")
#        f.write("             lot = lot)\n\n")
#        f.write("if not converged:\n")
#        f.write("    f = open('NOTFOUND','a')\n")
#        f.write("    f.write('ridft did not converge.')\n")
#        f.write("    f.close()\n")
#        f.write("else:\n")
#        if TS:
#            f.write("    found, reason = ST.transition_state_workflow(active_atoms, path_ref_educts, path_ref_products)\n")
#            f.write("    if found:\n")
#            f.write("        os.system('touch FOUND')\n")
#            f.write("    else:\n")
#            f.write("        if unimolecular:\n")
#            f.write("            ST.woelfling_workflow(unimolecular = True,  M = M, active_atoms = active_atoms, max_mem = max_mem, lot = lot, path_ref_ed = path_ref_educts[0], path_ref_prod = path_ref_products[0])\n")
#            f.write("        else:\n")
#            f.write("            ST.woelfling_workflow(unimolecular = False, M = M, active_atoms = active_atoms, max_mem = max_mem, lot = lot, path_ed_complex = path_ed_complex, path_prod_complex = path_prod_complex)\n")
#            f.write("        f = open('NOTFOUND','a')\n")
#
#            f.write("        f.write('%s' %reason)\n")
#            f.write("        f.close()\n")
#        else:
#            f.write("    if not atom:\n")
#            f.write("        found, reason = ST.minima_workflow()\n")
#            f.write("        if found:\n")
#            f.write("            os.system('touch FOUND')\n")
#            f.write("        else:\n")
#            f.write("            f = open('NOTFOUND','a')\n")
#            f.write("            f.write('%s' %reason)\n")
#            f.write("            f.close()\n")
#            f.close()
#        return
#


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





