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
            self.path = path
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

    def calculate_spin_multiplicity_from(nalpha, nbeta):
        """ Calculate the spin multiplicity from the number of alpha and beta electrons:
            M = 2S+1 = 2*(nalpha-nbeta)*(1/2)+1 = (naplha-nbeta)+1
        """
        M = abs(round(nalpha-nbeta))+1
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


    def make_tmole_dft_input(self, elems, xyz, M, max_mem, title, lot_DFT, fermi = True):
        """Creates a tmole input called 'turbo.in' with c1 symmetry, unrestricted DFT, and RI-J approximation.

        Parameters
        ----------
        elems   : the list of elements
        xyz     : numpy array of shape (len(elems),3)
        M       : the (initial) spin multiplicity
        max_mem : Maximum memory per core to use in the calculations.
        title   : title of the job
        lot_DFT : DFT level of theory, must be string
        fermi   : Boolean for Fermi smearing
        """
        turbo_in_path = os.path.join(self.path,"turbo.in")
        c = xyz*angstrom
        f = open(turbo_in_path,"w")
        f.write("%title\n")
        f.write("%s\n" % title)
        f.write("%method\n")
        f.write("ENRGY :: ri-u%s [gen_mult = %d, gen_symm = c1, scf_dsta = 1, scf_msil = 1000, scf_rico = %d, for_maxc = %d]\n" %
                (lot_DFT, M, max_mem, max_mem))
        f.write("%coord\n")
        for i in range(len(elems)):
            f.write("  %19.14f %19.14f %19.14f   %-2s\n" %
                    (c[i,0],c[i,1], c[i,2], elems[i]))
        f.write("%add_control_commands\n")
        f.write("$disp3\n")
        if fermi:
            f.write("$fermi tmstrt=300.00 tmend= 50.00 tmfac=0.900 hlcrt=1.0E-01 stop=1.0E-03\n")
        f.write("ADD END\n")
        f.write("%end\n")
        f.close()
        return

    def run_tmole(self):
        maindir = os.getcwd()
        os.chdir(self.path)
        os.system("tmole &>/dev/null")
        os.chdir(maindir)
        return

    def check_dscf_converged(self):
        converged = True
        if os.path.isfile(os.path.join(self.path,"dscf_problem")):
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
        maindir = os.getcwd()
        os.chdir(self.path)
        os.system("ridft > ridft.out")
        SPE =  self.get_energy_from_ridft_out()
        os.chdir(maindir)
        return SPE


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
            noise = GeometryTools._make_noise_without_active_atoms(mol, noise_ini, active_atoms) 
        else:
            noise = (numpy.random.rand(mol.natoms,3)-0.5)/10.0
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
            if i in active_atoms:
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


    def detect_and_change_point_group(path=os.getcwd()):
        ### write the define.in file ###
        old_control = os.path.join(path,'control')
        maindir = os.getcwd()
        newsym_path = os.path.join(path,"newsym")

        os.chdir(path)
        os.system('cpc %s' %newsym_path)

        os.chdir(newsym_path)
        ### write a define input ###
        define_in_path = os.path.join(newsym_path,"define.in")
        f = open(define_in_path,"w")
        f.write("\n")
        f.write("y\n")
        f.write("desy \n")
        f.write("ired\n")
        f.write("*\n")
        f.write("\n")
        f.write("use %s\n" %old_control)
        f.write("\n\n\nq")
        f.close()
        
        GeneralTools(newsym_path).invoke_define(define_out_name = "define-test.out")

        os.chdir(path)
        files = os.listdir(path)
        for f in files:
            if os.path.isfile(f) and f not in ['submit.py','submit.out','turbo.in']:
                os.remove(f)

        os.system('cp %s/* %s' %(newsym_path, path))
        shutil.rmtree(newsym_path)
        os.chdir(maindir)
        return

        

class Mol:

    def __init__(self):
        # The dictionary of number of electrons
        self.nelectrons = {'H':1,'He':2, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'S':16, 'Cl':17, 'Ar':18}
        return

    def count_number_of_electrons(self, mol, charge=0):
        """ Counts the number of electrons. """
        nel = 0
        ### count the number of electrons in the system ###
        for t in set(mol.elems):
           amount = mol.elems.count(t)
           nel += self.nelectrons[t.capitalize()]*amount
        ### account for the charge of the molecule ###
        nel -= charge
        return nel

    def make_molecular_graph(self, mol, thresh = 0.2):
        # if the connectivity information not defined before, detect it.
        if not any(mol.conn):
            mol.detect_conn(thresh)
        mol.addon("graph")
        mol.graph.make_graph()
        return

    def separate_molecules(self, mol, thresh = 0.2):
       """Returns a dictionary of mol objects."""
       self.make_molecular_graph(mol)
       mg = mol.graph.molg
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
                   mol_str += '%s %5.10f %5.10f %5.10f' %(mol.elems[j], mol.xyz[j,0], mol.xyz[j,1], mol.xyz[j,2])
                   counter += 1
                   if counter != n_atoms:
                       mol_str += '\n'
           mol_tmp = molsys.mol.from_string(mol_str, 'xyz')
           mol_tmp.detect_conn(thresh)
           mols.append(mol_tmp)
       return mols

 
class OptimizationTools:
    """Methods for the optimization of QM species with ridft."""
    def __init__(self, path=os.getcwd()):
        if not os.path.isdir(path):
            raise FileNotFoundError("The directory %s does not exist." % path)
        else:
            self.path = path
        return

    def freeze_atoms(self, active_atoms):
        a_a = ''
        for i in active_atoms:
            a_a += str(i+1)
            if i != active_atoms[-1]:
               a_a += ',' 
        """ writes a new define file to fix active atoms in internal coordinates. """
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
        os.system('define < %s > %s' %(define_in_path, os.path.join(self.path,'define.out')))
        return

    def ired_and_itvc_1(self):
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


    def check_geo_opt_converged(self):
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
        if ts:
            os.system("jobex -ri -c 150 -trans > jobex.out")
        else:
            os.system("jobex -ri -c 150 > jobex.out")
        return

    def aoforce(self):
        os.system("aoforce > aoforce.out")
        os.remove("dh")
        return

    def IRC(self):
        os.system("DRC -i -c 150 > IRC.out")
        return


    def common_workflow(self, fermi, M, active_atoms, max_mem, ref_struc_path, lot_DFT):
        ''' Performs single point calculation, adds a noise within 0.1 Angstrom to the reference structure. 
        fermi          : logical, to apply Fermi smearing or not
        M              : in case of Fermi smearing, the initial spin multiplicity; otherwise, the multiplicity
        max_mem        : maximum memory per core
        ref_struc_path : the path to the reference structure; e.g. from ReaxFF
        lot_DFT        : string, e.g. tpss/SVP 
        '''
        mol = molsys.mol.from_file(ref_struc_path)
        GT = GeneralTools(self.path)
        atom = False
        new_xyz = GeometryTools.add_noise(mol, active_atoms)
        if mol.natoms == 1:
            atom = True
        if fermi:
            # Determine the orbital occupations through Fermi smearing.
            if M == None:
                # if the multiplicity is not defined already assign the following:
                # M      : initial spin multiplicity
                #        = 3 for systems with even number of electrons
                #        = 2 for systems with odd number of electrons
                nel = Mol().count_number_of_electrons(mol = mol, charge = 0)
                if (nel % 2) == 0:
                    M = 3
                else:
                    M = 2
            GT.make_tmole_dft_input(
                    mol.elems, 
                    new_xyz, 
                    M, 
                    max_mem, 
                    ref_struc_path,
                    lot_DFT,
                    True) # True for Fermi smearing
            GT.run_tmole()
            # Remove the data group $fermi from the control file
            GT.kdg("fermi")
            # If there are partial occupations round them to integers
            GT.round_fractional_occupation()
        else:
            GT.make_tmole_dft_input(
                    mol.elems, 
                    new_xyz, 
                    M, 
                    max_mem, 
                    ref_struc_path,
                    lot_DFT,
                    False) # No Fermi smearing
            GT.run_tmole()
        energy = GT.ridft()
        converged = GT.check_dscf_converged()
        if not converged:
            f = open(os.path.join(self.path,'NOTFOUND'),'a')
            f.write('The ridft calculation did not converge.')
            f.close()
        elif atom:
            os.system('touch %s' %os.path.join(self.path,'FOUND'))
        return atom, converged, energy


    def find_end_points_from_IRC(self, M, max_mem, lot_DFT):
        maindir = os.getcwd()

        # MINUS
        path_minus = os.path.join(self.path, 'displaced_minus') 
        os.system("t2x %s/coord > %s/coord.xyz" %(path_minus, path_minus))
        mol_minus =  molsys.mol.from_file('%s/coord.xyz' %path_minus,'xyz')

        sub_path_minus = os.path.join(path_minus, 'minus')
        os.mkdir(sub_path_minus)
        GeneralTools(sub_path_minus).make_tmole_dft_input(mol_minus.elems, mol_minus.xyz, M, max_mem, 'minus', lot_DFT, False)
        GeneralTools(sub_path_minus).run_tmole()
        os.chdir(sub_path_minus)
        self.jobex()
        converged = self.check_geo_opt_converged()
        # TODO I should check if these are converged... 
        os.system("t2x %s/coord > %s/coord.xyz" %(sub_path_minus, sub_path_minus))
        opt_mol_minus =  molsys.mol.from_file('%s/coord.xyz' %sub_path_minus,'xyz')
        opt_mol_minus.detect_conn(thresh = 0.2)

        mols_minus = Mol().separate_molecules(opt_mol_minus)

        # PLUS
        path_plus  = os.path.join(self.path, 'displaced_plus')
        os.system("t2x %s/coord > %s/coord.xyz" %(path_plus, path_plus))
        mol_plus =  molsys.mol.from_file('%s/coord.xyz' %sub_path_plus,'xyz')

        sub_path_plus = os.path.join(path_plus, 'plus')
        os.mkdir(sub_path_plus)
        GeneralTools(sub_path_plus).make_tmole_dft_input(mol_plus.elems, mol_plus.xyz, M, max_mem, 'plus', lot_DFT, False)
        GeneralTools(sub_path_plus).run_tmole()
        os.chdir(sub_path_plus)
        self.jobex()
        converged = self.check_geo_opt_converged()
        # TODO I should check if these are converged... 
        os.system("t2x %s/coord > %s/coord.xyz" %(sub_path_plus, sub_path_plus))
        opt_mol_plus =  molsys.mol.from_file('%s/coord.xyz' %sub_path_plus,'xyz')
        opt_mol_plus.detect_conn(thresh = 0.2)

        mols_plus = Mol().separate_molecules(opt_mol_plus)

        os.chdir(maindir)
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
                Mol().make_molecular_graph(mol_ed)
                mg_ed = mol_ed.graph.molg
                print(mg_ed.list_properties)

                educt_minus_tmp = False
                for mol_minus in mols_minus:
                    Mol().make_molecular_graph(mol_minus)
                    mg_minus = mol_minus.graph.molg
                    is_equal = molsys.addon.graph.is_equal(mg_ed, mg_minus, use_fast_check=False)[0]
                    if is_equal: 
                        educt_minus_tmp = educt_minus_tmp or True
                    else:
                        educt_minus_tmp = educt_minus_tmp or False
                educt_minus_is_similar = educt_minus_is_similar and educt_minus_tmp

                educt_plus_tmp  = False
                for mol_plus in mols_plus:
                    Mol().make_molecular_graph(mol_plus)
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
                Mol().make_molecular_graph(mol_prod)
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

    def woelfling_workflow(self, unimolecular, M, max_mem, lot_DFT, path_ref_ed = '', path_ref_prod = '', mpfx_ed_complex = '', mpfx_prod_complex = ''):
        '''workflow to perform a TS search using woelfling
        '''
        maindir = os.getcwd()
        woefling_dir = os.path.join(maindir,'woelfling')
        os.mkdir(woefling_dir)
        os.chdir(woefling_dir)
        if not unimolecular:
            # optimize the ReaxFF complex structures
            path_ref_ed   = os.path.join(woefling_dir, 'educt_complex')
            path_ref_prod = os.path.join(woefling_dir, 'product_complex')
            converged_ed,   reason_ed   = self.optimize_rxn_complex(path_ref_ed,   mpfx_ed_complex,   M, max_mem, lot_DFT, active_atoms)
            converged_prod, reason_prod = self.optimize_rxn_complex(path_ref_prod, mpfx_prod_complex, M, max_mem, lot_DFT, active_atoms)
        os.chdir(path_ref_prod)
        os.system('cpc %s' % woefling_dir)
        os.chdir(woelfling_dir) 
        os.system('cp %s %s' %(os.path.join(path_ref_ed,   'coord'), os.path.join(woefling_dir, 'ini')))
        os.system('cp %s %s' %(os.path.join(path_ref_prod, 'coord'), os.path.join(woefling_dir, 'fin')))
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
        max_struc = self.getthemaxenergystruc(path, True)
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
        os.chdir(maindir)
        return

    def optimize_rxn_complex(path, mpfx_complex, M, max_mem, lot_DFT, active_atoms):
        ''' optimizes the reaction complex by constraining internal coordinates of the active atoms.
        path         : string  : the path to where the optimization should be made
        mfpx_complex : string  : the path to the mpfx structure of the complex
        M            : integer : the multiplicity as determined based on that of educts and products
        '''
        converged = False
        reason = ''
        maindir = os.getcwd()
        os.mkdir(path)
        os.chdir(path)
        mol = molsys.mol.from_file(mfpx_complex)
        GT = GeneralTools(path)
        GT.make_tmole_dft_input(elems = mol.elems, xyz = mol.xyz, M = M, max_mem = max_mem, title = mfpx_ed_complex, lot_DFT = lot_DFT, fermi = False)
        GT.run_tmole()
        converged = GT.check_dscf_converged()
        if not converged:
            reason = 'The self consistent field calculation did not converge for the educt complex.'
        else:
            self.freeze_atoms(active_atoms)
            self.jobex()
            converged = self.check_geo_opt_converged()
            if not converged:
                reason = 'The geometry optimization did not converge for the educt complex.' 
            else:
                os.system('kdg intdef')
                os.system('kdg redundant')
                self.ired()
                converged = True
        os.chdir(maindir)
        return converged, reason


    def transition_state_workflow(self, active_atoms, path_ref_educts, path_ref_products):
        found = False
        reason = ""    
        n_ed = len(path_ref_educts)
        n_prod = len(path_ref_products)
        # freeze the active atoms
        self.freeze_atoms(active_atoms)
        # pre-optimization
        self.jobex()
        if not self.check_geo_opt_converged():
            reason += 'Transition state pre-optimization did not converge.\n'
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
                GeometryTools.detect_and_change_point_group(self.path)

            # perform aoforce calculation     
            self.aoforce()
            
            # check the number of imaginary frequencies
            inum, imfreq = self.check_imaginary_frequency()
            if inum == 0:        
                reason += 'No imaginary frequency at the start structure.\n'
                print(reason)
            elif inum == 1:
                self.jobex(ts=True)
                if not self.check_geo_opt_converged():
                   reason += 'The transition state optimzation did not converge.\n'
                   print(reason)
                else:
                   self.aoforce()
                   inum, imfreq = self.check_imaginary_frequency()
                   if inum == 1: 
                      print('There is only one imaginary frequency. The intrinsic reaction coordinate is being calculated.\n')
                      self.IRC()
                      mols_minus, mols_plus = self.find_end_points_from_IRC()
                      found, reason = self.check_end_points(mols_minus, mols_plus, path_ref_educts, path_ref_products)
                   else:
                      reason += 'The final number of imaginary frequency is not 1.\n'
            elif inum > 1:
                reason += 'There are more than one imaginary frequencies. But it will try to optimize.\n'
                print(reason)
                self.jobex(ts=True)
                if not self.check_geo_opt_converged():
                   reason += 'The transition state optimzation did not converge.\n'
                   print(reason)
                else:
                   self.aoforce()
                   inum, imfreq = self.check_imaginary_frequency()
                   if inum == 1:
                      print('There is only one imaginary frequency. The intrinsic reaction coordinate is being calculated.\n')
                      self.IRC()
                      mols_minus, mols_plus = self.find_end_points_from_IRC()
                      found, reason = self.check_end_points(mols_minus, mols_plus, path_ref_educts, path_ref_products)
                   else:
                      reason += 'The final number of imaginary frequency is not 1.\n'
        return found, reason

    def minima_workflow(self):
        found = False
        reason = ""

        point_group_initial = GeometryTools.get_point_group_from_coord(self.path,'coord')
        print("point_group_initial", point_group_initial)
        self.jobex()
        if self.check_geo_opt_converged():
            point_group_final = GeometryTools.get_point_group_from_coord(self.path,'coord')
            print("point_group_final", point_group_final)
            if point_group_final != "c1":
                GeometryTools.detect_and_change_point_group(self.path)
                GeneralTools(self.path).ridft()
                self.jobex()
                if not self.check_geo_opt_converged():
                    reason += "The geometry optimization is failed.\n"
                else:
                    self.aoforce()
            else:
                self.aoforce()
            if os.path.isfile('aoforce.out'):
                inum, imfreq = self.check_imaginary_frequency()
                if inum == 0:
                    found = True
                    print("The equilibrium structure is found succesfully!")
                else:
                    reason += "There are imaginary frequencies. This is not an equilibrium structure.\n"
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

    def write_submit_py(self, TS, M, active_atoms, max_mem, ref_struc_path, lot_DFT, unimolecular = False, path_ref_educts = [], path_ref_products = [],  mfpx_ed_complex = '',  mfpx_prod_complex = ''):
        f_path = os.path.join(self.path,"submit.py")
        f = open(f_path,"a")
        f.write("import os\n")
        f.write("import molsys\n")
        f.write("from molsys.util import turbomole\n\n")
        f.write("max_mem           = %d\n" % max_mem)
        f.write("ref_struc_path    = '%s' \n" %ref_struc_path)
        f.write("lot_DFT           = '%s'\n" % lot_DFT)
        if TS:
            f.write("M                 = %d\n" %M)
            f.write("fermi             = False\n")
            f.write("active_atoms      = %s\n" %str(active_atoms))
            f.write("unimolecular      = %s\n" %str(unimolecular))
            f.write("path_ref_educts   = %s\n" %str(path_ref_educts))
            f.write("path_ref_products = %s\n" %str(path_ref_products))
            f.write("mfpx_ed_complex   = '%s'\n" %mfpx_ed_complex)
            f.write("mfpx_prod_complex = '%s'\n\n" %mfpx_prod_complex)
        else:
            f.write("M                 = None\n")
            f.write("fermi             = True\n")
            f.write("active_atoms      = []\n\n")
        f.write("# Optimization Tools\n")
        f.write("ST = turbomole.OptimizationTools()\n\n")
        f.write("atom, converged, init_energy = ST.common_workflow(")
        f.write("fermi = fermi,\n")
        f.write("                   M = M,\n")     
        f.write("        active_atoms = active_atoms,\n")
        f.write("             max_mem = max_mem,\n")
        f.write("      ref_struc_path = ref_struc_path,\n")
        f.write("             lot_DFT = lot_DFT)\n\n")
        f.write("if not converged:\n")
        f.write("    f = open('NOTFOUND','a')\n")
        f.write("    f.write('ridft did not converge.')\n")
        f.write("    f.close()\n")
        f.write("else:\n")
        if TS:
            f.write("    found, reason = ST.transition_state_workflow(active_atoms, path_ref_educts, path_ref_products)\n")
            f.write("    if found:\n")
            f.write("        os.system('touch FOUND')\n")
            f.write("    else:\n")
            f.write("        if unimolecular:\n")
            f.write("            ST.woelfling_workflow(unimolecular = True,  M = M, max_mem = max_mem, lot_DFT = lot_DFT, path_ref_ed = path_ref_educts[0], path_ref_prod = path_ref_products[0])\n")
            f.write("        else:\n")
            f.write("            ST.woelfling_workflow(unimolecular = False, M = M, max_mem = max_mem, lot_DFT = lot_DFT, mpfx_ed_complex = mfpx_ed_complex, mpfx_prod_complex = mpfx_prod_complex)\n")
            f.write("        f = open('NOTFOUND','a')\n")

            f.write("        f.write('%s' %reason)\n")
            f.write("        f.close()\n")
        else:
            f.write("    if not atom:\n")
            f.write("        found, reason = ST.minima_workflow()\n")
            f.write("        if found:\n")
            f.write("            os.system('touch FOUND')\n")
            f.write("        else:\n")
            f.write("            f = open('NOTFOUND','a')\n")
            f.write("            f.write('%s' %reason)\n")
            f.write("            f.close()\n")
            f.close()
        return



class Slurm:

    def get_avail_nodes():
        sinfo = os.popen('sinfo --format="%n %t %c %m"').read()
        n_t_c_m = sinfo.split("\n")
        avail_nodes = []
        for lines in n_t_c_m:
            line = lines.split()
            if len(line)==4 and line[1] == "idle":
                avail_nodes.append((line[0],int(line[2]),int(line[3])))
        return avail_nodes


    def get_avail_nodes_of_(partition="normal"):
        sinfo = os.popen('sinfo --format="%n %t %P"').read()
        n_t_P = sinfo.split("\n")
        avail_nodes = []
        for lines in n_t_P:
            line = lines.split()
            if len(line)==3 and line[1] == "idle" and line[2].startswith(partition):
                avail_nodes.append(line[0])
        return avail_nodes

    def get_partition_info(partition="normal"):
        sinfo = os.popen('sinfo --format="%P %c %m"').read()
        P_c_m = sinfo.split("\n")
        for lines in P_c_m:
            line = lines.split()
            if len(line)==3 and line[0].startswith(partition):
                CPUS = int(line[1])
                MEM  = int(line[2])
        return CPUS, MEM

    def write_submission_script(path=os.getcwd(), TURBODIR="", ntasks=8, partition="normal"):
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





