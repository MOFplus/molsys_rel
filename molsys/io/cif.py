import numpy
import string
import txyz
import logging

def write(mol,fname, name=''):
    """
    Routine, which writes a cif file in P1
    :Parameters:
        -fname  (str) : name of the cif file
        -mol    (obj) : instance of a molclass
    """
    f = open(fname, 'w')
    f.write("data_mofplus.org:%s\n" % name)
    f.write("_symmetry_cell_setting           triclinic \n")
    f.write("_symmetry_space_group_name_H-M   'P 1' \n")
    f.write("_symmetry_Int_Tables_number      1 \n")
    f.write("loop_ \n")
    f.write("_symmetry_equiv_pos_site_id \n")
    f.write("_symmetry_equiv_pos_as_xyz \n")
    f.write("1 x,y,z \n")
    f.write("_cell_length_a          %12.6f        \n" % (mol.cellparams[0]))
    f.write("_cell_length_b          %12.6f        \n" % (mol.cellparams[1]))
    f.write("_cell_length_c          %12.6f        \n" % (mol.cellparams[2]))

    f.write("_cell_angle_alpha       %12.6f        \n" % (mol.cellparams[3]))
    f.write("_cell_angle_beta        %12.6f        \n" % (mol.cellparams[4]))
    f.write("_cell_angle_gamma       %12.6f        \n" % (mol.cellparams[5]))

    f.write("loop_  \n")
    f.write("_atom_site_label   \n")
    f.write("_atom_site_type_symbol  \n")
    f.write("_atom_site_fract_x  \n")
    f.write("_atom_site_fract_y  \n")
    f.write("_atom_site_fract_z \n")
    mol.wrap_in_box()
    frac_xyz = mol.get_frac_xyz()
    for i in xrange(mol.natoms):
	f.write(" %s  %s %12.6f  %12.6f  %12.6f \n" % (string.upper(mol.elems[i]),string.upper(mol.elems[i]),
            frac_xyz[i,0],frac_xyz[i,1],frac_xyz[i,2],))

    f.write("  \n")
    f.write("#END  \n")
    f.close()
    return
