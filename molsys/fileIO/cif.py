import numpy
import string
from . import txyz
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
    for i in range(mol.natoms):
        f.write(" %s  %s %12.6f  %12.6f  %12.6f \n" % (string.upper(mol.elems[i]),string.upper(mol.elems[i]),\
            frac_xyz[i,0],frac_xyz[i,1],frac_xyz[i,2],))
    f.write("  \n")
    f.write("#END  \n")
    f.close()
    return

def read(mol,fname,make_P1=True,detect_conn=True):
    """BUG: cif instance cannot be deepcopied!"""
    """BUG: currently does not always support symmetry operations"""
    try: 
        import CifFile
    except ImportError:
        raise ImportError('pycifrw not installed, install via pip!')
    cf = CifFile.ReadCif(fname)
    if len(cf.keys()) != 1:
        for key in cf.keys(): print(key)
        raise IOError('Cif File has multiple entries ?!')
    cf = cf[cf.keys()[0]]
    cellparams=[]
    cellparams.append(cf['_cell_length_a'])
    
    elems = [str(i) for i in cf.GetItemValue('_atom_site_type_symbol')]
    elems = [i.lower() for i in elems]
    x = [format_float(i) for i in cf.GetItemValue('_atom_site_fract_x')]
    y = [format_float(i) for i in cf.GetItemValue('_atom_site_fract_y')]
    z = [format_float(i) for i in cf.GetItemValue('_atom_site_fract_z')]
    a = format_float(cf.GetItemValue('_cell_length_a'))
    b = format_float(cf.GetItemValue('_cell_length_b'))
    c = format_float(cf.GetItemValue('_cell_length_c'))
    alpha = format_float(cf.GetItemValue('_cell_angle_alpha'))
    beta = format_float(cf.GetItemValue('_cell_angle_beta'))
    gamma = format_float(cf.GetItemValue('_cell_angle_gamma'))
    mol.set_natoms(len(elems))
    mol.set_cellparams([a,b,c,alpha,beta,gamma])
    mol.set_xyz_from_frac(numpy.array([x,y,z]).T)
    #mol.wrap_in_box()
    mol.set_elems(elems)
    mol.set_atypes(['-1']*len(elems))
    mol.set_nofrags()
    mol.set_empty_conn()
    mol.cifdata = cf
    if make_P1: 
        mol.addon('spg')
        mol.proper_cif = mol.spg.make_P1()
    if detect_conn:
        mol.detect_conn()
    return

def format_float(data):
    if data.count('(') != 0:
        data = data.split('(')[0]
    return float(data)
            
