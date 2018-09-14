import numpy
import string
from . import txyz
from collections import Counter
import logging

logger = logging.getLogger("molsys.io")

def write(mol,fname, name='', write_bonds=True):
    """
    Routine, which writes a cif file in P1
    :Parameters:
        -fname  (str) : name of the cif file
        -mol    (obj) : instance of a molclass
    """
    try:
        f.readline ### do nothing
    except AttributeError:
        raise IOError, "%s is not readable" % f
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
        f.write(" %s%d  %s %12.6f  %12.6f  %12.6f \n" % (mol.elems[i].title(),i,mol.elems[i].title(),\
            frac_xyz[i,0],frac_xyz[i,1],frac_xyz[i,2],))
    if write_bonds:
        f.write("loop_  \n")
        f.write("_geom_bond_atom_site_label_1  \n")
        f.write("_geom_bond_atom_site_label_2  \n")
        mol.set_ctab_from_conn()
        for i,ctab in enumerate(mol.ctab):
            c1,c2=ctab[0],ctab[1]
            e1,e2 =mol.elems[c1].title(), mol.elems[c2].title()
            f.write('%s%d   %s%d \n' % (e1,c1,e2,c2) )
    f.write("  \n")
    f.write("#END  \n")
    f.close()
    return

def read(mol, fname, make_P1=True, detect_conn=True, disorder=None):
    """read CIF file
    :Arguments:
    - make_P1(bool): if True: make P1 unit cell from primitive cell
    - detect_conn(bool): if True: detect connectivity
    - disorder(dict or None): choose disorder group per each disorder assembly
        to consider. Use a dictionary of disorder assembly keys to disorder
        group items, e.g. {"A":"2", "B":"3", "C":2, ...}
        if None: first disordered group in lexical sort is taken per each
            disorder assembly (e.g. {"A":"1", "B":"1", ...}"""
    """BUG: cif instance cannot be deepcopied! (workaround w/ __mildcopy__?"""
    try: 
        import CifFile
    except ImportError:
        raise ImportError('pycifrw not installed, install via pip!')
    cf = CifFile.ReadCif(fname)
    if len(cf.keys()) != 1:
        for key in cf.keys(): print(key)
        raise IOError('Cif File has multiple entries ?!')
    cf = cf[cf.keys()[0]]
    
    try:
        occ = [format_float(i) for i in cf.GetItemValue("_atom_site_occupancy")]
        if any( [i!=1 for i in occ] ):
            disorder_assembly_full = [i for i in cf.GetItemValue("_atom_site_disorder_assembly")]
            disorder_group_full = [i for i in cf.GetItemValue("_atom_site_disorder_group")]
            logger.warning("fractional occupancies in cif file")
            select_disorder = [i for i,e in enumerate(disorder_group_full) if e != '.']
            # remove fully occupied positions (data could be polluted)
            disorder_group = [disorder_group_full[i] for i in select_disorder]
            disorder_assembly = [disorder_assembly_full[i] for i in select_disorder]
            # create disorder list for each disorder assembly
            if disorder is None: # first sorted as disordered
                disorder = {}
                disorder_couple = set(zip(disorder_assembly, disorder_group))
                disorder_dict = {}
                for a,g in disorder_couple:
                    try:
                        disorder_dict[a].append(g)
                    except KeyError:
                        disorder_dict[a] = [g]
                for a in disorder_dict:
                    disorder_dict[a].sort()
                    disorder[a] = disorder_dict[a][0] #take first (default)
            select = [i for i,e in enumerate(disorder_assembly_full) if i in select_disorder if disorder_group_full[i] == disorder[e]]
            select += [i for i,e in enumerate(disorder_group_full) if i not in select_disorder]
    except KeyError as e:
        disorder = None
        # re-raise only if the error message is different than the following
        # otherwise: go forward! no disorder, everything is fine!
        if e.message not in [
            "Itemname _atom_site_occupancy not in datablock",
            "Itemname _atom_site_disorder_assembly not in datablock",
            "Itemname _atom_site_disorder_group not in datablock"
        ]:
            raise(e)

    elems = [str(i) for i in cf.GetItemValue('_atom_site_type_symbol')]
    elems = [i.lower() for i in elems]
    x = [format_float(i) for i in cf.GetItemValue('_atom_site_fract_x')]
    y = [format_float(i) for i in cf.GetItemValue('_atom_site_fract_y')]
    z = [format_float(i) for i in cf.GetItemValue('_atom_site_fract_z')]

    if disorder is not None:
        # select according to given disorder
        elems = [e for i,e in enumerate(elems) if i in select]
        x = [e for i,e in enumerate(x) if i in select]
        y = [e for i,e in enumerate(y) if i in select]
        z = [e for i,e in enumerate(z) if i in select]

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

