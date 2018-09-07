from molsys.util.images import arr2idx, idx2arr
import numpy as np
import copy

### COLOR UTILITIES ###
# color dictionaries based on default molden element colors
ecolor2elem = [
    "b" ,"f" ,"n" ,"o" ,"c" ,"he","ne","ge","li","s" ,"cl","p" ,"al","si",
]
elem2ecolor = dict(ke[::-1] for ke in enumerate(ecolor2elem))
maxecolor = len(ecolor2elem)
vcolor2elem = [
    "n" ,"o" ,"b" ,"f" ,"c" ,"he","ne","ge","li","s" ,"cl","p" ,"si","al",
]
elem2vcolor = dict(kv[::-1] for kv in enumerate(vcolor2elem))
maxvcolor = len(vcolor2elem)

# string conversion tools: elem+atype+color <-> string #
def elematypecolor2string(elem, atype, color):
    """
    return formatted string from element, atomtype, and color
    """
    return "%s_%s/%s" % (elem, atype, color)

eac2str = elematypecolor2string #nickname

def string2elematypecolor(st):
    """
    return element, atomtype, and color from a given formatted string.

    Color is after the last slash
    Element is before the first underscore
    Atype is everything in the middle
    """
    self.assert_eacstr(st)
    colorsplit = st.split("/")
    rest, color = colorsplit[:-1], colorsplit[-1]
    rest = "".join(rest)
    elemsplit = rest.split("_")
    elem, atype = elemsplit[0], elemsplit[1:]
    return elem, atype, color

str2eac = string2elematypecolor #nickname

def assert_eacstr(st):
    """
    check format of element_atomtype/color string
    """
    assert st.index("_"), "string is not well-formatted: no \"_\" found"
    assert st.rindex("/"), "string is not well-formatted: no \"/\" found"
    assert st.index("_") < st.rindex("/"), "string is not well-formatted: first \"_\" after last \"/\""
    return

# make mol objects out of graph elements (edges and/or vertices) #

def make_emol(m, alpha, ecolors=None):
    """
    make mol object out of edge colors
    """
    if ecolors is None:
        vcolors = [0]*len(m.nbonds)
    etab = m.etab
    ralpha = 1./alpha #reverse alpha
    calpha = 1-ralpha #one's complement of reverse alpha
    xyz_a = []
    xyz_c = []
    new_etab = []
    if m.use_pconn:
        for ei,ej,p in etab: ### SELECTION TBI
            xyz_ai = m.xyz[ei]
            xyz_a.append(xyz_ai)
            xyz_ci = m.get_neighb_coords_(ei,ej,idx2arr[p])
            xyz_c.append(xyz_ci)
    else:
        for ei,ej in etab: ### SELECTION TBI
            xyz_ai = m.xyz[ei]
            xyz_a.append(xyz_ai)
            xyz_ci =  m.xyz[ej]
            xyz_c.append(xyz_ci)
    xyz_a = np.array(xyz_a)
    xyz_c = np.array(xyz_c)
    if alpha == 2:
        xyz_c = ralpha*(xyz_a + xyz_c)
        me = m.from_array(xyz_c, use_pconn=m.use_pconn)
    else:
        xyz_c1 = calpha*xyz_a + ralpha*xyz_c
        xyz_c2 = ralpha*xyz_a + calpha*xyz_c
        me = m.from_array(np.vstack([xyz_c1,xyz_c2]), use_pconn=m.use_pconn)
    me.is_topo = True
    if hasattr(m,'cell'):
        me.set_cell(m.cell)
    if hasattr(m,'supercell'):
        me.supercell = m.supercell[:]
    if m.use_pconn:
        me.use_pconn = True
    if alpha == 2:
        me.elems = [ecolor2elem[v] for v in ecolors] # N.B.: no connectivity
        if m.use_pconn:
            pimg = me.get_frac_xyz()//1
            me.xyz -= np.dot(pimg,me.cell)
            for k,(i,j,p) in enumerate(etab):
                newe1 = i,k,arr2idx[pimg[k]]
                newe2 = j,k,arr2idx[idx2arr[p]-pimg[k]]
                new_etab.append(newe1)
                new_etab.append(newe2)
        else:
            new_etab = etab[:]
    else:
        me.elems = [ecolor2elem[v] for v in ecolors*2] # with connectivity
        ctab = [[i,i+me.natoms/2] for i in range(me.natoms/2)]
        me.set_ctab(ctab, conn_flag=True)
        if m.use_pconn:
            pimg = me.get_frac_xyz()//1
            me.xyz -= np.dot(pimg,me.cell)
            ptab = [pimg[i+me.natoms/2]-pimg[i] for i in range(me.natoms/2)]
            me.set_ptab(ptab, pconn_flag=True)
            for k,(i,j,p) in enumerate(etab):
                newe1 = i,k,arr2idx[pimg[k]]
                newe2 = j,k+len(etab),arr2idx[idx2arr[p]-pimg[k+len(etab)]]
                new_etab.append(newe1)
                new_etab.append(newe2)
        else:
            new_etab = ctab[:]
    me.new_etab = new_etab ### MOVE TO make_mol
    return me

def make_vmol(m, vcolors=None):
    """
    make mol object out of graph vertices
    """
    if vcolors is None:
        vcolors = [0]*len(m.natoms)
    mv = copy.deepcopy(m)
    for i in range(mv.natoms):
        mv.atypes[i] = elematypecolor2string(
            mv.elems[i],
            mv.atypes[i],
            vcolors[i]
        )
        mv.elems[i] = vcolor2elem[vcolors[i]]
    if hasattr(m,'cell'): mv.set_cell(m.cell)
    if hasattr(m,'supercell'): mv.supercell = m.supercell[:]
    return mv

def make_mol(m, alpha, ecolors=None, vcolors=None, use_edge=True, use_vertex=True):
    """
    make mol object out of graph elements (edges and/or vertices)
    if both edges and vertices, it takes care of the connectivity too
    """
    if use_edge and use_vertex:
        mm = make_emol(m, alpha, ecolors=ecolors)
        ne = mm.natoms
        mv = make_vmol(m, vcolors=vcolors)
        mm.add_mol(mv) # N.B.: in THIS EXACT ORDER, otherwise KO connectivity
        ### connectivity ###
        ctab = []
        if m.use_pconn:
            ptab = []
        if m.use_pconn:
            for i,j,p in mm.new_etab:
                ctab.append((i+ne,j))
                ptab.append(idx2arr[p])
        else:
            for i,j in mm.new_etab:
                ctab.append((i+ne,j))
        mm.set_ctab(ctab, conn_flag=True)
        if m.use_pconn:
            mm.set_ptab(ptab, pconn_flag=True)
    elif use_edge:
        mm = make_emol(m, alpha, ecolors=ecolors)
    elif use_vertex:
        mm = m.make_vmol(vcolors=vcolors)
    if hasattr(m,'cell'): mm.set_cell(m.cell)
    if hasattr(m,'supercell'): mm.supercell = m.supercell[:]
    return mm

