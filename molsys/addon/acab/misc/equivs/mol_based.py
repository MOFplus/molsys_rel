from molsys.util.color import make_mol, make_emol, make_vmol
from itertools import permutations
from collections import Counter

def cheap_check(_acab, cola, colb):
    ma = make_structure(_acab, ecolors=cola[0], vcolors=cola[1], alpha=_acab.alpha)
    mb = make_structure(_acab, ecolors=colb[0], vcolors=colb[1], alpha=_acab.alpha)
    ea = Counter(ma.elems)
    eb = Counter(mb.elems)
    if sorted(ea.keys()) != sorted(eb.keys()):
        return False
    else:
        elems = tuple(ea.keys())
    perms = permutations(elems)
    for perm in perms:
        if elems != perm:
            scrambled = dict(zip(elems,perm))
            scrambled_elems = [scrambled[i] for i in ma.elems]
            if mb.elems == scrambled_elems:
                return True
    else:
        return False

def expensive_check(_acab, cola, colb):
    raise NotImplementedError
    if not cheap_check(_acab, cola, colb):
        return False
    if True:
        return True
    else:
        return False

def make_structure(_acab, ecolors=None, vcolors=None, alpha=2):
    if ecolors and vcolors:
        m = make_mol(_acab._mol, alpha, ecolors=ecolors, vcolors=vcolors,
            use_vertex=_acab.use_vertex, use_edge=_acab.use_edge)
    elif ecolors:
        m = make_emol(_acab._mol, alpha, ecolors=ecolors)
    elif vcolors:
        m = make_vmol(_acab._mol, vcolors=vcolors)
    else:
        logger.error("step solutions are not constraint: infinite loop!")
        raise TypeError("unconstraint step solutions: infinite loop!")
    return m

