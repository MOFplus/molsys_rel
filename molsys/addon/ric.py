"""
This file implements a ric_fit class.
It is inheriting all features from RedIntCoords in ric.py
and adds features to load a respective reference structure and Hessian.
In addition, a weight matrix is held to evaluate various kinds of weighted mean
square deviations to be used as ingredients to fittness values
"""

import string
import numpy as np
from ff_gen.ric import RedIntCoords 
import copy

ricmapping = {"bnd": "str",
        "ang": "ibe",
        "oop": "obe",
        "dih": "tor"}


class ric(RedIntCoords):
    """
    class to compute redundant internal coordinates (ric)
    by using the inherited ric module and to compute deviations from a corresponding
    reference.
    """
    
    def __init__(self, mol, lindict = {}):
        self._mol = mol
        ### check if ff is already initialized, else do
        if hasattr(self._mol,"ff") == False:
            self._mol.addon("ff")
            self._mol.ff.find_rics()
        elif hasattr(self._mol.ff.ric, "bnd") == False:
            self._mol.ff.ric.find_rics()
        ### init RedIntCoords
        RedIntCoords.__init__(self)
        self.val_dict = {"str": self.get_val_stretches,
                "ibe": self.get_val_in_bends,
                "obe": self.get_val_out_bends,
                "lbe": self.get_val_lin_bends,
                "tor": self.get_val_torsions,
                "eck": self.get_val_eckarts,
                "hes": self.get_ric_hessian}
        self.lindict = lindict
        return

    def setup_rics(self, full = True):
        self._mol.ff.ric.compute_rics()
        ### bonds
        for i,r in enumerate(self._mol.ff.ric.bnd):
            if full:
                self.add_stretch(np.array(list(r))+1)
            elif r.used:
                self.add_stretch(np.array(list(r))+1)
        ### angles
        # TODO implement linear angles
        for i,r in enumerate(self._mol.ff.ric.ang):
            if abs(r.value-180.0) < 2.0:
                r.lin = True                
            else:
                r.lin = False
            if full:
                if r.lin:
                    self.add_lin_bend_mod(r)
                else:
                    self.add_in_bend(np.array(list(r))+1)
            elif r.used:
                if r.lin:
                    self.add_lin_bend_mod(r)
                else:
                    self.add_in_bend(np.array(list(r))+1)
        ### oop
        for i, r in enumerate(self._mol.ff.ric.oop):
            a,b,c,d = r
            if full:
                self.add_out_bend(np.array([a,b,c,d])+1)
                self.add_out_bend(np.array([a,c,b,d])+1)
                self.add_out_bend(np.array([a,d,c,b])+1)
            elif r.used:
                self.add_out_bend(np.array([a,b,c,d])+1)
                self.add_out_bend(np.array([a,c,b,d])+1)
                self.add_out_bend(np.array([a,d,c,b])+1)
        ### dihedrals
        # TODO so far no fragtors are available
        for i, r in enumerate(self._mol.ff.ric.dih):
            if full:
                self.add_torsion(np.array(list(r))+1)
            elif r.used:
                self.add_torsion(np.array(list(r))+1)
        self.add_eckart()
        self._mol.set_real_mass()
        self.setup(masses = np.array(self._mol.get_mass()))
        self.report_rics("rics.dat")
        return
    
    def add_lin_bend_mod(self, indices):
        if indices in self.lindict.keys():
            self.add_lin_bend(np.array(list(indices))+1, ref = self.lindict[indices]+1)
        else:
            if len(self._mol.conn[indices[0]])>1:
                lref = copy.copy(self._mol.conn[indices[0]])
                lref.pop(lref.index(indices[1]))
                self.add_lin_bend(np.array(list(indices))+1, ref = lref[0])
            elif len(self._mol.conn[indices[2]])>1:
                lref = copy.copy(self._mol.conn[indices[2]])
                lref.pop(lref.index(indices[1]))
                self.add_lin_bend(np.array(list(indices))+1, ref = lref[0])
            else:
                raise ValueError("No reference atom found for linear bend %s" % indices)
        

    @property
    def first_str(self):
        return 0

    @property
    def first_ibe(self):
        return self.first_str + self.num_stretch

    @property
    def first_obe(self):
        return self.first_ibe + self.num_in_bend

    @property
    def first_tor(self):
        return self.first_obe + self.num_out_bend

    @property
    def first_lbe(self):
        return self.first_tor + self.num_torsion

    @property
    def all_rics(self):
        return {"str": self._stretches,
                "ibe": self._in_bends,
                "obe": self._out_bends,
                "tor": self.tor2dih(self._torsions),
                "lbe": self._lin_bends}

    @property
    def active_rics(self):
        all    = self.all_rics
        active = []
        for k in ["str", "ibe", "obe", "tor", "lbe"]:
            if len(all[k]) > 0: active.append(k)
        return active

    def map_ric(self, ric_type, ind, reversed = False):
        mapper = {"str": (self._stretches, self.first_str),
                "ibe": (self._in_bends, self.first_ibe),
                "obe": (self._out_bends, self.first_obe),
                "tor": (self._torsions, self.first_tor),
                "lbe": (self._lin_bends, self.first_lbe)
                }
        if (ric_type == "tor") and (len(ind) != 12): ind = self.dih2tor(ind)
        if mapper[ric_type][0].count(list(np.array(list(ind))+1)) != 0:
            i = mapper[ric_type][0].index(list(np.array(list(ind))+1))
            iglob = i + mapper[ric_type][1]
            return i, iglob
        elif reversed == False and ric_type in ["str", "ibe", "tor", "lbe"]:
            return self.map_ric(ric_type, ind[::-1], reversed = True)
        else:
            return None, None
    
    def dih2tor(self,ind):
        assert len(ind) == 4
        mapping = [0,5,6,7]
        new_ind = 12*[-1]
        for i, idx in enumerate(ind): new_ind[mapping[i]] = idx
        return new_ind

    def tor2dih(self, ind):
        if type(ind[0]) == int:
            return np.array(ind)[0,5,6,7].tolist()
        else:
            return np.array(ind)[:,[0,5,6,7]].tolist()

    def report_rics(self, file = None):
        buffer = ""
        count = 0
        rics = self.all_rics
        for k in self.active_rics:
            for i,r in enumerate(rics[k]):
                buffer+= ("%4d %3d %5s " % (count,i,k)+len(r)*" %4d")  % tuple(r)+"\n"
                count += 1
        if file == None:
            print buffer
        else:
            with open(file, "w") as f:
                f.write(buffer)
