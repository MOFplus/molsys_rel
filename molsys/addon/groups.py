import logging
import re  # regexp python module
logger = logging.getLogger("molsys.graph")

""" Group Addon

The Group addon is meant for the identification of sub units
    * get sub units from atom types
    * get sub units from mol object (which stores the sub-unit)
"""


class groups(object):

    def __init__(self, mol):
        """init groups addon

        for the moment simply pass the mol instance the addon was loaded

        Args:
            mol (molsys.mol): mol object
        """
        self._mol = mol
        logger.debug("generated the groups addon")
        self._groups = {}
        return

    def __call__(self, groupname, occurences=None):
        return self._groups[groupname].get_occurences(idx=occurences)

    def register_group(self, name, group_def):
        pass

    def register_group_atype(self, name, atypes=[''], ends_at=[]):
        self._groups.update({name: group(name)})
        g = self._groups[name]
        g.atypes = atypes
        g.ends_at = ends_at
        g.mode = 'atype'
        self.find_group_occurence_atype(g)
        return

    def find_group_occurence_atype(self, g):
        assert g.mode == 'atype'
        # get the list of atoms corresponding to the atype constraints
        atype_idxs = [i for i, e in enumerate(self._mol.get_atypes()) if (any([e.count(
            y) > 0 for y in g.atypes]) and not any([e.count(y) > 0 for y in g.ends_at]))]
        done = False
        while not done:
            idx = atype_idxs[0]
            group_idxs = self.walk_bond(
                self._mol, idx, inds=[idx], continue_at=g.atypes, stop_at=g.ends_at)
            # check if all atypes are contained in the current group
            group_atypes = list(set([self._mol.atypes[x] for x in group_idxs]))

            #print idx, group_idxs
            # add the group only if all atypes are contained in the group_atypes
            #print [[x.count(y) != 0 for x in group_atypes] for y in g.atypes], [any([x.count(y) != 0 for x in group_atypes]) for y in g.atypes]
            if all([any([x.count(y) != 0 for x in group_atypes]) for y in g.atypes]):
                g.occurences.append(group_idxs)
                #if len(group_idxs) > 0: g.occurences.append(group_idxs)
            try:
                for x in group_idxs:
                    atype_idxs.remove(x)
            except:
                pass
            if len(atype_idxs) == 0:
                done = True
        return

    def walk_bond(self, m, start_ind, inds=[], continue_at=[], stop_at=[],stop_at_idx=[]):
        """get group by walking over atoms using the connectivity stopping at certain atom types

        [description]

        Args:
            m ([type]): [description]
            start_ind ([type]): [description]
            inds (list, optional): Defaults to []. [description]
            continue_at (list, optional): Defaults to []. [description]
            stop_at (list, optional): Defaults to []. [description]

        Returns:
            [type]: [description]
        """
        for i, c in enumerate(m.conn[start_ind]):
            # check if proper atype
            if not any([m.atypes[c].count(x) > 0 for x in continue_at]):
                continue
            if any([m.atypes[c].count(x) != 0 for x in stop_at]):
                continue
            if stop_at_idx.count(c) != 0:
                continue
            if inds.count(c) == 0:
                inds.append(c)
                inds = self.walk_bond(
                    m, c, inds=inds, continue_at=continue_at, stop_at=stop_at,stop_at_idx=stop_at_idx)
            else:
                pass
        return inds

    def __repr__(self):
        s = '%d registered groups\n' % (len(self._groups.keys()),)
        for i, k in enumerate(self._groups.keys()):
            s += 'group %d (%20s) : %d occurences\n' % (i, k, self._groups[k].get_noccurences())
            s += self._groups[k].__repr__()
        return s

    def get_group_names(self):
        return self._groups.keys()


class group(object):
    def __init__(self, name):
        self.name = name
        self.mode = None
        self.occurences = []
        return

    def __repr__(self):
        s = ''
        if self.mode == 'atype':
            s += 'containing atypes: ' + \
                len(self.atypes)*' %12s ' % tuple(self.atypes) + '\n'
            s += 'restricted atypes: ' + \
                len(self.ends_at)*' %12s ' % tuple(self.ends_at) + '\n'
        return s

    def get_occurences(self, idx=None):
        if idx is None:
            return self.occurences
        elif hasattr(idx, '__iter__'):
            return [self.occurences[x] for x in idx]
        else:
            print('print provide valid idx')
            return []
    
    def get_noccurences(self):
        return len(self.occurences)
