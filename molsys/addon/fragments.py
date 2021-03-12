# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 21:30:26 2016

@author: rochus
"""
import string
import logging
logger = logging.getLogger("molsys.fragments")


class fragments:
    """
    fragments is an addon class to support advanced fragment
    handling
    """


    def __init__(self, mol):
        self._mol = mol
        self.setup = False
        self.check()
        if self.setup:
            self.make_frag_conn()
        return

    def check(self):
        """
        check if all atoms are in a fragment and all is consistent
        """
        self.fragnames = []        # this a list of all existing fragment names
        self.frag_atoms = []       # a list of the fragments with their atoms
        self.nfrags =  max(self._mol.fragnumbers)+1
        self.fraglist  = [None]*(self.nfrags) # additional list to hold the fragments with their name
        self.frag_conn = []        # fragment connectivity (indices of fragments)
        self.frag_conn_atoms = []  # atoms making the fragemnt connectivity (tuples of atom indices: first in frgmnt, second in other fragmnt)
        for i in range(self.nfrags):
            self.frag_atoms.append([])
        if self.nfrags == 0:
            self.setup = False
            return
        self.setup = True
        for i in range(self._mol.natoms):
            ft = self._mol.fragtypes[i]
            fn = self._mol.fragnumbers[i]
            if ft == "0":
                logger.error("atom %d (%s) is not in fragment" % (i, self._mol.atypes[i]))
                self.setup=False
            else:
                if self.fraglist[fn] is None:
                    # set name of fragment
                    self.fraglist[fn] = ft
                else:
                    # check if this is the same name
                    assert self.fraglist[fn] == ft, \
                         "The fragmentname %s of atom %d does not match with the prior definition %s" % (ft, i, self.fraglist[fn])
                if ft not in self.fragnames:
                    self.fragnames.append(ft)
                self.frag_atoms[self._mol.fragnumbers[i]].append(i)
        # in the end make sure that all fragments have been named
        if None in self.fraglist:
            raise ValueError("A fragment name is missing")
        return

    def get_occurence_of_frag(self,name):
        '''
            returns the fragment count of the fragment with the given name
        '''
        return self.fraglist.count(name)

    def get_fragnames(self):
        return self.fragnames


    def make_frag_conn(self):
        """
        generate a fragment connectivity
        """
        assert self.setup
        # prepare the atom list for the fragments
        for i in range(self.nfrags):
            self.frag_conn.append([])
            self.frag_conn_atoms.append([])
        for i,f in enumerate(self.fraglist):
            # determine all external bonds of this fragment
            for ia in self.frag_atoms[i]:
                for ja in self._mol.conn[ia]:
                    j = self._mol.fragnumbers[ja]
                    if i != j:
                        # this is an external bond
                        self.frag_conn[i].append(j)
                        self.frag_conn_atoms[i].append((ia,ja))
                        logger.debug("fragment %d (%s) bonds to fragment %d (%s) %s - %s" %\
                                      (i, f, j, self.fraglist[j], self._mol.atypes[ia], self._mol.atypes[ja]))
        return

    # TODO(RS): put the "validaor in here" ... no extra class. source is like in fragmentizer either "file" or "mofp"
    #           remove the reset_atypes stuff. we do not need this anymore
    def analyze_frag_conn(self, validator, reset_atypes=False):
        """
        detect all types of inter-fragment bonds out of frag_conn and frag_conn_atoms

        :Parameter:

            - validator: a validator obejct that validates a frag connection
            - reset_atypes (boolean): if True the fragment name (or the merged fragment name) will be appended to the atype

        """
        self.rev_fraglist = [None]*(self.nfrags)
        self.frag_bond_types= {}
        errors = 0
        # iterate over all frags and check their frag connections (if conencted frag is higher index)
        for i in range(self.nfrags):
            for nj, j in enumerate(self.frag_conn[i]):
                if j>i:
                    atom_pair = self.frag_conn_atoms[i][nj]
                    fbond = [self._mol.atypes[atom_pair[0]]+"_"+self.fraglist[i], self._mol.atypes[atom_pair[1]]+"_"+self.fraglist[j]]
                    fbond.sort()
                    fbond = ":".join(fbond)
                    if not fbond in self.frag_bond_types.keys():
                        self.frag_bond_types[fbond] = ""
                    # now check if this bond is allowed and if we need to change the fragtype in rev_fraglist
                    response = validator(fbond)
                    if response == False:
                        errors += 1
                        logger.error("No fragment connection for %s" % fbond)
                    else:
                        if response != "":
                            self.rev_fraglist[i] = response
                            self.rev_fraglist[j] = response
        if errors == 0:
            if reset_atypes:
                for i in range(self._mol.natoms):
                    f = self._mol.fragnumbers[i]
                    if self.rev_fraglist[f] is None:
                        ft = self.fraglist[f]
                    else:
                        ft = self.rev_fraglist[f]
                    self._mol.atypes[i] += "_"+ft
        return

    def make_frag_graph(self):
        """
        generate a graph of the frag_conn in analogy to the graph addon on the molecular level
        using the graph addons util_graph method
        """
        self._mol.addon("graph")
        # create here a second list of vertex types, with the aryl substituted species
        # for example a naph fragment is substituted by a ph
        phenyl_like = ["naph"]
        vtypes2 = []
        for i,t in enumerate(self.fraglist):
            if t in phenyl_like:
                vtypes2.append("ph")
            else:
                vtypes2.append(t)
        self.frag_graph = self._mol.graph.util_graph(self.fraglist, self.frag_conn, vtypes2=vtypes2)
        # DEBUG here just for debug reasons
        #self._mol.graph.plot_graph("frag_conn", g=self.frag_graph)
        return self.frag_graph
        
    def plot_frag_graph(self, fname, **kwargs):
        self._mol.graph.plot_graph(fname, g=self.frag_graph, **kwargs)
        return

    def upgrade(self, se, rep, rep_n = None):
        """
        upgrades the vertex labels in a frag graph
        :Parameters:
            - se  (str): vertex label to be replaced
            - rep (str): new vertex label
            - rep_n (int): optional. the amount of vertices to be replaced
        """
        assert type(se) == type(rep) == str
        assert hasattr(self, "frag_graph")
        nreplaced = 0
        for v in self.frag_graph.vertices():
            if self.frag_graph.vp.type[v] == se: 
                self.frag_graph.vp.type[v] = rep
                nreplaced += 1
                if rep_n is not None:
                    if nreplaced == rep_n:
                        return
        return


    def frags2atoms(self, frags):
        # TODO improve speed
        #assert type(frags) == list
        idx = []
        for i in range(self._mol.natoms):
            if self._mol.fragnumbers[i] in frags:
                idx.append(i)
        return idx





