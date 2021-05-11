"""MOF+id

function to generate a MOF+id from a MOF mol object

should this be a class?

"""





class MOFplusid:

    def __init__(self, mol):
        self.mol = mol

        self.mol.addon("graph")

        # decompose the MOF into its BBs and underlying topology
        self.t, self.v, self.vmap, self.e, self.emap = self.mol.graph.decompose()
        # NOTE: decompose returns the mol objects of topo, vertex and edges (unique) and their mapping
        #       but after calling, the graph class contains a number of further important info
        #       one of the most important is decomp_vbb_cons which contains the vertices connected to in the 
        #       order as the connectors appear in the unique returned BB

        # get the systrekey for topo
        self.t.addon("topo") 
        self.systrekey = self.t.topo.systrekey
        self.RCSRname  = self.t.topo.RCSRname

        # now test the vertices to make sure that the reduced set in the systrekey spans the system
        # => all vertices of the same systrekey have to have the same BB type
        # NOTE we assume that this then holds for edges, too (maybe need to be implemented as a separate test at some point?)
        reduce_OK = True
        # self.sysv_to_bb = self.systrekey.
        self.nsysv = 0
        for i, vbb in enumerate(self.vmap):
            for bb in vbb:
                sysv = self.t.fragnumbers[bb]
                if sysv not in self.sysv_to_bb:
                    self.sysv_to_bb[sysv] = i
                    self.nsysv += 1
                else:
                    if self.sysv_to_bb[sysv] != i:
                        reduce_OK = False
        # what to do if reduction is not ok? => TBI
        assert reduce_OK, "The reduction of vertices by systre is not accepable => implement 2x2x2 supercell!!"

        # now let us determine the edges each vertex is connected to


        


