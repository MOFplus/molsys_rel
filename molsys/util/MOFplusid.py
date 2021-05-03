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

        # get the systrekey for topo
        self.t.addon("topo") 
        self.systrekey = self.t.topo.systrekey
        self.RCSRname  = self.t.topo.RCSRname

        # now test the vertices to make sure that the reduced set in the systrekey spans the system
        # => all vertices of the same systrekey have to have the same BB type
        # NOTE we assume that this then holds for edges, too (maybe need to be implemented as a separate test at some point?)
        reduce_OK = True
        self.sysv_to_bb = {}
        self.nsysv = 0
        for i, bb in enumerate(self.vmap):
            for b in bb:
                sysv = self.t.fragnumbers[b]
                # print ("bb %d  vertex %d is systre_vertex %d" % (i, b, sysv))
                if sysv not in self.sysv_to_bb:
                    self.sysv_to_bb[sysv] = i
                    self.nsysv += 1
                else:
                    if self.sysv_to_bb[sysv] != i:
                        reduce_OK = False
        # what to do if reduction is not ok? => TBI
        assert reduce_OK, "The reduction of vertices by systre is not accepable => implement 2x2x2 supercell!!"
        


