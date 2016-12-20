# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 18:29:19 2016

@author: rochus

          Fragmentizer class

          depends on graph addon (this means graph_tool must be installed)

"""

import os
import numpy
import logging
import glob
import molsys

import logging

logger = logging.getLogger("molsys.fragmentizer")

class fragmentizer:

    def __init__(self):
        """
        fragmentizer loads a catalog of fragments with ending .frag from either
        the current directory or from $MOLSYS_FRAGS

        """
        if os.environ.has_key("MOLSYS_FRAGS"):
            self.frag_path = os.environ["MOLSYS_FRAGS"]
        else:
            self.frag_path = "."
        self.read_frags()
        return
        
    def read_frags(self):
        logger.info("Fragmentizer reading frags from %s" % self.frag_path)
        frag_files = glob.glob(self.frag_path + "/*.mfpb")
        self.fragments = {}
        frag_natoms,fragnames=[],[]
        for f in frag_files:
            m = molsys.mol()
            m.read(f, ftype="mfpx")
            frag_natoms.append(m.natoms)
            m.addon("graph")
            m.graph.make_graph()
            # m.graph.plot_graph(f)
            logger.info("read %s" % f)
            fragname = f.split("/")[-1].split(".")[0]
            fragnames.append(fragname)
            self.fragments[fragname] = m
        ### okay, we want the large fragments to be tested first
        self.frag_order = [fragnames[i] for i in numpy.argsort(frag_natoms)[::-1]]#[::-1]]
        print self.frag_order 

    def __call__(self, mol):
        """
        tries to assign all fragmnets in the catalog to the mol object

        :Parameters:

            - mol : mol object to be fragmentized

        """
        mol.addon("graph")
        mol.graph.make_graph()
        mol.set_nofrags()
        # mol.graph.plot_graph("mol")
        fi = 0
        for f in self.frag_order:
            fidx = mol.graph.find_fragment(self.fragments[f],add_hydrogen=False)
            for flist in fidx:
                exists=False
                for i in flist:
                    if mol.fragtypes[i] != '-1':
#                        print mol.fragtypes[i],i,flist,f
                        exists=True
                        break
                    mol.fragtypes[i]   = f
                    mol.fragnumbers[i] = fi
                if exists==False: fi += 1
                    
        # print mol.fragtypes
        # print mol.fragnumbers
        return 
