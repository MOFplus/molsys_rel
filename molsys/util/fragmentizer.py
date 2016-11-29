# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 18:29:19 2016

@author: rochus

          Fragmentizer class

          depends on graph addon (this means graph_tool must be installed)

"""

import os
import logging
import glob
import molsys

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
        logging.info("Fragmentizer reading frags from %s" % self.frag_path)
        frag_files = glob.glob(self.frag_path + "/*.mfpb")
        self.fragments = {}
        for f in frag_files:
            m = molsys.mol()
            m.read(f, ftype="mfpx")
            m.addon("graph")
            m.graph.make_graph()
            logging.info("read %s" % f)
            fragname = f.split("/")[-1].split(".")[0]
            print fragname
            self.fragments[fragname] = m
        return

    def __call__(self, mol):
        """
        tries to assign all fragmnets in the catalog to the mol object

        :Parameters:

            - mol : mol object to be fragmentized

        """
        mol.addon("graph")
        mol.graph.make_graph()
        fi = 0
        for f in self.fragments.keys():
            fidx = mol.graph.find_fragment(self.fragments[f])
            for flist in fidx:
                for i in flist:
                    mol.fragtypes[i]   = f
                    mol.fragnumbers[i] = fi
                fi += 1
        print mol.fragtypes
        print mol.fragnumbers
        return