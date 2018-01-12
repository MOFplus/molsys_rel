.. molsys documentation master file, created by
   sphinx-quickstart on Mon Aug 21 14:29:21 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MM Program Interfaces
#####################

Introduction
============
As already explained it is possible to assign force field parameters to a molecular system of
interest by using the FF addon. In addition MOLSYS features interfaces to several molecular
mechanis programms which makes it easy to run molecular dynamics simulations of the system
at the basis of the assigned FF. In this chapter it is explained how to use these interfaces.

Pydlpoly Interface
==================
The Pydlpoly interface is implemented in the ff2pydpoly module. By its help a Pydlpoly instance
for the system of interest can be setup easily. Of course it is required, that Pydlpoly is
installed. In the following example the setup is demonstrated using HKUST-1 as system of
interest.


.. code-block:: python
    :linenos:

    >> import pydlpoly
    >> import molsys
    >> from molsys.util import ff2pydlpoly
    >> m = molsys.mol.fromFile("HKUST-1.mfpx")
    >> m.addon("ff")
    >> m.ff.read("HKUST-1")
    >> wm = ff2pydpoly.wrapper(,)
    >> pd = pydlpoly.pydlpoly("HKUST-1")
    >> pd.setup(web = wm)


LAMMPS Interface
================
TBI by Rochus


Stand-Alone functions and utilities
###################################
molsys features a lot of tools and data to be used in Molecular Modelling

Unit Cell Operations (unitcell.py)
==================================
the unitcell.py subpackage lets you convert cell parameters into cell vectors and vice versa

.. code-block:: python
    :linenos:

    import molsys.util.unit_cell as unit_cell
    
    ### assume you want to convert some given cellparams to the cell vectors
    cellparams = [12.0,12.5,10.0,90.0,90.0,120.0]
    cellvect   = unit_cell.vectors_from_abc(cellparams)
    print cellvect
    >> [[ 12.        0.        0.     ]
    >>  [ -6.25     10.82532   0.     ]
    >>  [  0.        0.       10.     ]]
    
    ### and vice versa
    cellparams = unit_cell.abc_from_vectors(cellvect)
    print cellparams
    >> [12.0, 12.5, 10.0, 90.0, 90.0, 119.99999999999999]

Matplotlib-based structure display (plotter.py)
===============================================


The Units package (units.py)
============================
The units subpackage of molsys contains natural constants and conversion factors for various generic units.
The definitions are all done in atomic units, conversion factors can be constructed in the following ways

.. code-block:: python
    :linenos:
    
    from molsys.util.units import *
    # --- input ---
    # some values with unit 
    E1 = 4 * kcalmol
    print 'E1_in_kjmol', E1/kjmol#
    >> E1_in_kjmol 16.736
    
    l1 = 1.54 *angstrom
    print 'l1 in meter', l1/meter
    >> l1 in meter 1.54e-10

    P= 766.0 * torr
    print 'P in bar', P/bar
    >> P in bar 1.02124652

Alignment of molecules' principal axes of inertia to the coordinate axes
========================================================================
Molecules can be rotated in such a way as to put their principal axes of inertia (which are the eigenvectors of their inertia tensor)
to perform this on a file, simply use the x2b shell script in molsys/scripts (should be in your path if molsys is properly installed) as

.. code-block:: ruby 
    
    align_pax input_file.ftype [optional: output_file.ftype]

if no output file is given, the input file is converted and overwritten.
File types are being taken from the filename

It can also be used in a python program. Use as

.. code-block:: python
    :linenos:

    import molsys
    import molsys.util.rotations as rotations
    m = molsys.mol(); m.read('some_molecule.ftype')
    ### not sure if u need that, but i suggest to use the bb addon to center the molecule first
    m.addon('bb')
    m.center_point = 'com'
    m.bb.center()
    m.xyz = rotations.align_pax(m.xyz,m.get_mass)
    ### the pax of m are now aligned to x,y and z. 


