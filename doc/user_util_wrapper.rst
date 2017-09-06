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

