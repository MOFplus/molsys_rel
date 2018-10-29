INTRO
This README file explains how to test ACAB. Basic features, such as edge color ratio per vertex and angles btw. given colored edges, are tested.

USAGE
run tests:
    pytest
clean output of tests:
    ./clean_tests.sh
(clean after each test otherwise run directory index increments...)

OUTPUT
Each directory contains:
- the "grey" symmetry structure, defining the space group symmetry in which colors are searched
- the colors to be used to weave frameworks
- a pretty view of the colors, cutting out the edges with periodic boundary conditions
File mfpx and txyz are given for each structure.
N.B.: the pretty view may be misleading since not every edge is visible.


TBI
Self-explanatory name of output directories. Runs are ordered, but nets are not specified in the directory... and each test increases the directory index.
