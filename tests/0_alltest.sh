#!/bin/bash
#just run this executable script:
#   ${DIRECTORY_OF_THIS_SCRIPT}/0_alltest.sh
#or simplier, if your working directory is already the directory of this script:
#   ./0_alltest.sh
#TODO: refactoring, continous integration

SCRIPT_FOLDER=$(realpath $(dirname $0))
cd $SCRIPT_FOLDER

    mol/run.sh
    mol/clean.sh

    ptg/run.sh
    ptg/clean.sh


cd - >/dev/null
