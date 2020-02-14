#!/bin/bash
#TODO: refactoring, continous integration

SCRIPT_FOLDER=$(realpath $(dirname $0))
cd $SCRIPT_FOLDER

    mol/run.sh
    mol/clean.sh

    lqg/run.sh
    lqg/clean.sh

    ptg/run.sh
    ptg/clean.sh

    toper/run.sh
    toper/clean.sh

cd - >/dev/null
