OWD=$PWD
SCRIPT_FOLDER=$(realpath $(dirname $0))
cd $SCRIPT_FOLDER
    for i in */
    do
        cd $i
        ./clean_examples.sh
        cd ..
    done
cd $OWD