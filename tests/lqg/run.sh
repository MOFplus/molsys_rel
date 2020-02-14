SCRIPT_FOLDER=$(realpath $(dirname $0))
cd $SCRIPT_FOLDER
    rm -r __pycache__ .pytest_cache molsys.log .hypothesis 2>/dev/null
    python3 --disable-warnings lqg_test.py
cd - >/dev/null
