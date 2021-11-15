# molsys

Molsys is the library used to read, write and manipulate atomistic models.
It can be used in python programs and features also lots of scripts to run on the command line to transform structure files.
For a rationale why such a class is needed please have a look at the documentation (see below)

### Installing

In order to install molsys, clone this repository into the destination of your choosing (we always use /home/%USER/sandbox/, also the installation instructions below use this path)

```
https://github.com/MOFplus/molsys.git
```
or if you have put an ssh key to github
```
git@github.com:MOFplus/molsys.git
```

Afterwards the PATH and PYTHONOATH have to be updated. Add to your .bashrc :
```
export PYTHONPATH=/home/$USER/sandbox/molsys:$PYTHONPATH
export PATH=/home/$USER/sandbox/molsys/scripts:$PATH
```

There are several dependencies needed for certain addons. If you want to use them you need to install these dependencies, which can be painful.
We strongly recommend to install molsys via the [cmc-tools](https://github.com/MOFplus/cmc-tools) repository and follow the installation docu in the [README](https://github.com/MOFplus/cmc-tools#readme) there. It expalins to use a conda environment, which makes this process rather easy. 


## Building the Documentation

Mandatory dependencies to built the documentation can be obtained via pip:
```
pip install Sphinx
pip install sphinx-rtd-theme
```
The Documentation can be compiled by running
```
make html
```
in the doc folder.
A Built directory containing
```
/built/html/index.html
```
was created. It can be opened with the browser of your choice

## License

MIT License

