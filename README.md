# molsys

Molsys is the library used to read, write and manipulate atomistic models.
It can be used in python programs and features also lots of scripts to run on the command line to transform structure files.

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

Mandatory dependencies are:

* numpy (pip install numpy)
* 

Some addons and utility modules require additional packages to be installed. These are:

* [graph-tool](https://git.skewed.de/count0/graph-tool/wikis/installation-instructions#installation-via-package-managers) (molsys.addon.graph) (ppa via apt-get)
* spglib (molsys.addon.spg) (pip)
* pandas (molsys.addon.zmat) (pip)
* scipy (molsys.addon.zmat) (pip)
* [ff_gen](https://github.com/MOFplus/ff_gen) (molsys.addon.ric) (github)

## Running the tests

There will soon be a testing framework framework available.


## Contributing

Any changes to the main mol class (mol.py) have to be assured by [Rochus Schmid](https://github.com/rochusschmid)

## License

TBA

## Acknowledgments

TBA
