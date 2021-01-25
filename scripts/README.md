# Instructions on how the scripts work could go here ...

## moldenx / vmdx
moldenx/vmdx can be used as a .mfpx file viewer. It preprocesses the file to a tinker xyz file, calls molden and does all the tmpfile cleanup ( thx to roberto ;) )
#### MOLDEN installation
On recent debian-based linux distros, our good old molden binary does not work any more ... 
to compile [molden 5.8](https://www3.cmbi.umcn.nl/molden/howtoget.html) you may need the following packages:

```
sudo apt-get install libxmu-dev libxmu-headers
sudo apt-get install xutils-dev
```
#### VMD installation
Go to [the webpage](http://www.ks.uiuc.edu/Research/vmd/), create an account if you don't have one already, download the recent version and follow the installation instructions.

