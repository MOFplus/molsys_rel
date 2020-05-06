"""
    thumbnail

    takes a mol object and generates a png file

    prerequisites: needs vmd and imagemagick "convert" to be installed

    please help to improve options
"""


import tempfile
import subprocess
import os

def thumbnail(mol, size=400, scale=1.3, transparent=True, fname=None, debug=False):
    """
    generate a thumbnail from a mol object
    by default a png is returned.

    we generate the tcl commands on the fly in order to change stuff there
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        if fname is None:
            fname = os.getcwd() + "/" + mol.name + ".png"
        else:
            fname = os.getcwd() + "/" + fname
        # files
        xyzf = tmpdir+"/mol.xyz"
        tclf = tmpdir+"/vmd.tcl"
        tgaf = tmpdir+"/vmd.tga"
        # store structure as xyz
        mol.write(xyzf)
        # write vmd commands
        f = open(tclf, "w")
        vmd_settings = """
        color Display Background white
        color Name C black
        axes location Off
        mol modstyle 0 0 CPK 0.70000 0.300000 12.000000 12.000000
        scale by %f
        """ % scale
        # add additional stuff here ... optional
        f.write(vmd_settings)
        f.write("render TachyonInternal %s\n" % tgaf)
        f.write("exit\n")
        f.close()
        # now run vmd
        vmdoutf = open(tmpdir+"/vmd.out", "w")
        subprocess.run(["vmd", "-xyz", xyzf, "-size", str(size), str(size), "-e", tclf, "-dispdev", "text"], stdout=vmdoutf)
        # now convert to png
        convert = ["convert"]
        if transparent:
            convert += ["-transparent", "white"]
        convert += [tgaf, fname]
        subprocess.run(convert)
        if debug==True:
            print ("In DEBUG mode")
            print ("go with another shell to tempdir %s to see intermediate files")
            input ("press ENTER to end an delete everything")
    return
