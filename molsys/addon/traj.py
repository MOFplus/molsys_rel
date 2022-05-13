"""
             traj addon

RS 2020: an addon to attach a trajectory source to a mol object

you can advance a frame or pick a selected frame, define a source (should be essentially a numpy xyz array with proper size
and potentially a cell array)

supports also the writing to a trajectory file (just a concatenation ... use pdb if possible .. mfpx does not support trajectories)

"""

from molsys.util import mfp5io
from molsys.fileIO import formats
import numpy as np

# this is a fake force field engine ... just to have something that owns a mol object
class fake_ffe:

    def __init__(self, mol):
        self.mol = mol
        return


class traj:

    def __init__(self, mol, source="mfp5", **kwargs):
        self._mol = mol
        self.ffe = fake_ffe(mol) 
        # set defaults
        self.open_mfp5 = None
        self.fid = 0
        self.fmax = 0
        self.variable_cell = False
        # start up
        self.set_source(source, **kwargs)
        return

    def set_source(self, source, **kwargs):
        assert source in ["mfp5", "array", "xyz"]
        if source == "mfp5":
            # define mfp5 file by fname and stage and set arrays from there
            if self.open_mfp5 is not None:
                # there was a mfp5 as a source before  (maybe another stage?) close it!
                self.open_mfp5.close()
                self.variabale_cell = False
                self.fmax = 0
            # we compare with the mol object attached to the fake ffe to make sure that natoms etc matches
            self.open_mfp5 = mfp5io.mfp5io(kwargs["fname"], ffe=self.ffe, restart=kwargs["stage"])
            # get data from traj group
            traj = self.open_mfp5.get_traj_from_stage(kwargs["stage"])
            assert traj is not False, "Stage %s does not extist" % kwargs["stage"]
            assert len(traj) > 0 , "No trajectory info in that stage"
            assert "xyz" in traj, "Trajectory does not contain xyz coordinate info"
            self.xyz = traj["xyz"]
            self.fmax = self.xyz.shape[0]
            if "cell" in traj:
                self.cell = np.array(traj["cell"])
                self.variable_cell = True
        elif source == "array":
            self.xyz = kwargs["array"]
            self.fmax = self.xyz.shape[0]
            if "cell_array" in kwargs:
                self.cell = kwargs["cell_array"]
                self.variable_cell = True
        elif source == "xyz":
            # in thsi case we reopen the xyz file and read only the xyz positions. no variable cell
            # also elements are skipped
            # TBI make parallel
            axyz = []
            fxyz = open(kwargs["fname"], "r")
            stop = False
            while not stop:
                line = fxyz.readline()
                if len(line) > 0:
                    natoms = int(line.split()[0])
                    assert natoms == self._mol.get_natoms()
                    fxyz.readline()
                    xyz = []
                    for i in range(natoms):
                        line = fxyz.readline()
                        xyz.append([float(v) for v in line.split()[1:4]])
                    axyz.append(xyz)
                else:
                    stop = True
            fxyz.close()
            self.xyz = np.array(axyz)
            self.fmax = self.xyz.shape[0]
        # check if the number of atoms in the provided xyz array is correct
        assert self.xyz.shape[1] == self._mol.get_natoms()
        assert self.xyz.shape[2] == 3
        # TBI : check shapes of cell arrays
        self.set_frame(self.fid)
        return

    def set_frame(self, fid):
        assert fid >= 0
        assert fid < self.fmax
        self.fid = fid
        self._mol.set_xyz(np.array(self.xyz[fid]).astype("float64")) # in case of a mfp5 file -> make sure we get the right dtype
        if self.variable_cell:
            self._mol.set_cell(self.cell[fid], cell_only=True)
        return

    def next_frame(self, step=1):
        self.fid += step
        self.set_frame(self.fid)
        return

    def write(self, fname, first=0, last=None, stride=1):
        ftype = fname.rsplit(".", 1)[-1]
        assert ftype in ["pdb", "xyz"] , "Filetype %s not allowed for trajectory writing" % ftype
        if self._mol.mpi_rank == 0:
            f = open(fname, "w")
            if last is None:
                last = self.fmax
            for i in range(first, last, stride):
                self.set_frame(i)
                formats.write[ftype](self._mol, f)
            f.close()
        return


