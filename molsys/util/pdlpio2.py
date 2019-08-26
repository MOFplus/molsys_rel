"""

                            pdlpio2

    implements a hdf5/h5py based class to handle storage of pydlpoly
    molecular systems (used for restart and trajectory data)

    the class object is a frontend to a hdf5/h5py file, and it can be used
    outside of pydlpoly to read and postprocess trajectory data

    This is a rewrite of the V1.0 that ships with molsys and is based on a 
    molsys type mol object. It is meant to work with any MM engine that provides
    a molsys mol object and provides a certain API. It can be used with pydlpoly
    (using the new startup mode with molsys) or pylmps currently. Note that in pylmps 
    trajectory info is written in C++ in lammps using dump_pdlp.cpp (which needs to be
    compiled in)

    one part of the data is "fixed"
    - number of atoms
    - elements
    - atom types
    - connectivity
    - boundary conditions
    - forcefield parameters
    in a restart this data is used as a replacement to the usual file input

    the other part is "staged"
    in other words if needed multiple stages can be defined like "equilibration",
    "sample1", "sample2" etc.
    the initial and default stage is called "default" by default ... yes!
    for each stage current restart information (position/velocitites/cell) is
    automatically stored every N steps. restart data is stored in double precision
    in addition, trajectory info is stored every M steps. (both N and M can be changed)
    trajectory data (pos or pos/vel or pos/vel/force) can be stored both in single and
    double precision on request.
    
    Versions (in file attribute "version"
    - no version info == 1.0 : original
    - version 1.1            : nstep of trajectories is stored per dataset and not for the stage
    - version 1.2            : pdlpio2 version to work with new mol objects (can write par info)
    
"""

import h5py
import numpy as np

import types
import os

from molsys import mpiobject
import molsys


class pdlpio2(mpiobject):

    def __init__(self, fname, ffe=None, restart=None, mpi_comm=None, out=None):
        """generates a pdlp file object

        This object knows three states: conencted to a ffe (force field engine) (with a mol object),
        connected to a ffe but without mol object, which implies a restart (will generate mol from pdlp file)
        or the analysis mode. In the analysis mode writing to trajectory data is forbidden.
        In the connected mode one can request appending to a stage. 
        Otherwise a new stage is required.

        NOTE: always use self.open() as a first command in a method that will access the file. it will check if the file
              is already open and only open if it is not open. in certain cases (e.g. dumping with lammps) the file needs to be 
              closed before passing to lammps. To make sure the file is connected call open().
        
        Args:
            fname (string): filename (including .pdlp extension)
            ffe (object, optional): Defaults to None. Force Field engine (must contain a ffe.mol instance of type molsys)
            restart (string, optional): Defaults to None. Name of the stage from which a restart should be read
        """ 
        super(pdlpio2, self).__init__(mpi_comm,out)
        self.verbose = 0
        # helper object for hdf5 variable length strings
        self.str_dt = h5py.special_dtype(vlen=str)
        #
        self.fname = fname
        # check if this file exists or should be generated
        self.fexists = os.path.isfile(self.fname)
        self.h5file_isopen = False
        # check if we run with an FFE or not
        self.ffe  = ffe
        if ffe is not None:
            self.ffe_name = type(ffe).__name__
            if restart is None:
                # a force field engine is connected => check if mol object is available and of the right type
                try:
                    self.mol = ffe.mol
                except AttributeError:
                    self.pprint("PDLP ERROR: Force Field Engine has no mol object!")
                    raise
                # ok, now let us check if mol is actually a molsys object (require it to have a is_topo attribute)
                try:
                    is_topo = self.mol.is_topo
                except AttributeError:
                    self.pprint("PDLP ERROR: mol object seems not be a molsys object! This is required for pdlpio2")
                    raise
                # last but not least check that mol is not a topo ... we do not want that
                if is_topo is True:
                    self.pprint("PDLP ERROR: mol object is a topo object .. that should not happen!")
                    raise Exception("pdlp error")
                self.mode = "ffe"
                # now all is fine .. set defaults for this runmode
                self.stage = "default" # default stage which can be overwritten (contains no trajectory data, only restarts)
            else:
                # this is a restart, which means the file must exist and also the stage name
                assert self.fexists, "PDLP ERROR: For restart the file must exist. Check filename and path!"
                self.restart_stage = restart
                self.mode = "restart"
        else:
            # if ffe is None we are in analysis mode (means the file must exist)
            assert self.fexists, "PDLP ERROR: For analysis the file must exist!"
            self.mode = "analysis"
        # now open the file
        self.open()
        #
        self.file_version = 1.2
        if self.fexists is True:
            if self.mode is "ffe" or self.mode is "restart":
                # this is an exisiting file and we are in ffe/restart mode and want to add data => make sure all is matching
                if self.is_master:
                    file_version = self.h5file.attrs["version"]
                else:
                    file_version = None
                file_version = self.mpi_comm.bcast(file_version)  
                assert self.file_version == file_version, "PDLP ERROR: Exisiting file has a different version! Can not add data"
                if self.mode is "ffe":
                    # check the system if it contains the same molecule (we could be more clever here, but if the list of atomtypes is identical we should be safe)
                    # TBI mor consistent tests like paramters etc .. we could read the complete mol out and compare
                    self.compare_system()
        else:
            # file is new .. set version
            if self.is_master:
                self.h5file.attrs["version"] = self.file_version
            # this must be ffe mode so we can start generating the system group
            self.write_system()
            # new file so we write the initital structure to the restart
            self.add_stage(self.stage)
            self.write_restart()
        return
        
    def __del__(self):
        self.close()
        return

    def open(self):
        """ method to make sure that h5file is open """
        if self.h5file_isopen:
            return
        self.h5file_isopen = True
        if self.is_master:
            self.h5file = h5py.File(self.fname)
        else:
            self.h5file = None
        return

    def close(self):
        """ method to make sure that the file is closed """
        if not self.h5file_isopen:
            return
        self.h5file_isopen = False
        if self.is_master:
            self.h5file.close()
        return         

    # system IO
    #
    # The system group contains all the time independent information of the system (connectivity, atomtypes, parameter)
    # Note: in the old pdlpio system this data was prepared by the mol (assign_FF) class and just written here
    #       Now in pdlpio2 we expect a certain interface, namely the mol object from molsys and collect all the 
    #       data in an active way

    def write_system(self):
        """ Writes all the info from a mol object to the system group in the pdlp file
        
        The following info is written:
            - elems       (strings)
            - atypes      (strings)
            - fragtypes   (strings)
            - fragnumbers (int)
            - bcd         (int)
        
        Further info is stored in subgroups for different addons
            - ff addon (parameter and rics)
            - molecules addon (TBI)

        """
        self.open()
        if self.is_master:
            system = self.h5file.require_group("system")
            # elems
            na = self.ffe.mol.get_natoms()
            elems = self.ffe.mol.get_elems()
            pdlp_elems = system.require_dataset("elems", shape=(na,), dtype=self.str_dt)
            pdlp_elems[...] = elems
            # atypes
            atypes = self.ffe.mol.get_atypes()
            pdlp_atypes = system.require_dataset("atypes", shape=(na,), dtype=self.str_dt)
            pdlp_atypes[...] = atypes
            # fragtypes
            fragtypes = self.ffe.mol.get_fragtypes()
            pdlp_fragtypes = system.require_dataset("fragtypes", shape=(na,), dtype=self.str_dt)
            pdlp_fragtypes[...] = fragtypes
            # fragnumbers
            fragnumbers = self.ffe.mol.get_fragnumbers()
            pdlp_fragnumbers = system.require_dataset("fragnumbers", shape=(na,), dtype="i")
            pdlp_fragnumbers[...] = fragnumbers
            # get connectivity table
            cnc_table = self.ffe.mol.get_ctab()
            if len(cnc_table) > 0:
                # catch if there are no bonds at all: h5py does not like zero size selections
                cnc_table = np.array(cnc_table, dtype="i")
                pdlp_cnc_table = system.require_dataset("cnc_table", shape=cnc_table.shape, dtype=cnc_table.dtype)
                pdlp_cnc_table[...] = cnc_table
            system.attrs["bcd"] = self.ffe.mol.get_bcond()
            if "ff" in self.ffe.mol.loaded_addons:
                # a force field is loaded: we assume it is also initialized. TBI: a switch in the ff addon to verify this
                data = self.ffe.mol.ff.pack()
                ff = system.require_group("ff")
                ric = ff.require_group("ric")
                par = ff.require_group("par")
                if data["FF"] is not None:
                    par.attrs["FF"] = data["FF"]
                else:
                    par.attrs["FF"] = ""
                for r in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
                    # write ric arrays to pdlp file
                    if r in data:
                        d = data[r]
                        fd = ric.require_dataset(r, shape=d.shape, dtype="i")
                        fd[...] = d
                    # write par data to pdlp file
                    #  order is ptype, names, npars, pars
                    if r+"_par" in data:
                        p = data[r+"_par"]
                        fp = par.require_group(r)
                        n = len(p[0]) # number of paramters for this ric type
                        fptypes = fp.require_dataset("ptypes", shape=(n,), dtype=self.str_dt)
                        fptypes[...] = p[0]
                        fnames = fp.require_dataset("names", shape=(n,), dtype=self.str_dt)
                        fnames[...] = p[1]
                        fnpars = fp.require_dataset("npars", shape=p[2].shape, dtype="i")
                        fnpars[...] = p[2]
                        fpars = fp.require_dataset("pars", shape=p[3].shape, dtype="float64")
                        fpars[...] = p[3]
        return

    def compare_system(self):
        """ method to verify that the mol object in the ffe is consistent with the stuff on file 
            NOTE: this can get arbitrary complex .. do a simple test of atypes here """
        self.open()
        OK = None
        if self.is_master:
            system = self.h5file["system"]
            pdlp_atypes = system["atypes"]
            atypes = self.ffe.mol.get_atypes()
            OK = True
            if len(pdlp_atypes)==len(atypes):
                for a1,a2 in zip(pdlp_atypes, atypes):
                    if a1 != a2:
                        OK = False
            else:
                OK = False
        OK = self.mpi_comm.bcast(OK)
        assert OK, "PDLP ERROR: The system in the pdlp file is not equivalent to your actual system. Aborting!"
        return

    def get_mol_from_system(self):
        """ read mol info from system group and generate a mol object 

        in parallel this is done on the master only and the data is broadcasted to the other nodes
        """
        self.open()
        # read the data from file on the master
        if self.is_master:
            system = self.h5file["system"]
            elems  = list(system["elems"])
            atypes = list(system["atypes"])
            fragtypes = list(system["fragtypes"])
            fragnumbers = list(system["fragnumbers"])
            if "cnc_table" in list(system.keys()):
                cnc_table = list(system["cnc_table"])
            else: 
                cnc_table = []
            bcd = system.attrs["bcd"]
        # broadcast if we run parallel
        if self.mpi_size > 1:
            if self.is_master:
                data = (elems, atypes, fragtypes, fragnumbers, cnc_table, bcd)
            else:
                data = None
            elems, atypes, fragtypes, fragnumbers, cnc_table, bcd = self.mpi_comm.bcast(data)
        # generate mol object
        na = len(elems)
        mol = molsys.mol()
        mol.set_natoms(na)
        mol.set_elems(elems)
        mol.set_atypes(atypes)
        mol.set_fragnumbers(fragnumbers)
        mol.set_fragtypes(fragtypes)
        mol.set_ctab(cnc_table, conn_flag=True)
        mol.bcond = bcd
        # now read the restart info that needs to be passed to the mol instance (Note: no velocities etc are read here)
        if self.is_master:
            try:
                rstage = self.h5file[self.restart_stage]
            except KeyError:
                self.pprint("PDLP ERROR: The requested restart stage %s does not exist in the file!" % self.restart_stage)
                raise
            restart = rstage["restart"]
            xyz = np.array(restart["xyz"])
            cell = np.array(restart["cell"])
        if self.mpi_size>1:
            if self.is_master:
                self.mpi_comm.Bcast(xyz)
                self.mpi_comm.Bcast(cell)
            else:
                xyz = np.empty([na, 3], dtype="float64")
                cell = np.empty([3, 3], dtype="float64")
                self.mpi_comm.Bcast(xyz)
                self.mpi_comm.Bcast(cell)
        mol.set_xyz(xyz)
        mol.set_cell(cell)
        # new check if addon data is present in the system group
        # start with ff
        ff_data = None
        if self.is_master:
            if "ff" in list(system.keys()):
                # ok, there is force field data in this file lets read it in as a packed directory
                ff_data = {}
                ff = system["ff"]
                ric = ff["ric"]
                par = ff["par"]
                ff_data["FF"] = par.attrs["FF"]
                for r in ["bnd", "ang", "dih", "oop", "cha", "vdw"]:
                    # read in ric integer arrays
                    if r in ric:
                        ff_data[r] = np.array(ric[r])
                    # now get the parameter
                    if r in par:
                        rpar = par[r]
                        ff_data[r+"_par"] = (\
                                        list(rpar["ptypes"]),\
                                        list(rpar["names"]),\
                                        np.array(rpar["npars"]),\
                                        np.array(rpar["pars"]))
        if self.mpi_size>1:
            ff_data = self.mpi_comm.bcast(ff_data)
        if ff_data is not None:
            # yes there was ff data .. set up mol addon ff
            mol.addon("ff")
            mol.ff.unpack(ff_data)
        return mol

    def add_stage(self, stage):
        """add a new stage 
        
        Args:
            stage (string): stage name to add

        Returns:
            (bool) True if it worked and False if the stage already existed
        """
        if stage == "system":
            self.pprint("PDLP ERROR: Stage name system is not allowed!")
            raise IOError
        OK = True
        if self.is_master:
            if stage in list(self.h5file.keys()):
                OK = False
            else:
                sgroup = self.h5file.create_group(stage)
                rgroup = sgroup.create_group("restart")
                # generate restart arrays for xyz, cell and vel
                na = self.ffe.get_natoms()
                rgroup.require_dataset("xyz",shape=(na,3), dtype="float64")
                rgroup.require_dataset("vel",shape=(na,3), dtype="float64")
                rgroup.require_dataset("cell",shape=(3,3), dtype="float64")
        OK = self.mpi_comm.bcast(OK)
        return OK

    def write_restart(self, stage=None, velocities=False):
        """write restart info to file

        Args:
            stage (string, optional): Defaults to None. Name of the stage to write restart (must exist). Uses current stage by default
            velocities (bool, optional): Defaults to False. When True writes also velocities
        """
        self.open()
        if stage is None:
            stage = self.stage
        OK = True
        if self.is_master:
            if stage not in list(self.h5file.keys()):
                OK = False
            else:
                restart = self.h5file[stage+"/restart"]
                xyz = self.ffe.get_xyz()
                cell = self.ffe.get_cell()
                rest_xyz = restart.require_dataset("xyz",shape=xyz.shape, dtype=xyz.dtype)
                rest_xyz[...] = xyz
                rest_cell = restart.require_dataset("cell", shape=cell.shape, dtype=cell.dtype)
                rest_cell[...] = cell
                if velocities:
                    vel = self.ffe.get_vel()
                    rest_vel = restart.require_dataset("vel", shape=vel.shape, dtype=vel.dtype)
        OK = self.mpi_comm.bcast(OK)
        if not OK:
            self.pprint("PDLP ERROR: writing restart to stage %s failed. Stage does not exist" % stage)
            raise IOError
        return

    def prepare_stage(self, stage, traj_data, traj_nstep, data_nstep=1, prec="float64", tstep=0.001):
        """prepare a stage for trajectory writing
        
        Args:
            stage (string): name of the stage (must exist)
            traj_data (list of strings): name of the data to be written
            traj_nstep (int): frequency in steps to write
            data_nstep (int or list of ints), optional): Defaults to 1. freq to write each datatype
            prec (str, optional): Defaults to "float64". precision to use in writing traj data
            tstep (float, optional): Defaults to 0.001. Timestep in ps
        """ 
        self.open()       
        if type(data_nstep)== type([]):
            assert len(data_nstep)==len(traj_data)
        else:
            data_nstep = len(traj_data)*[data_nstep]
        OK=True
        if self.is_master:
            if stage not in list(self.h5file.keys()):
                OK = False
            else:
                traj = self.h5file.require_group(stage+"/traj")
                traj.attrs["nstep"] = traj_nstep
                traj.attrs["tstep"] = tstep
                for dname,dnstep in zip(traj_data, data_nstep):
                    assert dname in list(self.ffe.data_funcs.keys())
                    data = self.ffe.data_funcs[dname]()
                    dshape = list(data.shape)
                    pdlp_data = traj.require_dataset(dname, 
                                    shape=tuple([1]+dshape),
                                    maxshape=tuple([None]+dshape),
                                    chunks=tuple([1]+dshape),
                                    dtype=prec)
                    pdlp_data[...] = data
                    pdlp_data.attrs["nstep"] = dnstep
        OK = self.mpi_comm.bcast(OK)
        if not OK:
            self.pprint("PDLP ERROR: preparing stge %s failed. Stage does not exist!" % stage)
            raise IOError
        return

        
    def __call__(self, force_wrest=False):
        """ this generic routine is called to save restart info to the hdf5 file
            if force_wrest is True then restart is written in any case """
        if not self.pd: raise IOError("No pydlpoly instance")
        data_written = False
        if (self.counter%self.rest_nstep == 0) or force_wrest:
            # restart is written .. fetch all data 
            for d in self.rest_data:
                if self.verbose:
                    self.pprint("Writing restart data %s to pdlp file" % d)
                data = self.data_funcs[d]()
                self.rest_datasets[d][...] = data
            data_written = True
        if self.track_data != None:
            for i,d in enumerate(self.track_data):
                if (self.counter%self.traj_nstep[i] == 0):
                    self.traj_frame[i] += 1
                    data = self.data_funcs[d]()
                    tds = self.traj_datasets[d]
                    tds.resize(self.traj_frame[i], axis=0)
                    tds[self.traj_frame[i]-1,...] = data
                    data_written = True
        # now we are done 
        if data_written : self.h5file.flush()
        self.counter += 1
        return
        
            
        