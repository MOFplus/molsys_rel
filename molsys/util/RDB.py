"""Reaction Data Base

This is RDB which implements a pyDAL database wrapper for a reaction database

First version .. most likely the database strcutre is going to change
currently we use a sqlite database and store files in the storage folder

"""

from pydal import DAL, Field
import io
import os
from collections import OrderedDict

import molsys

# DB typekeys
typekeys = {
    "s" : "string",
    "t" : "text",
    "i" : "integer",
    "b" : "boolean",
    "u" : "upload",
    "r" : "reference",
    "d" : "double",
    "li": "list:integer",
    "ls": "list:string"
}

class RDB:

    def __init__(self, db_path, mode="a"):
        """generate a RDB access object
        
        Args:
            db_path (string): path to the sqlite database and the storage folder
            mode (str, optional): mode for opening database: a is append (must exist) n is new (must not exist). Defaults to "a".
        """
        db_path = os.path.abspath(db_path)
        print (db_path)
        # check access mode
        if mode == "a":
            assert os.path.isdir(db_path)
            assert os.path.exists(db_path+"/RDB.db")
        elif mode == "n":
            assert not os.path.exists(db_path)
            os.mkdir(db_path)
        self.db = DAL("sqlite://RDB.db", folder=db_path)
        self.db_path = db_path
        ###########################################################################################
        # define databse structure
        # 
        #           core of RDB -> here the tabel layout of the sql database is defined
        #
        dbstruc = OrderedDict()
        dbstruc["md"] = [       
            "s:path",         # filename of the pdlp file
            "s:stage",        # name of the stage
            "i:nframes",      # number of frames
            "d:timestep",     # time in fs between two frames (MD timestep times frame rate)
            "d:temp"          # temperature in Kelvin
        ]         # TBI !!! add more info here ... this is is just for testing
                  #         username, datetime, method, .... more
        dbstruc["revent"] = [
            "r:md",           # ref to table md
            "b:uni",          # true if unimolecular
            "i:frame",        # frame number
            "li:ed",          # educt species (sorted)
            "li:ts",          # transition state species
            "li:pr",          # product species
            "i:tr_ed",        # number of traced educt species
            "i:tr_pr",        # number of traced product species
            "li:rbonds",      # reactive bonds (list with 2*nbonds atom ids of the TS)
        ]
        dbstruc["md_species"] = [
            "r:revent",       # ref to revent
            "t:smiles",       # smiles
            "d:energy",       # ReaxFF energy
            "i:spec",         # species ID in the frame
            "i:foffset",      # frame offset (-1 educt, 0 TS, 1 product)
            "u:mfpx",         # upload mfpx file
            "b:tracked",      # is tracked?
        ]
        dbstruc["lot"] = [
            "s:name"          # level of theory
        ]
        dbstruc["opt_species"] = [
            "r:md_species",   # ref to md_species
            "r:lot",          # ref to lot
            "d:energy",       # energy (in kcal/mol)
            "u:xyz",          # upload xyz file
        ]
        dbstruc["react"] = [
            "r:revent:from_rev",        # source reaction event ref
            "r:revent:to_rev",          # target reaction event ref
            "r:md_species:from_spec",   # source species ref
            "r:md_species:to_spec",     # target species ref
        ]
        #############################################################################################
        self._define_tables(dbstruc)
        # define soem defaults
        self.current_md = None
        return

    def _define_tables(self, dbstruc, dryrun=False):
        """helper function to set up a database from a dictionary

        keys are teble names, vlaues are lists corresponding to fields with names split by :
        """ 
        # define soem defaults here
        db=self.db
        for t in dbstruc:
            fields = []
            for f in dbstruc[t]:
                kwargs = {}
                d = f.split(":")
                assert d[0] in typekeys
                # catch special field types
                if d[0] == "r":
                    # make sure that ref is in tables
                    assert d[1] in dbstruc
                    if len(d)>2:
                        # if a specific field name is requestd
                        fieldname = d[2]
                    else:
                        # else use tablenameID
                        fieldname = d[1] + "ID"
                    fieldtype = "reference %s" % d[1]
                elif d[0] == "u":
                    # generate an upload field .. by default the folder is in /storage/<fieldname>
                    fieldname = d[1]
                    fieldtype = typekeys[d[0]]
                    kwargs["uploadfolder"] = "%s/storage/%s" % (self.db_path, fieldname)
                else:
                    fieldname = d[1]
                    fieldtype = typekeys[d[0]]
                # add field
                field = Field(fieldname, type=fieldtype, **kwargs)
                if dryrun:
                    print (str(field))
                else:
                    fields.append(field) 
            # now all fields are defined ... generate the table
            if dryrun:
                print (" define table %s with the above fields" % t)
            else:
                db.define_table(t, fields)
        db.commit()
        return

    ######### API methods #####################################################################

    def commit(self):
        self.cd.commit()
        return

    def set_md_run(self, pdlp_fname, stage, **kwargs):
        # make sure that the pdlp file really exists
        assert os.path.exists(pdlp_fname)
        # get the absolute pathname to have a reliable identifier
        pdlp_fname = os.path.abspath(pdlp_fname)
        # find out if pdlp_fname is in database
        rows = self.db((self.db.md.path == pdlp_fname) & (self.db.md.stage == stage)).select()
        assert len(rows) < 2, "The pdlp file %s is twice in the database" % pdlp_fname
        if len(rows) == 0:
            # this is a new entry
            for k in ["nframes", "timestep", "temp"]:
                assert k in kwargs, "for a new md entry %s must be specified" % k
            mdID = self.db.md.insert(
                path     = pdlp_fname,
                stage    = stage,
                nframes  = kwargs["nframes"],
                timestep = kwargs["timestep"],
                temp     = kwargs["temp"]
            )
            self.db.commit()
        else:
            mdID = rows[0].id
        self.current_md = mdID
        return mdID

    def get_lot(self, lot):
        row = self.db(self.db.lot.name==lot).select().first()
        id = row.id
        if row is None:
            id = self.db.lot.insert(name=lot)
        return id

    def register_revent(self, frame, ed, ts, pr, tr_ed, tr_pr, rbonds, uni=False, postpone_commit=False):
        reventID = self.db.revent.insert(
            mdID       = self.current_md,
            uni        = uni,
            frame      = frame,
            ed         = ed,
            ts         = ts,
            pr         = pr,
            tr_ed      = tr_ed,
            tr_pr      = tr_pr,
            rbonds     = rbonds,
        )
        if not postpone_commit:
            self.db.commit()
        return reventID

    def get_revent(self, frame):
        event = self.db(self.db.revent.frame == frame).select().first()
        return event

    def add_md_species(self, reventID, mol, spec, foff, tracked=True, energy=0.0):
        # generate smiles
        mol.addon("obabel")
        smiles = mol.obabel.cansmiles
        # generate the file stream
        mfpxf = io.BytesIO(bytes(mol.to_string(), "utf-8"))
        # generate a filename
        fname = "%s_%d.mfpx" % (("ed", "ts", "pr")[foff+1], spec)
        # register in the database
        specID = self.db.md_species.insert(
            reventID     = reventID.id,
            smiles      = smiles,
            energy      = energy,
            spec        = spec,
            foffset     = foff,
            tracked     = tracked,
            mfpx        = self.db.md_species.mfpx.store(mfpxf, fname)
        )
        self.db.commit()
        return specID

    # TBI .. this is really stupid because we have to get revent for each species .. for DEBUG ok
    #        but merge these methods and make it more clever
    def get_revent_species(self,frame):
        # get requested revent of current md
        revent = self.db((self.db.revent.mdID == self.current_md) & (self.db.revent.frame == frame)).select().first()
        assert revent is not None, "No reaction event for frame %d" % frame
        return (revent.ed, revent.ts[0], revent.pr)


    def get_md_species(self, frame, spec, foff):
        # get requested revent of current md
        revent = self.db((self.db.revent.mdID == self.current_md) & (self.db.revent.frame == frame)).select().first()
        assert revent is not None, "No reaction event for frame %d" % frame
        # now get the species
        mdspec = self.db(
            (self.db.md_species.reventID == revent.id) &
            (self.db.md_species.foffset == foff) & 
            (self.db.md_species.spec == spec)
        ).select().first()
        assert mdspec is not None, "No species %d for this reaction event" % spec
        # DEBUG DEBUG
        print ("Frame %d Species %d frame offset %d" % (frame, spec, foff))
        print (mdspec.smiles)
        # get the file and convert to a mol object
        fname, mfpxf = self.db.md_species.mfpx.retrieve(mdspec.mfpx)
        mfpxs = mfpxf.read().decode('utf-8')
        mfpxf.close()
        mol = molsys.mol.from_string(mfpxs)
        return mol, fname

    def set_react(self, from_fid, to_fid, from_spec, to_spec):
        """connect reaction events -> make an edge in the reaction graph
        
        Args:
            from_fid (int): frame id of intial Product frame
            to_fid (int): frame id of final Educt frame
            from_spec (int): spec id of initial species (products)
            to_spec (int): spec id of final species  (educts)
        """
        from_ev = self.get_revent(from_fid-1)
        to_ev   = self.get_revent(to_fid+1)
        from_smd = self.db(
            (self.db.md_species.reventID == from_ev.id) &
            (self.db.md_species.foffset == 1) &
            (self.db.md_species.spec == from_spec)
        ).select().first()
        assert from_smd is not None, "no species %d in frame %d to connect" % (from_fid, from_spec)
        to_smd = self.db(
            (self.db.md_species.reventID == to_ev.id) &
            (self.db.md_species.foffset == -1) &
            (self.db.md_species.spec == to_spec)
        ).select().first()
        assert to_smd is not None, "no species %d in frame %d to connect" % (to_fid, to_spec)
        # now we can add a new edge into the reaction graph
        reactID = self.db.react.insert(
            from_rev  = from_ev.id,
            to_rev    = to_ev.id,
            from_spec = from_smd.id,
            to_spec   = to_smd.id 
        )
        self.db.commit()
        return

    def add_opt_species(self, mol, lot, energy, mdspecID):
        """add an optimized structure to the DB
        
        Args:
            mol (mol object): structure to be stored
            lot (string or int): name or id of the level of theory
            energy (float): energy of the system (unit is defiend by lot)
            mdspecID (int): reference id of the md_species entry
        """
        if type(lot) != type(1):
            lot = self.get_lot(lot)
        xyzf = io.BytesIO(bytes(mol.to_string(ftype="xyz"), "utf-8"))
        optID = self.db.opt_species.insert(
            md_speciesID = mdspecID,
            lotID        = lot,
            energy       = energy,
            xyz          = self.db.opt_species.xyz.store(xyzf, "opt.xyz")
        )
        self.db.commit()
        return
        


if __name__ == "__main__":
    rdb = RDB("./test", mode="n")


