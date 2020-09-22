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

    def __init__(self, db_path, mode="a", do_commit=True):
        """generate a RDB access object
        
        Args:
            db_path (string): path to the sqlite database and the storage folder
            mode (str, optional): mode for opening database: a is append (must exist) n is new (must not exist). Defaults to "a".
        """
        db_path = os.path.abspath(db_path)
        # check access mode
        if mode == "a":
            assert os.path.isdir(db_path)
            assert os.path.exists(db_path+"/RDB.db")
        elif mode == "n":
            assert not os.path.exists(db_path)
            os.mkdir(db_path)
        self.db = DAL("sqlite://RDB.db", folder=db_path)
        self.db_path = db_path
        self.do_commit = do_commit
        ###########################################################################################
        # define databse structure
        # 
        #           core of RDB -> here the tabel layout of the sql database is defined
        #

        # GS: Note maybe obvious, but as reminder for references the order in dbstruc matters
        dbstruc = OrderedDict()
        dbstruc["md"] = [       
            "s:path",         # filename of the pdlp file
            "s:stage",        # name of the stage
            "i:nframes",      # number of frames
            "d:timestep",     # time in fs between two frames (MD timestep times frame rate)
            "d:temp"          # temperature in Kelvin
        ]         # TBI !!! add more info here ... this is is just for testing
                  #         username, datetime, method, .... more
        dbstruc["unique_revent"] = [
            "b:uni",          # true if unimolecular
            "b:change",       # change in graph for educt and products (used for screening)
            "i:frame",        # frame number
            "li:ed",          # educt species (sorted)
            "li:ts",          # transition state species
            "li:pr",          # product species
            "i:tr_ed",        # number of traced educt species
            "i:tr_pr",        # number of traced product species
        ]
        dbstruc["revent"] = [
            "r:unique_revent",# ref to table unique_revent
            "b:reversed",     # is it the back reaction in unique_revent?
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
            "s:sumform",      # sum formula
            "d:energy",       # ReaxFF energy
            "i:spec",         # species ID in the frame
            "i:foffset",      # frame offset (-1 educt, 0 TS, 1 product)
            "u:mfpx",         # upload mfpx file
            "b:tracked",      # is tracked?
            "u:png",          # thumbnail
        ]
        dbstruc["lot"] = [
            "s:name"          # level of theory
        ]
        dbstruc["opt_species"] = [
            "r:md_species",   # ref to md_species
            "r:lot",          # ref to lot
            "d:energy",       # energy (in kcal/mol)
            "u:xyz",          # upload xyz file
            "u:png",          # thumbnail
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
        # define some defaults here
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
        self.db.commit()
        return

    def set_md_run(self, pdlp_fname, stage, **kwargs):
        # get the absolute pathname to have a reliable identifier
        pdlp_fname = os.path.abspath(pdlp_fname)
        # find out if pdlp_fname is in database
        rows = self.db((self.db.md.path == pdlp_fname) & (self.db.md.stage == stage)).select()
        assert len(rows) < 2, "The pdlp file %s is twice in the database" % pdlp_fname
        if len(rows) == 0:
            # this is a new entry
            for k in ["nframes", "timestep", "temp"]:
                assert k in kwargs, "for a new md entry %s must be specified" % k
            # make sure that the pdlp file really exists
            assert os.path.exists(pdlp_fname)
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
        if row is None:
            id = self.db.lot.insert(name=lot)
        else:
            id = row.id
        print (type(id))
        return id
    
    def register_unique_revent(self, frame, ed, ts, pr, tr_ed, tr_pr, uni=False, change=True):

        reventID = self.db.unique_revent.insert(
            uni        = uni,
            change     = change,
            frame      = frame,
            ed         = ed,
            ts         = ts,
            pr         = pr,
            tr_ed      = tr_ed,
            tr_pr      = tr_pr
        )

        if self.do_commit:
            self.db.commit()
        return reventID

    def register_revent(self, frame, ed, ts, pr, tr_ed, tr_pr, rbonds, uni=False):

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

        if self.do_commit:
            self.db.commit()
        return reventID

    def get_revent(self, frame):
        event = self.db(self.db.revent.frame == frame).select().first()
        return event

    def add_md_species(self, reventID, mol, spec, foff, tracked=True, energy=0.0):
        # generate smiles
        mol.addon("obabel")
        smiles = mol.obabel.cansmiles
        sumform = mol.get_sumformula()
        # generate the file stream
        mfpxf = io.BytesIO(bytes(mol.to_string(), "utf-8"))
        # generate a filename
        fname = "%s_%d.mfpx" % (("ed", "ts", "pr")[foff+1], spec)
        # register in the database
        specID = self.db.md_species.insert(
            reventID     = reventID.id,
            smiles      = smiles,
            sumform     = sumform,
            energy      = energy,
            spec        = spec,
            foffset     = foff,
            tracked     = tracked,
            mfpx        = self.db.md_species.mfpx.store(mfpxf, fname)
        )
        if self.do_commit:
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
        # print ("Frame %d Species %d frame offset %d" % (frame, spec, foff))
        # print (mdspec.smiles)
        return self.get_md_species_mol(mdspec)

    def get_md_species_mol(self, mdspec):
        """

        Args:
            - mdspec : a md_species database row
        """
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
        # assert from_smd is not None, "no species %d in frame %d to connect" % (from_spec, from_fid)
        to_smd = self.db(
            (self.db.md_species.reventID == to_ev.id) &
            (self.db.md_species.foffset == -1) &
            (self.db.md_species.spec == to_spec)
        ).select().first()
        # assert to_smd is not None, "no species %d in frame %d to connect" % (to_spec, to_fid)
        # now we can add a new edge into the reaction graph
        if (from_smd is not None) and (to_smd is not None):
            reactID = self.db.react.insert(
                from_rev  = from_ev.id,
                to_rev    = to_ev.id,
                from_spec = from_smd.id,
                to_spec   = to_smd.id 
            )
            if self.do_commit:
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
        if type(lot) == type(""):
            lot = self.get_lot(lot)
        xyzf = io.BytesIO(bytes(mol.to_string(ftype="xyz"), "utf-8"))
        optID = self.db.opt_species.insert(
            md_speciesID = mdspecID,
            lotID        = lot,
            energy       = energy,
            xyz          = self.db.opt_species.xyz.store(xyzf, "opt.xyz")
        )
        if self.do_commit:
            self.db.commit()
        return
        

################################################################################################

# reaction graph generation


    def view_reaction_graph(self, start=None, stop=None, browser="firefox", only_unique_reactions=False):
        """ generate a reaction graph

        we use the current md (must be called before)

        if you use svg the database must be locally available (only links to the figures are stored)
        png can get very big

        Args:
            name (string, optional): name of the output file, default = rgraph
            format (string, optional): format of the output (either png or svg), default = png
            start (int, optional) : first frmae to consider
            staop (int, optional) : last frame to consider
        """
        import pydot
        import tempfile
        import webbrowser
        import warnings
        # set the path to the images
        img_path = self.db_path + "/storage/png/"
        # get all relevant revents
        if start is None:
            start = -1

        if only_unique_reactions:
            revents = self.db((self.db.unique_revent.frame >= start)).select(orderby=self.db.unique_revent.frame)

            all_revents = self.db((self.db.revent.mdID == self.current_md) & \
                                  (self.db.revent.frame >= start)).select(orderby=self.db.revent.frame)
        else:
            revents = self.db((self.db.revent.mdID == self.current_md) & \
                              (self.db.revent.frame >= start)).select(orderby=self.db.revent.frame)


        rgraph = pydot.Dot(graph_type="digraph")
        rgnodes = {} # store all generated nodes by their md_speciesID
        # start up with products of first event
        cur_revent = revents[0]
        if only_unique_reactions:
            # Find a reaction event of this reaction class
            reventID = self.db(self.db.revent.unique_reventID == cur_revent).select()[0]

            mds = self.db((self.db.md_species.reventID == reventID) & \
                          (self.db.md_species.foffset == 1)).select()
        else:
            mds = self.db((self.db.md_species.reventID == cur_revent) & \
                          (self.db.md_species.foffset == 1)).select()
        for m in mds:
            new_node = pydot.Node("%d_pr_%d" % (cur_revent.frame, m.spec),
                                       image = img_path+m.png,
                                       label = "",
                                       shape = "box")
            rgraph.add_node(new_node)
            rgnodes[m.id] = new_node
        # now loop over revents
        for (i, cur_revent) in enumerate(revents[1:]):
            if (stop is not None) and (cur_revent.frame > stop):
                break
            # make the nodes of the revent
            educ = []
            prod = []
            if only_unique_reactions:

                if not cur_revent["change"]:
                    continue

                reventIDs = self.db(self.db.revent.unique_reventID == cur_revent).select()

                if len(reventIDs) > 0:
                    reventID = reventIDs[0]
                else:
                    #GS TODO this is a quick hack to make it run through. There is something fishy here
                    continue

                mds = self.db((self.db.md_species.reventID == reventID)).select()
            else:
                mds = self.db((self.db.md_species.reventID == cur_revent)).select()

            for m in mds:
                if m.foffset == -1:
                    new_node = pydot.Node("%d_ed_%d" % (cur_revent.frame, m.spec),\
                                       image = img_path+m.png,\
                                       label = "%s" % m.sumform,\
                                       labelloc = "t", \
                                       shape = "box")
                    educ.append(new_node)
                elif m.foffset == 1:
                    new_node = pydot.Node("%d_pr_%d" % (cur_revent.frame, m.spec),\
                                       image = img_path+m.png,\
                                       label = "%s" % m.sumform,\
                                       labelloc = "t", \
                                       shape = "box")
                    prod.append(new_node)
                else:
                    if cur_revent.uni:
                        label = "%10d (unimol)" % cur_revent.frame
                    else: 
                        label = "%10d" % cur_revent.frame
                    new_node = pydot.Node("%d_ts_%d" % (cur_revent.frame, m.spec),\
                                       image = img_path+m.png,\
                                       label = label,\
                                       labelloc = "t",\
                                       shape = "box",\
                                       style = "rounded")
                    ts = new_node
                rgraph.add_node(new_node)
                rgnodes[m.id] = new_node
            # now add edges
            for e in educ:
                rgraph.add_edge(pydot.Edge(e, ts))
            for p in prod:
                rgraph.add_edge(pydot.Edge(ts, p))
            # now connect from the previous events
            if only_unique_reactions:
                warnings.warn("Using only the unique reactions I can not generate a fully connected reaction graph!")
            else:
                concts = self.db(self.db.react.to_rev == cur_revent).select()
                for c in concts:
                    rgraph.add_edge(pydot.Edge(rgnodes[c.from_spec], rgnodes[c.to_spec], color="blue"))
        # done
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.curdir
            #os.chdir(tmpdir)
            rgraph.write_svg("rgraph.svg")
            webbrowser.get(browser).open_new("rgraph.svg")
            #os.chdir(cwd)
        
    



            

        





if __name__ == "__main__":
    rdb = RDB("./test", mode="n")

