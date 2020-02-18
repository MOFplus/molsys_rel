"""Reaction Data Base

This is RDB which implements a pyDAL database wrapper for a reaction database

First version .. most likely the database strcutre is going to change
currently we use a sqlite database and store files in the storage folder

"""

from pydal import DAL, Field
import io
import os

class RDB:

    def __intit__(self, db_path, mode="a"):
        """generate a RDB access object
        
        Args:
            db_path (string): path to the sqlite database and the storage folder
            mode (str, optional): mode for opening database: a is append (must exist) n is new (must not exist). Defaults to "a".
        """
        # check access mode
        if mode == "a":
            assert os.path.isdir(db_path)
            assert os.path.exists(db_path+"/RDB.db")
        elif mode == "n":
            assert not os.path.exists(db_path)
            os.mkdir(db_path)
        self.db = DAL("sqlite://RDB.db", folder=db_path)
        #
        self.define_tables()
        return

    def define_tables(self):
        db=self.db
        # MD runs
        db.define_table("MD",\
            Field("path", type="string"),\
            )
        # reaction events
        # MD structures
        # optimized structures
        # level of theory
        return

        

