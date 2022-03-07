import os
import time
import argparse
from config import Config

"""Databse Management System.
This utility is for the MOLBAR project.
The goal is to retain metadata about every compute and every output file.
MOLBAR programs should use this system to log their results.
This utility is not thread safe.
The output files do not have unique names (potential for overwrite on rerun).
"""

class ScienceDB:
    """Databse Management System class.
    """
    VERSION = 2.0
    FILENAME = "science.db"

    def format_date(self,year,month,day):
        show = "%04d/%02d/%02d" % (year,month,day)
        return show

    def format_file_date (self,filename):
        year=0
        month=0
        day=0
        try:
            last_modified_time = os.stat(filename).st_mtime
            time_tuple = time.localtime(last_modified_time)
            year = time_tuple[0]
            month = time_tuple[1]
            day = time_tuple[2]
            return self.format_date(year,month,day)
        except Exception:
            print("WARNING: logger cannot determine date of this file: "+filename)
        return self.format_date(year,month,day)

    def __init__(self,filename=FILENAME):
        self.filename = filename
        exists = os.path.isfile(self.filename)
        self.theDB = Config(self.filename)

    def set_version(self):
        """This is a no-op."""
        return True

    def add_data_file_auto_date(self,filename,filetype,description):
        """Record file metadata.
        """
        filedate = self.format_file_date(filename)
        self.add_data_file(filename,filetype,filedate,description)
        return True

    def add_data_file (self,filename,filetype,filedate,description):
        """Record file metadata.
        """
        resource=filename
        #
        attribute="filetype"
        value=filetype
        self.theDB.add_value(resource,attribute,value)
        #
        attribute="description"
        value=description
        self.theDB.add_value(resource,attribute,value)
        datestamp = self.format_file_date(filename)
        #
        attribute="datestamp"
        value=datestamp
        self.theDB.add_value(resource,attribute,value)
        return True

    def add_process (self,procname,command,infiles,outfiles):
        """Record process metadata."""
        resource=procname
        #
        attribute="command"
        value=command
        self.theDB.add_value(resource,attribute,value)
        #
        attribute="infiles"
        value=infiles
        self.theDB.add_value(resource,attribute,value)
        #
        attribute="outfiles"
        value=outfiles
        self.theDB.add_value(resource,attribute,value)
        return True

