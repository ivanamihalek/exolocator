#!/usr/bin/python
'''
Created on Apr 15, 2012

@author: Mario OT
'''
# Python imports
import os, sys, re
from mysql import connect_to_mysql, connect_to_db, search_db

       
class ConfigurationReader:
    '''
    Loads configuration files from the cfg database - which is assumend to be named '%exoloc%config%'
    '''

    def __init__ (self):
        
        self.util_path = {}
        self.dir_path  = {}
        self.parameter_value  = {}
        self.cfg_db_name = ""
        self.get_cfg_db()
        self.load_cfg()
        
    def get_cfg_db(self):
        db     = connect_to_mysql()
        cursor = db.cursor()
        qry    = "show databases like'%exoloc%config%'" 
        rows   = search_db (cursor, qry, verbose=False)
        if (not rows):
            print "config database not found (by guessing)"
            print "pls fix config_reader.py"
            cursor.close()
            db.close()
            exit (1)
        self.cfg_db_name = rows[0][0]
        print "reading config from ", self.cfg_db_name
        cursor.close()
        db.close()


    def load_cfg(self):
        '''
        Load a configuration file and add the
        key, value pairs to the internal configuration dictionary
        '''
        db     = connect_to_db(self.cfg_db_name)
        cursor = db.cursor()
        
        # utils
        qry = "select name, path from util_path"
        rows   = search_db (cursor, qry, verbose=False)
        if (not rows):
            print "no util_path info in %s" % self.cfg_db_name
            cursor.close()
            db.close()
            exit (1)
        for row in rows:
            [name, path] = row
            self.util_path[name] = path

        # utils
        qry  = "select name, path from dir_path"
        rows = search_db (cursor, qry, verbose=False)
        if (not rows):
            print "no dir_path info in %s" % self.cfg_db_name
            cursor.close()
            db.close()
            exit (1)
        for row in rows:
            [name, path] = row
            self.dir_path[name] = path


        # utils
        qry  = "select name, value from parameter"
        rows = search_db (cursor, qry, verbose=False)
        if (not rows):
            print "no dir_path info in %s" % self.cfg_db_name
            cursor.close()
            db.close()
            exit (1)
        for row in rows:
            [name, value]        = row
            self.parameter_value[name] = value

        cursor.close()
        db.close()
   
    def spill_all(self):
        for name, path in self.dir_path.iteritems():
            print name, path
        for name, path in self.util_path.iteritems():
            print name, path
        for name, value in self.parameter_value.iteritems():
            print name, path


    def get_path (self, name):
        if name in self.dir_path.keys():
            return self.dir_path[name]
        elif name in self.util_path.keys():
            return self.util_path[name]
        else:
            return ""
        
    def get_value (self, name):
        if name in self.parameter_value.keys():
            return self.parameter_value[name]
        else:
            return ""
        
        
if __name__ == '__main__':
    cr = ConfigurationReader()
    cr.spill_all()
    print " ** ", cr.get_path('mafft')
