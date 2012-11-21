#!/usr/bin/python
'''
Created on Apr 15, 2012

@author: Mario OT
'''
# Python imports
import os, sys, re
from mysql import connect_to_mysql, search_db, switch_to_db

       
class ConfigurationReader:
    '''
    Loads configuration files from the cfg database - which is assumend to be named '%exoloc%config%'
    '''

    def __init__ (self, user=None, passwd=None, host=None, port =None, check=True):
        
        self.util_path        = {}
        self.dir_path         = {}
        self.parameter_value  = {}
        self.cfg_db_name      = ""
        self.user             = user
        self.passwd           = passwd
        self.host             = host
        self.port             = port
        self.check            = check
        self.get_cfg_db()
        self.load_cfg()
        
        
    def get_cfg_db(self):
        db     = connect_to_mysql(self.user, self.passwd, self.host, self.port)
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
        #print "reading config from ", self.cfg_db_name
        cursor.close()
        db.close()


    def load_cfg(self):
        '''
        Load a configuration file and add the
        key, value pairs to the internal configuration dictionary
        '''
        db     = connect_to_mysql(self.user, self.passwd, self.host, self.port)
        cursor = db.cursor()
        switch_to_db (cursor, self.cfg_db_name)
        
        # utils
        qry = "select name, path from util_path"
        rows   = search_db (cursor, qry, verbose=False)
        if (not rows):
            print "no util_path info in %s" % self.cfg_db_name
            cursor.close()
            db.close()
            sys.exit (1)
        for row in rows:
            [name, path] = row
            if (self.check and not os.path.exists(path)):
                print path, " not found "
                sys.exit (1)
            self.util_path[name] = path

        # utils
        qry  = "select name, path from dir_path"
        rows = search_db (cursor, qry, verbose=False)
        if (not rows):
            print "no dir_path info in %s" % self.cfg_db_name
            cursor.close()
            db.close()
            sys.exit (1)
        for row in rows:
            [name, path] = row
            if (self.check and not os.path.exists(path)):
                print path, " not found "
                sys.exit (1)
            self.dir_path[name] = path


        # utils
        qry  = "select name, value from parameter"
        rows = search_db (cursor, qry, verbose=False)
        if (not rows):
            print "no dir_path info in %s" % self.cfg_db_name
            cursor.close()
            db.close()
            sys.exit (1)
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
            print name, value


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

    local_db = False
    if local_db:
        cfg      = ConfigurationReader()
    else:
        cfg    = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cfg.spill_all()
    print " ** ", cfg.get_path('mafft')
