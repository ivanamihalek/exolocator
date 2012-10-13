#!/usr/bin/python

import MySQLdb
import os
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.mysql   import  check_table_exists, store_or_update
'''
  Set paths to Ensembl directories and to various utility programs.
'''
#########################################
def  make_path_table (cursor, table):

    print "making ", table

    qry  = "create table " + table + "  (id int(10) primary key auto_increment)"
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        return False

    # make the columns
    column  = 'name'
    qry = "alter table  %s  add  %s  varchar (20)" % (table, column)
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        return False

    column  = 'path'
    qry = "alter table  %s  add  %s  blob" % (table, column)
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        return False


#########################################
def  make_seqregion2file_table (cursor):

    table = 'seqregion2file'

    
    qry  = "create table " + table + "  (seqregion_id int(10) primary key)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    # make the columns
    column  = 'seq_name'
    qry = "alter table  %s  add  %s  varchar (100)" % (table, column)
    rows = search_db (cursor, qry)
    if (rows):
        return False

    column  = 'file_name'
    qry = "alter table  %s  add  %s  blob" % (table, column)
    rows = search_db (cursor, qry)
    if (rows):
        return False


#########################################
def make_table (cursor, table):
    
    if   table == 'util_path':
        make_path_table (cursor, table)
    elif table == 'dir_path':
        make_path_table (cursor, table)
    elif table == 'seqregion2file':
        make_seqregion2file_table (cursor)
       
    else:
        print "I don't know how to make table '%s'" % table




#########################################
def main():

    
    db_name   = "exolocator_config"
    util_path = {}
    util_path['mafft']    = '/usr/local/bin/mafft'
    util_path['blastall'] = '/usr/bin/blastall'
    util_path['fastacmd'] = '/usr/bin/fastacmd'
    
    dir_path = {}
    dir_path['ensembl_fasta'] = '/home/ivanam/databases/ensembl/fasta'

    # check if the paths are functioning (at this point at least)
    for util in util_path.values():
        if (not os.path.exists(util)):
            print util, " not found "
            exit (1)

    for dir in dir_path.values():
        if (not os.path.exists(dir)):
            print dir, " not found "
            exit (1)
        if (not os.path.isdir (dir)):
            print dir, " is not a directory "
            exit (1)
            
    db     = connect_to_mysql()
    cursor = db.cursor()

    # check if the config db exists -- if not, make it
    qry  = "show databases like'%s'" % db_name
    rows = search_db (cursor, qry)
    if (not rows):
        print db_name, "database not found"
        qry = "create database %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            print "some problem creating the database ..."
            rows = search_db (cursor, qry, verbose = True)
    else:
        print db_name, "database  found"

    qry = "use %s " % db_name
    search_db (cursor, qry)
        
    # make tables
    for table in ['util_path', 'dir_path', 'seqregion2file']:
        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            make_table (cursor,  table)
   
    # fill util and and dir path tables
    for [name, path] in util_path.iteritems():
        fixed_fields  = {}
        update_fields = {}
        fixed_fields['name'] = name
        fixed_fields['path'] = path
        store_or_update (cursor, 'util_path', fixed_fields, update_fields)

    for [name, path] in dir_path.iteritems():
        fixed_fields  = {}
        update_fields = {}
        fixed_fields['name'] = name
        fixed_fields['path'] = path
        store_or_update (cursor, 'dir_path', fixed_fields, update_fields)

    cursor.close()
    db.close()



#########################################
if __name__ == '__main__':
    main()
