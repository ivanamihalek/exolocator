#!/usr/bin/python

import MySQLdb, glob
import os, commands, sys
from   time import clock, time 
from   el_utils.mysql         import *
from   el_utils.config_reader import ConfigurationReader
from   el_utils.utils         import erropen
from   el_utils.threads       import  parallelize


########################################
def make_name_resolution_table (cursor):


    qry  = "CREATE TABLE name_resolution  (id INT(10) PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False
        

    for column_name in ['synonym',]:
        qry = "ALTER TABLE name_resolution  ADD %s blob" %  column_name
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['stable_id']:
        qry = "ALTER TABLE name_resolution  ADD %s VARCHAR(50)" %  column_name
        rows = search_db (cursor, qry)
        if (rows):
            return False

    if 0:
        qry = "create index key_id on  " + table + "(exon_key)";
        rows = search_db(cursor, qry)
        if rows:
            print rows
            exit(1)


#########################################
def store(cursor, in_path, infile):

    table = 'name_resolution'
    inf   = erropen (in_path+"/"+infile, "r")

    print "storing contents of ", infile
    for line in inf:


        fixed_fields    = {}
        update_fields   = {}


        line   = line.rstrip()
        fields = line.split("\t")
        if  'ENSG' in fields[-1]: 
            ensembl_gene_id = fields[-1]
        else:
            continue
        
        # check we are tracking that gene (for example, if it is pseudo, we are not)
        qry  = "select count(1) from exon_homo_sapiens where ensembl_gene_id  = '%s'" % ensembl_gene_id
        rows = search_db(cursor, qry)

        if not rows or not rows[0][0]: continue

        for field in fields[:-1]:

            if not field.replace (' ',''): continue

            fixed_fields ['synonym']   = field.replace("'", "").upper()
            update_fields['stable_id'] = ensembl_gene_id
            store_or_update (cursor, table, fixed_fields, update_fields)

    inf.close()
    
    
#########################################
def main():

    
    no_threads = 1
    
    db_name =  "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)

    cfg     = ConfigurationReader(user="marioot", passwd="tooiram", check=False)
    in_path = cfg.get_path('resources')
    if (not os.path.exists(in_path)):
        print in_path, "not found"

    
    ###############
    if not check_table_exists (cursor, db_name, 'name_resolution'):
        make_name_resolution_table (cursor)
   
    ###############
    os.chdir(in_path)
    filenames = glob.glob("*name_resolution.txt")
    
    for infile in filenames:
        store (cursor, in_path, infile)

    ###############
    cursor.close()
    db    .close()


#########################################
if __name__ == '__main__':
    main()

