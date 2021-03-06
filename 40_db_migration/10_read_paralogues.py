#!/usr/bin/python -u

import MySQLdb, glob
import os, commands, sys
from   time import clock, time 
from   el_utils.mysql         import connect_to_mysql, search_db, switch_to_db
from   el_utils.mysql         import store_or_update, check_table_exists, create_index
from   el_utils.config_reader import ConfigurationReader
from   el_utils.utils         import erropen


#########################################
def store(cursor, infile):

    inf = erropen(infile, "r")

    total        = 0
    id_not_found = 0
    for line in inf:
        line.rstrip()
        total += 1
        if not total%1000: print "\t", total
        if ( len(line.split()) !=  2 or not 'ENS' in line):
            continue
        [stable_id1, stable_id2] = line.split()
        fixed_fields    = {}
        update_fields   = {}
        
        fixed_fields['gene_id1'] = stable_id1
        fixed_fields['gene_id2'] = stable_id2

        store_or_update (cursor, 'paralog', fixed_fields, update_fields)

    print "done with ", infile, "total ",  total

    inf.close ()


#########################################
def main():

    db_name = "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)

    cfg      = ConfigurationReader (user="marioot", passwd="tooiram", check=False)
    in_path  = cfg.get_path('afs_dumps')
    in_path += "/para_dump"
    if (not os.path.exists(in_path)):
        print in_path, "not found"
        sys.exit(1) # exit on non-existent outdir

    
    ###############
    if 1:
        qry = "drop table paralog"
        search_db (cursor, qry)
        qry = "create table paralog (id int(10) primary key auto_increment) "
        search_db (cursor, qry)
        qry = "alter table paralog  ADD gene_id1 varchar(30) " 
        search_db (cursor, qry)
        qry = "alter table paralog  ADD gene_id2 varchar(30) " 
        search_db (cursor, qry)
        create_index (cursor, db_name,'gene_id_index', 'paralog', ['gene_id1', 'gene_id2'])
        

    ###############
    os.chdir(in_path)
    filenames = glob.glob("*_para_dump.txt")

    ###############
    for infile in filenames:
        print infile
        store(cursor, infile)

    cursor.close()
    db    .close()



#########################################
if __name__ == '__main__':
    main()
