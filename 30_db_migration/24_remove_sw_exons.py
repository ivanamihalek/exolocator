#!/usr/bin/python

import MySQLdb, glob
import os, commands, sys
from   time import clock, time 
from   el_utils.mysql         import *
from   el_utils.config_reader import ConfigurationReader
from   el_utils.utils         import erropen
from   el_utils.threads       import  parallelize



#########################################
def main():

    
    no_threads = 1
    
    db_name =  "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)

    
    qry  = "show tables like 'exon_%'"
    rows = search_db(cursor, qry)
    print qry

    for row in rows:
        table = row[0]
        print table

        qry  = "delete from %s where source='sw_sharp'" % table
        rows = search_db(cursor, qry, verbose=True)
 
    
    cursor.close()
    db    .close()
    
 


#########################################
if __name__ == '__main__':
    main()

