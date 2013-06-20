#!/usr/bin/python

import MySQLdb, glob
import os, commands, sys
from   time import clock, time 
from   el_utils.mysql         import *
from   el_utils.config_reader import ConfigurationReader
from   el_utils.utils         import erropen
from   el_utils.threads       import  parallelize


########################################
def make_exon_table (cursor, table):

    qry  = "CREATE TABLE " + table + "  (id INT(10) PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    for column_name in ['exon_key', 'ensembl_gene_id', 'ensembl_exon_id']:
        qry = "ALTER TABLE %s  ADD %s VARCHAR(50)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False
        

    for column_name in ['start_in_gene',  'end_in_gene']:
        qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['strand', 'is_known', 'is_coding', 'is_canonical', 'is_constitutive']:
        qry = "ALTER TABLE %s  ADD %s tinyint" %  (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['species', 'source', 'protein_seq', 'left_flank', 'right_flank', 'dna_seq' ]:
        qry = "ALTER TABLE %s  ADD %s blob" %  (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False


#########################################
def check_exon_table(cursor, db_name, species, verbose = False):
    table =  'exon_' + species
    
    if ( check_table_exists (cursor, db_name, table)):
        if verbose: print table, " found in ", db_name
        #qry = "drop table "+table
        #rows = search_db(cursor, qry)
        #make_exon_table (cursor, table)
        #if rows:
        #    return rows[0][0]
        #else:
        #   return 0
    else:
        if verbose: print table, " not found in ", db_name
        make_exon_table (cursor, table)



#########################################
def check_exon_table_size(cursor, db_name, species):
    table =  'exon_' + species

    qry  = "select count(1) from " + table
    rows = search_db(cursor, qry)

    if rows:
        return rows[0][0]
    else:
        return 0


#########################################
def store(cursor, in_path, infile):

    inf   = erropen (in_path+"/"+infile, "r")

    print "storing contents of ", in_path, " file ", infile
    
    ct = 0
    start = time()
    for line in inf:
        ct += 1
        if (not ct%1000):
            print "     %5d    %8.3f" % (ct,  time()-start);
            start = time()

        fixed_fields    = {}
        update_fields   = {}


        line   = line.rstrip()
        field = line.split("\t")
        if len(field) < 4: continue
        [human_stable_id, cognate_stable_id, species, common_name]  =  field
  
        fixed_fields ['ensembl_gene_id'] = human_stable_id  
        fixed_fields ['species']         = species  
  
        update_fields['cognate_gene_id'] = cognate_stable_id
        update_fields['common_name']     = common_name
 
        store_or_update (cursor, 'ortholog', fixed_fields, update_fields)

    inf.close()
    
#########################################
def load_from_infiles (infiles, in_path):
    
    db_name =  "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)
    
    ###############
    infiles.reverse()
    for infile in infiles:
         
        store  (cursor, in_path, infile)
        #print "\t done in  %8.3f sec" % (time()-start) 
        
    
#########################################
def main():

    
    no_threads = 1
    
    db_name =  "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)

    cfg      = ConfigurationReader(user="marioot", passwd="tooiram", check=False)
    in_path = cfg.get_path('afs_dumps')
    if (not os.path.exists(in_path)):
        print in_path, "not found"


    
    cursor.close()
    db    .close()
    
    ###############
    os.chdir(in_path)
    filenames = glob.glob("orthologue_dump.txt")
    
    parallelize (no_threads, load_from_infiles, filenames, in_path)



#########################################
if __name__ == '__main__':
    main()

'''
    if (0):
    qry = "show columns from " + table + " like 'dna_seq'"
    rows      =  search_db(cursor, qry)
    seq_space =  rows[0][1]
    print seq_space
    if (not seq_space == 'blob'):
        print "drop"
        qry = "drop table "+table
        print qry
        rows      =  search_db(cursor, qry)
        print rows
        make_exon_table (cursor, table)
    print
'''
