#!/usr/bin/python

import MySQLdb
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

    for column_name in ['exon_key', 'ensembl_gene_id']:
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
def get_exon_dump_files(in_path):
    infiles = []
    for dirname, dirnames, filenames in os.walk(in_path):
        for filename in filenames:
            if ('exon_dump' in filename): infiles.append(filename)
    return infiles

#########################################
def store(cursor, in_path, infile, species):

    table =  'exon_' + species
    inf   = erropen (in_path+"/"+infile, "r")

    print "storing contents of ", infile

    ct = 0
    start = time()
    for line in inf:
        ct += 1
        if (not ct%1000):
            print "   %s   %5d    %8.3f" % (species, ct,  time()-start);
            start = time()
        line   = line.rstrip()
        field  = line.split("\t")

        exon_id         = int(field[0])
        ensembl_gene_id =     field[1]
        start_in_gene   = int(field[2])
        end_in_gene     = int(field[3])
        strand          = int(field[4])
        is_known        = int(field[5])
        is_coding       = int(field[6])
        is_canonical    = int(field[7])
        is_constitutive = int(field[8])
        species         =     field[9]
        source          =     field[10]
        protein_seq     =     field[11]
        left_flank      =     field[12]
        right_flank     =     field[13]
        dna_seq         =     field[14]

        fixed_fields    = {}
        update_fields   = {}

        exon_key = ensembl_gene_id +"_"+str(exon_id)+"_"+str(is_known)
        fixed_fields['exon_key']         = exon_key  
  
        update_fields['ensembl_gene_id'] =  ensembl_gene_id
        update_fields['start_in_gene']   =  start_in_gene 
        update_fields['end_in_gene']     =  end_in_gene  
        update_fields['strand']          =  strand      
        update_fields['is_known']        =  is_known    
        update_fields['is_coding']       =  is_coding    
        update_fields['is_canonical']    =  is_canonical 
        update_fields['is_constitutive'] =  is_constitutive
        update_fields['species']         =  species         
        update_fields['source']          =  source      
        update_fields['protein_seq']     =  protein_seq  
        update_fields['left_flank']      =  left_flank   
        update_fields['right_flank']     =  right_flank  
        update_fields['dna_seq']         =  dna_seq     

        store_or_update (cursor, table, fixed_fields, update_fields)


    inf.close()
    
#########################################
def load_from_infiles (infiles, in_path):
    
    db_name =  "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)
    
    ###############
    for infile in infiles:
        #if (not 'sus_scrofa' in infile): continue
        
        fields = infile.split ("_")
        species = fields[0] + "_" + fields[1] 
        if ('mustela') in fields[0]:
            species += "_" + fields[2]
        check_exon_table(cursor, db_name, species)
        
        table_size = check_exon_table_size (cursor, db_name, species)
        if table_size > 0: continue
        
        store  (cursor, in_path, infile, species)

    
#########################################
def main():

    
    no_threads = 12
    
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
    infiles = get_exon_dump_files(in_path)


    parallelize (no_threads, load_from_infiles, infiles, in_path)



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
