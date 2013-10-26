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
def check_exon_table(cursor, db_name, table, verbose = False):
    
    date = table_create_time (cursor, db_name, table)
    print date
    exit(1)

    if 0:
        if ( check_table_exists (cursor, db_name, table)):
            if verbose: print table, " found in ", db_name
            #qry = "drop table "+table
            rows = search_db(cursor, qry)
            make_exon_table (cursor, table)

            if rows:
                return rows[0][0]
            else:
               return 0
            qry = "create index key_id on  " + table + "(exon_key)";
            rows = search_db(cursor, qry)
            if rows:
                print rows
                return 0
        else:
            if verbose: print table, " not found in ", db_name
            make_exon_table (cursor, table)

    return

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
def store(cursor, table,  in_path, infile, species):

    inf   = erropen (in_path+"/"+infile, "r")
    if not inf: exit(1) 


    print "storing contents of ", in_path, " file ", infile
    
    ct = 0
    start = time()
    for line in inf:
        ct += 1
        if (not ct%1000):
            print "   %s   %5d    %8.3f" % (species, ct,  time()-start);
            start = time()

        fixed_fields    = {}
        update_fields   = {}


        line   = line.rstrip()
        field  = line.split("\t")
        if len(field) < 18: continue

        exon_id         = int(field[0])
        ensembl_gene_id =     field[1]
        ensembl_exon_id =     field[2]
        start_in_gene   = int(field[3])
        end_in_gene     = int(field[4])
        strand          = int(field[5])
        is_known        = int(field[6])
        is_coding       = int(field[7])
        is_canonical    = int(field[8])
        is_constitutive = int(field[9])
        species         =     field[10]
        source          =     field[11]
        #if source == 'sw_sharp' or source=='usearch':
        #    human_exon      = field[12]
        #    protein_seq     = field[13]
        # here I have two fields showing where the peptide translation starts and where it ends
        #   left_flank      = field[16]
        #    right_flank     = field[17]
        #    dna_seq         = field[18]
        #    fixed_fields['maps_to_human_exon_id'] = human_exon
        #else:
        protein_seq     =     field[12]
        # here I have two fields showing where the peptide translation starts and where it ends
        left_flank      =     field[15]
        right_flank     =     field[16]
        dna_seq         =     field[17]


        exon_key = ensembl_gene_id + "_" + str(exon_id) + "_" + str(is_known)
        fixed_fields ['exon_key']        = exon_key  
  
        update_fields['ensembl_gene_id'] = ensembl_gene_id
        update_fields['ensembl_exon_id'] = ensembl_exon_id
        update_fields['start_in_gene']   = start_in_gene 
        update_fields['end_in_gene']     = end_in_gene  
        update_fields['strand']          = strand      
        update_fields['is_known']        = is_known    
        update_fields['is_coding']       = is_coding    
        update_fields['is_canonical']    = is_canonical 
        update_fields['is_constitutive'] = is_constitutive
        update_fields['species']         = species         
        update_fields['source']          = source      
        update_fields['protein_seq']     = protein_seq  
        update_fields['left_flank']      = left_flank   
        update_fields['right_flank']     = right_flank  
        update_fields['dna_seq']         = dna_seq     

        store_or_update (cursor, table, fixed_fields, update_fields)

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
        #if 'pong_abelii' in infile: continue
        start = time()
        print "reading ", infile
        fields  = infile.split ("_")
        species = fields[0] + "_" + fields[1] 
        if ('mustela') in fields[0]:
            species += "_" + fields[2]
            

        table =  'exon_' + species
        check_exon_table (cursor, db_name, table, verbose = True)

        #store  (cursor, table, in_path, infile, species)
        #print "\t %s  done in  %8.3f sec" % (species, time()-start) 
       
    
#########################################
def main():

    
    no_threads = 1
    
    db_name =  "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)

    cfg      = ConfigurationReader(user="marioot", passwd="tooiram", check=False)
    # afs is killing me here ...
    in_path  = cfg.get_path('afs_dumps')+"/exons"
    if (not os.path.exists(in_path)):
        print in_path, "not found"


    
    cursor.close()
    db    .close()
    
    ###############
    os.chdir(in_path)
    filenames = glob.glob("*exon_dump.txt")
    
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
