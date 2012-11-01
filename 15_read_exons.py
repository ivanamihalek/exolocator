#!/usr/bin/python

import MySQLdb
import os, commands, sys
from   time import clock, time 
from   el_utils.mysql         import connect_to_mysql, search_db, switch_to_db
from   el_utils.mysql         import store_or_update, check_table_exists
from   el_utils.config_reader import ConfigurationReader
from   el_utils.utils         import erropen


########################################
def make_exon_table (cursor):

    table = 'exon'

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
def check_exon_table(cursor, db_name):
    table = 'exon'
    if ( check_table_exists (cursor, db_name, table)):
        print table, " found in ", db_name
    else:
        print table, " not found in ", db_name
        make_exon_table (cursor)
        create_index (cursor, db_name,'exon_key_idx', table, ['exon_key'])

#########################################
def get_exon_dump_files(in_path):
    infiles = []
    for dirname, dirnames, filenames in os.walk(in_path):
        for filename in filenames:
            if ('exon_dump' in filename): infiles.append(filename)
    return infiles

#########################################
def store(cursor, in_path, infile):

    fields  = infile.split("_")
    species = "_".join(fields[0:2])
    inf     = erropen (in_path+"/"+infile, "r")

    print
    print "storing contents of ", infile
    print species

    ct = 0
    start = time()
    for line in inf:
        ct += 1
        if (not ct%1000):
            print "   %s   %5d    %8.3f" % (species, ct,  time()-start);
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

        store_or_update (cursor, 'exon', fixed_fields, update_fields)


    inf.close()
    
#########################################
def main():

    db_name =  "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)

    cfg      = ConfigurationReader(user="marioot", passwd="tooiram", check=False)
    in_path = cfg.get_path('afs_dumps')
    if (not os.path.exists(in_path)):
        print in_path, "not found"

    ###############
    check_exon_table(cursor, db_name)
    
    ###############
    infiles = get_exon_dump_files(in_path)

    ###############
    for infile in infiles:
        store(cursor, in_path, infile)

    cursor.close()
    db    .close()



#########################################
if __name__ == '__main__':
    main()
