#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.mysql   import  create_index, check_table_exists, check_column_exists
from   el_utils.ensembl import  get_species, get_gene_ids

#########################################
def make_sw_exon_table (cursor):


    # if maps_to_human_exon_id is 0
    # the region was searched, but nothing was found
    # has NNN refers to the fact that the searched region  contains NNN stretch,
    #     indicating that the exon might not have been sequenced
    # 5p_ss and 3p_ss refer to canonical splice sites -r' and t'  refer to the intron

    table = 'sw_exon'

    qry  = "CREATE TABLE " + table + "  (exon_id INT(10) PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    for column_name in ['gene_id', 'start_in_gene', 
                        'end_in_gene', 'maps_to_human_exon_id', 'exon_seq_id']:
        qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['strand', 'phase',  'has_NNN', 'has_stop', 'has_3p_ss', 'has_5p_ss']:
        qry = "ALTER TABLE %s  ADD %s tinyint" %  (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

#########################################
def make_gene2exon_table (cursor):


    table = 'gene2exon'



    qry  = "CREATE TABLE " + table + " (gene2exon_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    for column_name in ['gene_id', 'exon_id', 'start_in_gene', 'end_in_gene', 
                        'canon_transl_start', 'canon_transl_end', 'exon_seq_id']:
        qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    
    for column_name in ['strand', 'phase',  'is_known', 
                        'is_coding', 'is_canonical', 'is_constitutive']:
        qry = "ALTER TABLE %s  ADD %s tinyint" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['covering_exon']:
        qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['covering_is_known']:
        qry = "ALTER TABLE %s  ADD %s tinyint" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['analysis_id']:
        qry = "ALTER TABLE %s  ADD %s INT" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

#########################################
def make_exon_seq_table (cursor):


    table = 'exon_seq'

    qry  = "CREATE TABLE " + table + "  (exon_seq_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    for column_name in ['exon_id']:
        qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    
    for column_name in ['is_known', 'is_sw']:
        qry = "ALTER TABLE %s  ADD %s tinyint" %  (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['dna_seq', 'left_flank', 'right_flank', 'protein_seq']:
        qry = "ALTER TABLE %s  ADD %s blob" %  (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    qry  = "alter TABLE " + table + "  add column pepseq_transl_start int(10)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    qry  = "alter TABLE " + table + "  add column pepseq_transl_end int(10)"
    rows = search_db (cursor, qry)
    if (rows):
        return False


#########################################
def make_coding_region_table(cursor):


    table = 'coding_region'

    qry  = "CREATE TABLE " + table + " (gene_id INT PRIMARY KEY)"
    rows = search_db (cursor, qry)
    if (rows):
        return False


    # make the columns
    columns = ['start', 'end']
    for column in columns:
        qry = "ALTER TABLE %s ADD %s  INT(10)" % (table, column)
        rows = search_db (cursor, qry)
        if (rows):
            return False
     

#########################################
def make_orthologue_table (cursor, table):
    

    # if congamet_gene_id is 0, and source is 'rbh'
    # means that the reciprocal-best-hit was attempted but nothing was found

    qry  = "CREATE TABLE " + table + "  (orth_pair_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    for column_name in ['gene_id', 'cognate_gene_id']:
        qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['cognate_genome_db_id']:
        qry = "ALTER TABLE %s  ADD %s INT" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    
    for column_name in ['source']:
        qry = "ALTER TABLE %s  ADD %s VARCHAR(20)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

#########################################
def make_exon_map_table (cursor):
    
    table = 'exon_map'

    qry  = "CREATE TABLE " + table + "  (exon_map_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    for column_name in ['exon_id', 'cognate_exon_id']:
        qry = "ALTER TABLE %s add %s INT(10)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['exon_known', 'cognate_exon_known']:
        qry = "ALTER TABLE %s add %s tinyint" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['cognate_genome_db_id']:
        qry = "ALTER TABLE %s add %s INT" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['cigar_line']:
        qry = "ALTER TABLE %s add %s blob" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False
    
    for column_name in ['similarity']:
        qry = "ALTER TABLE %s add %s float" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False
    
    for column_name in ['source']:
        qry = "ALTER TABLE %s add %s VARCHAR(20)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['msa_bitstring']:
        qry = "ALTER TABLE %s add %s varbinary(1000)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False


#########################################
def make_table (cursor, db_name, table):
    

    qry = "use %s" % db_name
    rows = search_db (cursor, qry, verbose=False)
    if (rows):
        return False


    if   table == 'gene2exon':
        make_gene2exon_table (cursor)
    elif table == 'exon_seq':
        make_exon_seq_table (cursor)
    elif table == 'sw_exon':
        make_sw_exon_table (cursor)
    elif table == 'coding_region':
        make_coding_region_table (cursor)
    elif table in ['orthologue', 'unresolved_ortho', 'paralogue']:
        make_orthologue_table (cursor, table)
    elif table == 'exon_map':
        make_exon_map_table (cursor)

    else:
        print "I don't know how to make table '%s'" % table


    
#########################################
def add_filename_column (cursor, db_name):
    
    qry = "alter table seq_region add file_name blob"
    rows = search_db (cursor, qry)
    if (rows):
        return False
  
    return True
   
    
#########################################
def modify_filename_column (cursor, db_name):
    
    qry = "alter table seq_region modify column  file_name blob"
    rows = search_db (cursor, qry)
    if (rows):
        return False
  
    return True
   

#########################################
def main():
    
    #db     = connect_to_mysql()
    db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)


    # add exon tables to all species
    for species in all_species:
        print species
        db_name = ensembl_db_name[species]
        for table in ['gene2exon', 'exon_seq', 'sw_exon', 'coding_region']:

            if ( check_table_exists (cursor, db_name, table)):
                print table, " found in ", db_name
            else:
                print table, " not found in ", db_name
                make_table (cursor, db_name, table)
 
        create_index (cursor, db_name, 'eg_index', 'gene2exon', ['exon_id', 'gene_id'])
        create_index (cursor, db_name, 'gene_id_idx', 'gene2exon', ['gene_id'])
        create_index (cursor, db_name, 'ek_index', 'exon_seq', ['exon_id', 'is_known'])

    # add file_name column to seq_region table
    for species in all_species:
        print species
        db_name = ensembl_db_name[species]

        if ( check_column_exists (cursor, db_name, "seq_region", "file_name")):
            print "file_name found in seq_region, ", db_name
            #modify_filename_column (cursor, db_name)
        else:
            print "file_name  not found in seq_region, ", db_name
            add_filename_column (cursor, db_name)
        
    
    # add orthologue table to human - we are human-centered here
    # ditto for map (which exons from other species map onto human exons)
    print "adding orthologue to human"
    species = 'homo_sapiens'
    db_name = ensembl_db_name[species]
    for table in ['orthologue', 'unresolved_ortho', 'paralogue', 'exon_map']:
        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            make_table (cursor, db_name, table)
        if table == 'exon_map':
            create_index (cursor, db_name,'gene_index', table, ['exon_id'])
        else:
            create_index (cursor, db_name,'gene_index', table, ['gene_id'])

    

    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
