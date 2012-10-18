#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db, 
from   el_utils.mysql   import  create_index, check_table_exists
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

    for column_name in ['gene_id', 'exon_id', 'start_in_gene', 'end_in_gene', 'exon_seq_id']:
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
    else:
        print "I don't know how to make table '%s'" % table


    
    

#########################################
def main():
    
    db     = connect_to_mysql()
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
        
        create_index (cursor, db_name,'eg_index', 'gene2exon', ['exon_id', 'gene_id'])
        create_index (cursor, db_name,'ek_index', 'exon_seq', ['exon_id', 'is_known'])


    # add orthologue table to human - we are human-centered here
    print "adding orthologue to human"
    species = 'homo_sapiens'
    db_name = ensembl_db_name[species]
    for table in ['orthologue', 'unresolved_ortho', 'paralogue']:
        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            make_table (cursor, db_name, table)

        create_index (cursor, db_name,'gene_index', table, ['gene_id'])

    

    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
