#!/usr/bin/python

import MySQLdb
from   mysql  import search_db


########
def get_species (cursor):

    ensembl_db_name = {}
    all_species     = []

    qry = "show databases like '%core%'"
    cursor.execute(qry)

    rows = cursor.fetchall()
    if (not rows):
        print "No databases with 'core' in the name found"
        return 1

    for row in rows:
        db_name = row[0]
        name_tokens = db_name.split ('_')
        species = name_tokens[0]+'_'+ name_tokens[1]
        ensembl_db_name[species] = db_name
        all_species.append(species)

    return all_species, ensembl_db_name


########
def get_gene_ids (cursor, ensembl_db_name, biotype):

    gene_ids = []
    
    qry  = "use %s " % ensembl_db_name
    rows = search_db (cursor, qry)
    
    if (rows):
        rows = search_db (cursor, qry, verbose = True)
        print rows
        exit (1)

    qry = "select gene_id from gene where biotype='%s'" % biotype
    rows = search_db (cursor, qry)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
    
    else:
        for row in rows:
            gene_ids.append(row[0])
    
    return gene_ids
