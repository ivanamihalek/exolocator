#!/usr/bin/python

import MySQLdb
from   mysql  import search_db

 

########
def stable2gene (cursor, db_name=None, stable_id=None):

    if (not stable_id):
        return ""

    if  (db_name):
        qry  = "use %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            rows = search_db (cursor, qry, verbose = True)
            print rows
            exit (1)


    qry = "select gene_id from gene where stable_id='%s'" % stable_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return int(rows[0][0])
    
########
def gene2stable (cursor, db_name=None, gene_id=None):

    if (not gene_id):
        return ""

    if  (db_name):
        qry  = "use %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            rows = search_db (cursor, qry, verbose = True)
            print rows
            exit (1)


    qry = "select stable_id from gene where gene_id=%d" % gene_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return rows[0][0]
    

########
def get_gene_ids (cursor, db_name=None, biotype = None, is_known = None):

    gene_ids = []
    
    if  (db_name):
        qry  = "use %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            rows = search_db (cursor, qry, verbose = True)
            print rows
            exit (1)

    qry = "select gene_id from gene"

    if ( biotype or is_known):
        qry +=  " where "
        if ( biotype):
           qry += "biotype='%s'" % biotype
        if (biotype and is_known):
            qry += " and "
        if (is_known):
           qry += "status='known'"

    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return []
    else:
        if ('Error' in rows[0]):
            print rows[0]
            return []

        for row in rows:
            if ( not type(row[0]) is long ):
                print row
                exit(1)
            gene_ids.append(int(row[0]))
    
    return gene_ids


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

