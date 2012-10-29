#!/usr/bin/python

import MySQLdb
from   mysql   import search_db, switch_to_db
from   exon    import  Exon

#########################################
def gene2exon_list (cursor, gene_id, db_name=None):

    exons = []

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    qry  = "select * from gene2exon where gene_id = %d " % gene_id
    rows = search_db(cursor, qry)
    if (not rows):
        rows = search_db(cursor, qrym, verbose = True)
        exit (1)

    for row in rows:
        exon = Exon()
        exon.load_from_gene2exon(row)
        exons.append(exon)

    return exons
########################################
def  get_canonical_exons (cursor, gene_id):

    exons = gene2exon_list (cursor, gene_id)
    if (not exons):
        print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
        exit(1) # shouldn't happen at this point

    # sorting exons in place by their start in gene:
    exons.sort(key=lambda exon: exon.start_in_gene)

    canonical_coding_exons = []
    reading = False
    for exon in exons:
        if (not exon.is_canonical): 
             continue
        if (not exon.canon_transl_start is None):
            reading = True
        if (reading):
            canonical_coding_exons.append(exon)
        if (not exon.canon_transl_end is None):
            break  

    return canonical_coding_exons


#########################################
def get_selenocysteines (cursor, gene_id):

    selenoC_pos = []

    canonical_transl_id = gene2canon_transl(cursor, gene_id)

    qry  = "select value from translation_attrib "
    qry += " where attrib_type_id = 12 and  translation_id = %d " % canonical_transl_id

    rows = search_db (cursor, qry)
    if (not rows):
       return []

    for row in rows:
        blah  = row[0].split (" ")
        start = int(blah[0])
        end   = int(blah[1])
        for pos in range(start, end+1):
            selenoC_pos.append(pos-1)

    return selenoC_pos


#########################################
def gene2stable_canon_transl(cursor, gene_id, db_name=None):

    if (db_name and not switch_to_db(cursor, db_name)):
            return False
 
    qry  = "select translation.stable_id  from translation, gene "
    qry += " where gene.canonical_transcript_id = translation.transcript_id "
    qry += " and gene.gene_id = %d " % gene_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return  rows[0][0]

#########################################
def gene2canon_transl(cursor, gene_id, db_name=None,):

    if  (db_name and not switch_to_db(cursor, db_name)):
            return False
 
    qry  = "select translation.translation_id  from translation, gene "
    qry += " where gene.canonical_transcript_id = translation.transcript_id "
    qry += " and gene.gene_id = %d " % gene_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return  rows[0][0]




########
def stable2gene (cursor, stable_id=None, db_name=None, ):

    if (not stable_id):
        return ""

    if (db_name and not switch_to_db(cursor, db_name)):
            return False

    qry = "select gene_id from gene where stable_id='%s'" % stable_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return 0

    return int(rows[0][0])
    
########
def gene2stable (cursor, gene_id=None, db_name=None, ):

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

########
def get_compara_name (cursor):

    qry = "show databases like '%compara%'"
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return rows[0][0]

########
def species2taxid (cursor, species):


    qry  = "select taxon_id from genome_db where name = '%s'" % species
    rows = search_db (cursor, qry)
    if (not rows):
        search_db (cursor, qry, verbose = True)
        return ""
    
    return rows[0][0]
