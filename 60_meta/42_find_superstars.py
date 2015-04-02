#!/usr/bin/python

import MySQLdb
import commands
from   el_utils.mysql   import  connect_to_mysql, search_db, store_or_update
from   el_utils.ensembl import  *

#########################################
def gene_name2gene_id(cursor, gene_name):
    
    qry  = "select ensembl_id from object_xref, external_synonym "
    qry += "where object_xref.ensembl_object_type = 'Gene'  "
    qry += "and object_xref.xref_id= external_synonym.xref_id "
    qry += "and external_synonym.synonym = '%s' " % gene_name
    qry += "group by synonym"
    rows = search_db (cursor, qry)
    if not rows:
        return ""
    else:
        return rows[0][0]
#########################################
def search_description (cursor, gene_name):
    qry  = "select gene_id, description from gene "
    qry += "where description like '%"+gene_name+"%'"
    rows = search_db (cursor, qry)
    if not rows:
        return ["", ""]
    else:
        return rows[0]


#########################################
#########################################
def main():

    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    switch_to_db(cursor, ensembl_db_name['homo_sapiens'])

    magical_list = ['APC', 'BUB1', 'BUB1B', 'BUB3', 'C11orf51', 
                    'CDC20', 'CDC27', 'CENPF', 'TERF1', 'TPR', 
                    'TTK', 'UBE2C', 'UBE2D1', 'UBE2E1', 'TP53', 'BCL',
                    ' RAS', ' MIC ', 'actin']
    for gene_name in magical_list:
        description = ""
        gene_id = gene_name2gene_id(cursor, gene_name)
        if (not gene_id): 
            [gene_id, description] = search_description (cursor, gene_name)
        if (not gene_id): continue

        print gene_name, " ** ",  gene_id, description


    cursor.close()
    db    .close()
    



#########################################
if __name__ == '__main__':
    main()
