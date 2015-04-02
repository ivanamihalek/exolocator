#!/usr/bin/python

import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl import *
from el_utils.special_gene_sets  import get_theme_ids, get_complement_ids
from el_utils.utils       import *



#########################################
def main():
    
    special    = None
    no_threads = 1
    db  = connect_to_mysql()

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    total = 0
    for species in all_species:
        print species
        switch_to_db (cursor, ensembl_db_name[species])
        qry  = "select count(1) from usearch_exon"
        rows = search_db (cursor, qry)
        count = int(rows[0][0])
        print "\t usearch exons: ", count 
        total += count
        qry  = "select count(1) from sw_exon"
        rows = search_db (cursor, qry)
        count = int(rows[0][0])
        print "\t sw exons: ", count 
        total += count
    print
    print 'total: ', total
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

