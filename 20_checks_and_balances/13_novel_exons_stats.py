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
    local_db = False
    

    if local_db:
        db  = connect_to_mysql()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    for species in all_species:
        print species
        switch_to_db (cursor, ensembl_db_name[species])
        qry  = "select count(1) from usearch_exon"
        rows = search_db (cursor, qry)
        print "\t usearch exons: ", rows[0]

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

