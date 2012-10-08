#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  get_species, get_gene_ids
from   el_utils.ojects  import  Exon
from   el_utils.threads import  parallelize


#########################################
def gene2exon():

    print "hello?"

    return True

#########################################
def main():

    no_threads = 1

    
    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    parallelize (no_threads, gene2exon, all_species)


    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
