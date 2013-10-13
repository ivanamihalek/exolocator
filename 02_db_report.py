#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  get_species, get_gene_ids


#########################################
def check_genome_sizes (cursor, all_species, ensembl_db_name):
    for species in all_species:
        print species
        gene_ids = get_gene_ids (cursor, ensembl_db_name[species], 'protein_coding')  
        print " protein coding genes:  %15d " %  len(gene_ids)
    print "there are %d core dbs available " % len(all_species)
  

#########################################
def check_table_sizes (cursor, all_species, ensembl_db_name):
    for species in all_species:
        print
        print "##########################"
        print species
        qry  = "use "+ensembl_db_name[species]
        search_db(cursor, qry)
        qry  = "show tables"
        rows = search_db(cursor, qry)
        for row in rows:
            table = row[0]
            qry = " select count(1) from "+table
            rows = search_db(cursor, qry)
            table_size = rows[0][0]
            print "\t ", table, table_size


#########################################
def main():
    
    #db    = connect_to_mysql()
    db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    if 0:
        check_genome_sizes (cursor, all_species, ensembl_db_name)

    if 1:
        check_table_sizes (cursor, all_species, ensembl_db_name)
      

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
