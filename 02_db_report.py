#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql
from   el_utils.ensembl import  get_species, get_gene_ids


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)


    for species in all_species:
        print species
        gene_ids = get_gene_ids (cursor, ensembl_db_name[species], 'protein_coding')  
        print " protein coding genes:  %15d " %  len(gene_ids)
    print "there are %d core dbs available " % len(all_species)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
