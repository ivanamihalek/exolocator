#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  get_species, get_gene_ids


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)


    species = 'homo_sapiens'
    gene_ids = get_gene_ids (cursor, ensembl_db_name[species], 'protein_coding')  
    print " protein coding genes in human:  %10d " %  len(gene_ids)

    for ortho_table in ['orthologue', 'unresolved_ortho', 'paralogue']:
        # how many orthologue pairs, in principle
        qry = "select count(1) from %s.%s" % (ensembl_db_name[species], ortho_table)
        rows = search_db (cursor, qry)
        if (not rows):
            rows = search_db (cursor, qry, verbose=True)
        print ortho_table, " table size: ", rows[0][0]

    # how many per gene
    for ortho_table in ['orthologue', 'unresolved_ortho', 'paralogue']:
        
        print "histogram for ", ortho_table

        histogram = {}
        for gene_id in gene_ids:

            qry  = "select count(1) from %s.%s" % (ensembl_db_name[species], ortho_table)
            qry += " where gene_id = %d " % gene_id
            rows = search_db (cursor, qry)
            if (not rows):
                rows = search_db (cursor, qry, verbose=True)
            if ( not histogram.has_key(rows[0][0])):
                histogram[rows[0][0]] = 0
            histogram[rows[0][0]] += 1
 
        for number_of_orthos in sorted (histogram.keys()):
            print " %4d  %4d " % (number_of_orthos, histogram[number_of_orthos])
        

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
