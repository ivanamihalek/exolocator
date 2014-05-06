#!/usr/bin/python

import MySQLdb
import commands, pdb
from   el_utils.mysql   import  *  
from   el_utils.ensembl import  * 
from   el_utils.utils   import  *

#########################################
def main():

    local_db   = False
    number_of_lists = 10

    if local_db:
        db  = connect_to_mysql()
     else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    gene_list = [] 
    for n in range(number_of_lists):
        gene_list.push([]) # make number_of_lists empty lists
        
    # distribute in round robin way
    for gene_id in gene_ids[:117]:
        stable_id = gene2stable (cursor, gene_id = gene_id)
        n = gene_ids.id(gene_id)%number_of_lists  # number (or index) of this gene_id in the whole list
        gene_list[n].push(stable_id)

    for n in range(number_of_lists):
        print
        print "list", n
        for stable_id in gene_list:
            print stable_id

    cursor.close()
    db    .close()



#########################################
if __name__ == '__main__':
    main()


