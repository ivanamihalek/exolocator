#!/usr/bin/python

import MySQLdb
import subprocess, pdb
from   el_utils.mysql   import  *  
from   el_utils.ensembl import  * 
from   el_utils.utils   import  *

#########################################
def main():

    number_of_lists = 10

    db     = connect_to_mysql()
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    gene_list = [] 
    for n in range(number_of_lists):
        gene_list.append([]) # make number_of_lists empty lists
        
    # distribute in round robin way
    for gene_id in gene_ids:
        stable_id = gene2stable (cursor, gene_id = gene_id)
        n = gene_ids.index(gene_id)%number_of_lists  # number (or index) of this gene_id in the whole list
        gene_list[n].append(stable_id)

    for n in range(number_of_lists):
        outf = open ("stable_ids_"+str(n)+".txt", "w")
        for stable_id in gene_list[n]:
            print(stable_id, file=outf)
        outf.close()

    cursor.close()
    db    .close()



#########################################
if __name__ == '__main__':
    main()


