#!/usr/bin/python

import MySQLdb
from el_utils.mysql   import  *
from el_utils.ensembl import  *

#########################################
def main():
    
    db  = connect_to_mysql()

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    species = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1, ref_only=True)
    for gene_id in gene_ids:
        if not is_mitochondrial(cursor, gene_id): continue
        print gene2stable(cursor, gene_id), get_description(cursor, gene_id)
        

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
