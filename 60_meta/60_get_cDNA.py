#!/usr/bin/python

import MySQLdb
import commands
from   random               import choice
from   el_utils.mysql       import  connect_to_mysql, search_db
from   el_utils.ensembl     import  *
from   el_utils.exon        import  Exon


#########################################
def main():

    local_db = False

    db     = connect_to_mysql()
    acg    = AlignmentCommandGenerator()

    cursor = db.cursor()
    
    [all_species, ensembl_db_name] = get_species (cursor)
    
    species = 'homo_sapiens'

    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    for gene_id in gene_ids:
        
        print gene2stable (cursor, gene_id = gene_id),

        # what is the length of the canonical transcript according to Ensembl
        canonical_translation = get_canonical_transl (acg, cursor, gene_id, species, strip_X=False)
        if ( not canonical_translation):
            print "no canonical transl found for ", gene2stable (cursor, gene_id = gene_id)
            continue

        # find all canonical coding exons associated with the gene id
        exons = get_canonical_coding_exons (cursor, gene_id)
        if (not exons):
            ct +=1
            print gene_id, gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot

        exit(1)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
