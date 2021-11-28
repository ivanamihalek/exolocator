#!/usr/bin/python
# why is git not pushing this?

import MySQLdb
from   el_utils.el_specific import  *
########################################
def main():

    db  = connect_to_mysql()
    # acg = AlignmentCommandGenerator()
    #
    # cursor = db.cursor()
    # [all_species, ensembl_db_name] = get_species (cursor)
    # species = 'homo_sapiens'
    # switch_to_db (cursor,  ensembl_db_name[species])
    # gene_list = get_gene_ids (cursor,  biotype='protein_coding')
    #
    # for gene_id in gene_list:
    #     # find stable
    #     stable_id = gene2stable(cursor,  gene_id=gene_id)
    #     canonical = get_canonical_transl (acg, cursor, gene_id, species, strip_X = False)
    #     if canonical:
    #         print stable_id, canonical
        
#    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
