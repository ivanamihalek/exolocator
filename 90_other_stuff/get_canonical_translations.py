#!/usr/bin/python
# why is git not pushing this?

import MySQLdb
from   el_utils.el_specific import  *
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator

#########################################
def main():

    local_db = False
    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    species = 'homo_sapiens'
    switch_to_db (cursor_species,  ensembl_db_name[species])
    gene_list = get_gene_ids (cursor_species, biotype='protein_coding')

    for gene_id in gene_list:
        # find stable
        stable_id = gene2stable(cursor_species, gene_id=gene_id)
        canonical = get_canonical_transl (acg, cursor, gene_id, species, strip_X = False)
        if canonical:
            print stable_id, canonical
        exit(1)
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
