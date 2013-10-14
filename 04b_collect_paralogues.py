#!/usr/bin/python
# why is git not pushing this?

import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, connect_to_db, search_db
from   el_utils.mysql   import  store_or_update,  switch_to_db
from   el_utils.ensembl import  *
from   el_utils.threads import  parallelize

#########################################
def store_paralogues (cursor_species,  gene_id, orthos):

    for ortho in orthos:
        [ortho_stable, species, cognate_genome_db_id] = ortho
        ortho_gene_id = stable2gene (cursor_species, ortho_stable)
        
        fixed_fields  = {}
        fixed_fields ['gene_id']              = gene_id
        fixed_fields ['cognate_genome_db_id'] = cognate_genome_db_id
        fixed_fields ['cognate_gene_id']      = ortho_gene_id

        update_fields = {}
        update_fields['source']          = 'ensembl'
                    
        store_or_update (cursor_species, 'paralogue', fixed_fields, update_fields)


#########################################
def collect_paralogues(species_list, db_info):
    
    [local_db, ensembl_db_name] = db_info


    if local_db:
        db_species  = connect_to_mysql()
    else:
        db_species  = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor_species = db_species.cursor()

    ensembl_compara_name = get_compara_name(cursor_species)
    print ensembl_compara_name
 
    if local_db:
        db_compara     = connect_to_mysql()
    else:
        db_compara     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor_compara     = db_compara.cursor()
    switch_to_db (cursor_compara, ensembl_compara_name)

    for species in species_list:
        switch_to_db (cursor_species,  ensembl_db_name[species])
        gene_list = get_gene_ids (cursor_species, biotype='protein_coding', is_known=1)
        ct = 0
        for gene_id in gene_list:
            ct += 1
            # find stable
            stable_id = gene2stable(cursor_species, gene_id=gene_id)
            # memebr id refers to entries in compara db
            member_id = stable2member(cursor_compara, stable_id)

            #print gene_id, stable_id, member_id
            if ( not ct%10):
                print species, ct , "out of ", len(gene_list) 
            # find all paralogue pairs suggested for this gene
            ortho_type = 'within_species_paralog'
            paralogues = get_orthologues(cursor_compara, ortho_type, member_id)
            if not paralogues: continue
            store_paralogues (cursor_species, gene_id, paralogues)
        
    cursor_species.close()
    db_species.close()
    cursor_compara.close()
    db_compara.close()

#########################################
def main():
    
    no_threads = 5
    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db.close()

    parallelize (no_threads, collect_paralogues, all_species, [local_db, ensembl_db_name])
    
    return True


#########################################
if __name__ == '__main__':
    main()
