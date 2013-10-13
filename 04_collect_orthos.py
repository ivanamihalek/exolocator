#!/usr/bin/python
# why is git not pushing this?

import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, connect_to_db, search_db
from   el_utils.mysql   import  store_or_update,  switch_to_db
from   el_utils.ensembl import  *
from   el_utils.threads import  parallelize

#########################################
def store_orthologues (cursor_human, ortho_table, cursor, all_species, 
                       ensembl_db_name,  gene_id, orthos):

    for ortho in orthos:
        [ortho_stable, species, cognate_genome_db_id] = ortho
        if (not species in all_species):
            continue
        ortho_gene_id = stable2gene (cursor, ortho_stable, ensembl_db_name[species])
        
        fixed_fields  = {}
        fixed_fields ['gene_id']              = gene_id
        fixed_fields ['cognate_genome_db_id'] = cognate_genome_db_id

        update_fields = {}
        update_fields['source']          = 'ensembl'
        
        if ( ortho_table == 'orthologue'):
            update_fields['cognate_gene_id'] = ortho_gene_id
        else:
            fixed_fields['cognate_gene_id'] = ortho_gene_id
            

        store_or_update (cursor_human, ortho_table, fixed_fields, update_fields)



#########################################
def collect_orthologues(gene_list, db_info):
    
    [local_db, ensembl_db_name] = db_info

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)


    if local_db:
        db_human  = connect_to_mysql()
    else:
        db_human  = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor_human = db_human.cursor()
    switch_to_db (cursor_human,  ensembl_db_name['homo_sapiens'])

    ensembl_compara_name = get_compara_name(cursor)
    print ensembl_compara_name
    exit(1);

    if local_db:
        db_compara     = connect_to_mysql()
    else:
        db_compara     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor_compara     = db_compara.cursor()
    switch_to_db (cursor_compara, ensembl_compara_name)


    ortho_table = {}
    ortho_table ['ortholog_one2one']          = 'orthologue'
    ortho_table ['possible_ortholog']         = 'orthologue'
    ortho_table ['apparent_ortholog_one2one'] = 'orthologue'
    ortho_table ['ortholog_one2many']         = 'unresolved_ortho'
    ortho_table ['ortholog_many2many']        = 'unresolved_ortho'
    ortho_table ['within_species_paralog']    = 'paralogue'
    ct = 0
    for gene_id in gene_list:
        ct += 1
        # find stable
        stable_id = gene2stable(cursor_human, gene_id=gene_id)
        # memebr id refers to entries in compara db
        member_id = stable2member(cursor_compara, stable_id)

        #print gene_id, stable_id, member_id
        if ( not ct%10):
            print ct , "out of ", len(gene_list)
        # in compara table, get everything that homology has to say about
        # the possible orthologues
        # find all orthologous pairs suggested for this gene
        for ortho_type in ['ortholog_one2one','possible_ortholog', 'apparent_ortholog_one2one',
                           'ortholog_one2many','ortholog_many2many','within_species_paralog']:
            orthos = get_orthologues(cursor_compara, ortho_type, member_id)

            if ( orthos):
                store_orthologues (cursor_human, ortho_table[ortho_type], cursor, all_species, 
                                   ensembl_db_name, gene_id, orthos)
        
 
    cursor.close()
    db.close()
    cursor_human.close()
    db_human.close()
    cursor_compara.close()
    db_compara.close()

#########################################
def main():
    
    no_threads = 1
    local_db = True

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    species                        = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_list                      = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
    cursor.close()
    db.close()

    parallelize (no_threads, collect_orthologues, gene_list, [local_db, ensembl_db_name])
    
    return True


#########################################
if __name__ == '__main__':
    main()
