#!/usr/bin/python
# why is git not pushing this?

import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, connect_to_db, search_db
from   el_utils.mysql   import  store_or_update,  switch_to_db
from   el_utils.ensembl import  get_species, get_gene_ids, gene2stable, stable2gene
from   el_utils.ensembl import  get_compara_name
from   el_utils.threads import  parallelize



########
def stable2member (cursor, stable_id):
    
    # member_id refers to compara db
    # of which we need to have one
    qry = "select  member_id from member where stable_id = '%s'" % stable_id
    rows = search_db (cursor, qry)
    if (not rows or 'ERROR' in rows[0]):
        rows = search_db (cursor, qry, verbose = True)
        exit(1)
        return ""
    
    return int(rows[0][0])

########
def member2stable (cursor, member_id):
    
    # member_id refers to compara db
    # of which we need to have one
    qry = "select  stable_id from member where member_id = %d" % member_id
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return rows[0][0]



########
def get_orthologues(cursor, ortho_type, member_id):

    orthos = []

    qry  = "select homology.homology_id from homology_member, homology "
    qry += " where homology_member.member_id =%d " % member_id
    qry += " and homology.homology_id = homology_member.homology_id "
    qry += " and  homology.description = '%s' " % ortho_type
    rows = search_db (cursor, qry)

    if (not rows):
        return [] # no orthologs here

    # for each homology id find the other member id
    for row in rows:
        homology_id = row[0]

        qry  = "select member_id from homology_member "
        qry += " where homology_id = %d"  % int(homology_id)
        qry += " and not  member_id = %d" % member_id
        rows = search_db (cursor, qry, verbose = True)
        if (not rows):
            rows = search_db (cursor, qry, verbose = True)
            return []
        ortho_id     = rows[0][0]
        
        qry  = "select  member.stable_id, genome_db.name, genome_db.genome_db_id "
        qry += " from member, genome_db "
        qry += " where member.member_id = %d " % ortho_id
        qry += " and genome_db.genome_db_id = member.genome_db_id"
        rows = search_db (cursor, qry, verbose = True)
        if (not rows):
            rows = search_db (cursor, qry, verbose = True)
            return []
        [ortho_stable, species, genome_db_id]     = rows[0]
        orthos.append([ortho_stable, species,  int(genome_db_id)])
        
    return orthos
        
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

    #db_compara    = connect_to_mysql()
    ensembl_compara_name = get_compara_name(cursor)
    if local_db:
        db_compara     = connect_to_mysql()
    else:
        db_compara     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor_compara     = db_compara.cursor()
    switch_to_db (cursor_compara, ensembl_compara_name)
    print ensembl_compara_name


    ortho_table = {}
    ortho_table ['ortholog_one2one']       = 'orthologue'
    ortho_table ['ortholog_one2many']      = 'unresolved_ortho'
    ortho_table ['ortholog_many2many']     = 'unresolved_ortho'
    ortho_table ['within_species_paralog'] = 'paralogue'
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
        for ortho_type in ['ortholog_one2one','ortholog_one2many','ortholog_many2many','within_species_paralog']:
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
    
    no_threads = 5
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
