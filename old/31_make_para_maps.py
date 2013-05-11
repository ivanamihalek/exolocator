#!/usr/bin/python


import MySQLdb, os
from   el_utils.mysql         import  connect_to_mysql, connect_to_db, search_db
from   el_utils.mysql         import  store_or_update,  switch_to_db
from   el_utils.ensembl       import  *
from   el_utils.threads       import  parallelize
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.utils         import  erropen
from   el_utils.custom        import  get_theme_ids


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
def member2homology (cursor, member_id):

    qry  = "select homology.homology_id from homology_member, homology "
    qry += " where homology_member.member_id =%d " % member_id
    qry += " and homology.homology_id = homology_member.homology_id "
    qry += " and (homology.description = 'within_species_paralog' "
    qry += " or homology.description   = 'other_paralog')"
    rows = search_db (cursor, qry)

    if not rows: return []

    hom_ids = []
    for row in rows:
        hom_ids.append(row[0])

    return hom_ids

################ 
def  homology2other_member (cursor, homology_id, member_id):

    qry  = "select member_id from homology_member "
    qry += " where homology_id = %d"  % int(homology_id)
    qry += " and not  member_id = %d" % member_id
    rows = search_db (cursor, qry)

    if (not rows):
        return "" # no orthologs here

 
    return rows[0][0]

################
def get_paralogues(cursor,  member_id):

    paralogues = [member_id]
    prev_new   = [member_id]
    done = False
    while not done:

        new = []

        for member_id in prev_new:
            # for each homology id find the other member id
            homology_ids  = member2homology(cursor, member_id)
            for homology_id in homology_ids:
                other_member = homology2other_member (cursor, homology_id, member_id)
                if not other_member in paralogues and not other_member in new:
                    new.append(other_member)
        if new:
            paralogues += new
            prev_new    = new
        else:
            done = True

    paralogues_stable = []
    for member_id in paralogues[1:]:
        paralogues_stable.append(member2stable(cursor, member_id))

    return paralogues_stable
        
#########################################
def store_paralogues (cursor,  gene_id, paralogues):

    for para_stable in paralogues:
        

        para_gene_id = stable2gene (cursor, para_stable)
        
        fixed_fields  = {}
        fixed_fields ['gene_id']              = gene_id
        fixed_fields[ 'cognate_gene_id'] = para_gene_id

        update_fields = {}
        update_fields['source']          = 'ensembl'
        
            
        store_or_update (cursor, 'paralogue', fixed_fields, update_fields)



#########################################
def make_orthologue_table (cursor, table):
    

    # if congamet_gene_id is 0, and source is 'rbh'
    # means that the reciprocal-best-hit was attempted but nothing was found

    qry  = "CREATE TABLE " + table + "  (orth_pair_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    for column_name in ['gene_id', 'cognate_gene_id']:
        qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    for column_name in ['cognate_genome_db_id']:
        qry = "ALTER TABLE %s  ADD %s INT" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

    
    for column_name in ['source']:
        qry = "ALTER TABLE %s  ADD %s VARCHAR(20)" % (table, column_name)
        rows = search_db (cursor, qry)
        if (rows):
            return False

#########################################
def collect_paralogues(species_list, db_info):
    

    verbose = False

    [local_db, ensembl_db_name] = db_info

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor  = db.cursor()

    ensembl_compara_name = get_compara_name(cursor)
    print ensembl_compara_name

    if local_db:
        db_compara = connect_to_mysql()
    else:
        db_compara = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor_compara = db_compara.cursor()

    switch_to_db (cursor_compara, ensembl_compara_name)

    for species in species_list:
        if species in ['homo_sapiens']: continue
        print
        print "############################"
        print  species
        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)
        para_table = 'paralogue'
        
        #qry = "drop table " + para_table
        #search_db(cursor, qry)
        #make_orthologue_table (cursor, para_table)


        gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        #gene_ids = get_theme_ids(cursor, cfg, 'wnt_pathway')

        gene_ct =  0
        seen    = []
        for gene_id in  gene_ids:
        #for gene_id in [378128, 398993]: # wnt6, wnt8
        #for gene_id in [378128, 398993]:
            gene_ct += 1
            if verbose:         
                stable_id = gene2stable(cursor, gene_id)
                print stable_id, get_description (cursor, gene_id)
            elif (not gene_ct%100): 
                print species, gene_ct, "out of ", len(gene_ids)

            # find stable
            stable_id = gene2stable(cursor, gene_id=gene_id)
            # memeber id refers to entries in compara db
            member_id = stable2member(cursor_compara, stable_id)

            if stable_id in seen: continue
                
            paralogues = get_paralogues(cursor_compara, member_id)

            #for para_stable in paralogues:
            #    para_gene_id = stable2gene (cursor, para_stable)
            #    print para_stable, get_description (cursor, para_gene_id)

            if (paralogues):
                store_paralogues (cursor,  gene_id, paralogues)
                seen += paralogues

    cursor.close()
    db.close()
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
    db    .close()

    parallelize (no_threads, collect_paralogues, all_species, [local_db, ensembl_db_name])

    
    return True


#########################################
if __name__ == '__main__':
    main()


'''
    # for each homology id find the other member id
    panic_ct = 0
    seen_homology_pairs = []
    while paralogues_prev:
        paralogues_new = []
        for  member_id in paralogues_prev:
            
            homology_ids  = member2homology(cursor, member_id)

            for homology_id in homology_ids:
                if homology_id in seen_homology_pairs: continue
                other_member = homology2other_member (cursor, homology_id, member_id)
                if not other_member in paralogues:
                    paralogues_new.append(other_member)
                seen_homology_pairs.append(homology_id)
        
        print " *  ", paralogues_prev
        print " ** ", paralogues_new
        paralogues     += paralogues_new
        paralogues_prev = paralogues_new

        panic_ct += 1
        if (panic_ct == 10000):
            return []

 '''
