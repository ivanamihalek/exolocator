#!/usr/bin/python


import MySQLdb, os
from   el_utils.mysql         import  connect_to_mysql, connect_to_db, search_db
from   el_utils.mysql         import  store_or_update,  switch_to_db
from   el_utils.ensembl       import  get_species, get_gene_ids, gene2stable, stable2gene
from   el_utils.ensembl       import  get_compara_name
from   el_utils.threads       import  parallelize
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.utils         import  erropen



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
def read_paralogues(cursor, gene_id):

    paralogues = []
    qry = "select gene.stable_id from gene, paralogue "
    qry += " where gene.gene_id = paralogue.cognate_gene_id "
    qry += " and paralogue.gene_id = %d " % gene_id
    rows = search_db (cursor, qry)


    if rows and not 'ERROR' in rows[0]:
        for row in rows:
            paralogues.append(row[0])
    
    return paralogues
        


#########################################
def dump_paralogues(species_list, db_info):
    
    [local_db, ensembl_db_name] = db_info

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    #for species in species_list:
    for species in ['homo_sapiens']:
        print
        print "############################"
        print  species
        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)
        
        outfile  = "{0}/{1}_para_dump.txt".format(cfg.dir_path['afs_dumps'], species)
        print outfile
        #if (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
        #    continue
        of       = erropen (outfile,"w")
        
        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        para_table = 'paralogue'

        ct   =  0
        seen = []
        for gene_id in [378128, 398993] + gene_ids:
            ct += 1
            if not ct%10: print "\t", ct, "out of", len(gene_ids)

            if gene_id in seen: continue

            stable_id  = gene2stable(cursor, gene_id)
                
            paralogues = read_paralogues(cursor, gene_id)
            
            if ( paralogues):
                # dump
                for para in paralogues:
                    print >> of,  stable_id, para
                seen += paralogues
                exit(1)
        of.close()
 
    cursor.close()
    db.close()

#########################################
def main():
    
    no_threads = 1
    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    all_species = ['bos_taurus', 'callithrix_jacchus']

    cursor.close()
    db    .close()

    parallelize (no_threads, dump_paralogues, all_species, [local_db, ensembl_db_name])

    
    return True


#########################################
if __name__ == '__main__':
    main()
