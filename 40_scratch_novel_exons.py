#!/usr/bin/python


import MySQLdb
from   el_utils.mysql       import  *
from   el_utils.ensembl     import  *
from el_utils.config_reader import  ConfigurationReader
from el_utils.special_gene_sets  import get_theme_ids, get_complement_ids
from el_utils.threads import  parallelize
from subprocess import Popen, PIPE, STDOUT

#########################################
def delete_novel_exons (human_gene_list, db_info):

    # 
    [local_db, ensembl_db_name] = db_info
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    # find db ids and common names for each species db
    all_species, ensembl_db_name = get_species (cursor)

    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])

    for human_gene_id in human_gene_list:

	human_exons = [e for e in gene2exon_list(cursor, human_gene_id) 
                       if e.covering_exon < 0 and e.is_canonical and e.is_known]

        for species in all_species:
            if species=='homo_sapiens': continue
            print species
            cognate_genome_db_id = species2genome_db_id(cursor, species) # moves the cursor!

            ###################################
            # usearch exons:
            db_name = ensembl_db_name[species]
            switch_to_db (cursor, ensembl_db_name[species])

            exons_to_delete = []
            for he in human_exons:
                qry = "select exon_id from usearch_exon"
                qry += " where maps_to_human_exon_id=%d " % he.exon_id
                rows = search_db(cursor, qry)
                if rows: 
                    exons_to_delete.append(int(rows[0][0]))

            for e in exons_to_delete:
                qry = "delete from  usearch_exon where exon_id=%d" % e
                search_db(cursor, qry)
                qry = "delete from exon_seq  where is_sw=2 and exon_id=%d" % e
                search_db(cursor, qry)

                qry = "delete from gene2exon where is_known=3 and exon_id=%d" % e
                search_db(cursor, qry)

            switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
            for e in exons_to_delete:
                qry  = "delete from exon_map  where cognate_exon_id=%d" % e
                qry += " and  cognate_exon_known=3 and cognate_genome_db_id = %d " % cognate_genome_db_id
                search_db(cursor, qry)

   
            ###################################
            # sw_sharp exons:
            db_name = ensembl_db_name[species]
            switch_to_db (cursor, ensembl_db_name[species])

            exons_to_delete = []
            for he in human_exons:
                qry = "select exon_id from sw_exon"
                qry += " where maps_to_human_exon_id=%d " % he.exon_id
                rows = search_db(cursor, qry)
                if rows: 
                    exons_to_delete.append(int(rows[0][0]))

            for e in exons_to_delete:
                qry = "delete from  usearch_exon where exon_id=%d" % e
                search_db(cursor, qry)
                qry = "delete from exon_seq  where is_sw=1 and exon_id=%d" % e
                search_db(cursor, qry)

                qry = "delete from gene2exon where is_known=2 and exon_id=%d" % e
                search_db(cursor, qry)

            switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
            for e in exons_to_delete:
                qry  = "delete from exon_map  where cognate_exon_id=%d" % e
                qry += " and  cognate_exon_known=2 and cognate_genome_db_id = %d " % cognate_genome_db_id
                search_db(cursor, qry)

   


#########################################
def main():

    special    = 'test'
    no_threads = 1


    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <set name> <number of threads> " % sys.argv[0]
        exit(1)
    elif len(sys.argv)==3:

        special = sys.argv[1]
        special = special.lower()
        if special == 'none': special = None

        no_threads = int(sys.argv[2])
        
 

    local_db   = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)


    print '======================================='
    print sys.argv[0]
    if special:
        print "using", special, "set"
        if special == 'complement':
            gene_list = get_complement_ids(cursor, ensembl_db_name, cfg)
        else:
            gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special )

    else:
        print "using all protein coding genes"
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    cursor.close()
    db.close()

    parallelize (no_threads, delete_novel_exons, gene_list, [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
    main()
