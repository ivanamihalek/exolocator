#!/usr/bin/python


import StringIO
import MySQLdb, commands, re, sys
from hashlib import sha1
from random  import random
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.utils   import  *
from   el_utils.threads import  parallelize
from   el_utils.map     import  *
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.special_gene_sets  import *


    
#########################################
def store (cursor, maps, ensembl_db_name, source = None):

 
    for map in maps:
        fixed_fields  = {}
        update_fields = {}
        fixed_fields ['exon_id']              = map.exon_id_1
        fixed_fields ['exon_known']           = map.exon_known_1
        fixed_fields ['cognate_genome_db_id'] = species2genome_db_id(cursor, map.species_2)
        fixed_fields ['cognate_exon_id']      = map.exon_id_2
        fixed_fields ['cognate_exon_known']   = map.exon_known_2
        fixed_fields ['source']               = map.source 
        update_fields['cigar_line']           = map.cigar_line
        update_fields['similarity']           = map.similarity
        ################################
        switch_to_db(cursor,ensembl_db_name['homo_sapiens']) 
        store_or_update (cursor, 'exon_map', fixed_fields, update_fields)

    return True

#########################################
def  map_cleanup (cursor, ensembl_db_name, human_exons):
    
    switch_to_db(cursor,ensembl_db_name['homo_sapiens']) 
    for exon in human_exons:
        qry  = "delete from exon_map where exon_id = %d " % exon.exon_id
        qry += " and exon_known = %d " % exon.is_known
        qry += " and cognate_exon_known > 1 " 
        rows = search_db (cursor, qry, verbose=False)


    return True


#########################################
def gene_has_a_map (cursor, ensembl_db_name, human_exons):

    has_a_map = False
    for human_exon in human_exons:
        if ( not human_exon.is_canonical or  not human_exon.is_coding): continue
        maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
        if maps:
            has_a_map = True
            break

    return has_a_map

#########################################
def maps_for_gene_list(gene_list, db_info):
    

    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        cfg    = ConfigurationReader()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    missing_exon_info = 0
    missing_seq_info  = 0
    ct                = 0
    no_maps           = 0

    #######################################
    for gene_id in gene_list:

        ct += 1
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        if not ct%10: print ct, "out of ", len(gene_list) 
        
        # get _all_ exons
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        print  gene_id,  gene2stable(cursor, gene_id), get_description (cursor, gene_id)

        human_exons = gene2exon_list(cursor, gene_id)
        if (not human_exons):
            print 'no exons for ', gene_id
            continue

        # get rid of the old maps # can't do that here bcs this script is only updating sw exons
        # map_cleanup (cursor, ensembl_db_name, human_exons)

        orthologues  = get_orthos (cursor, gene_id, 'orthologue') # get_orthos changes the db pointer
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        orthologues += get_orthos (cursor, gene_id, 'unresolved_ortho')

        ##########
        for [ortho_gene_id, ortho_species] in orthologues:

            #if not ortho_species=='ochotona_princeps': continue
            #print ortho_species, ortho_gene_id, species2genome_db_id (cursor, ortho_species)

            switch_to_db (cursor, ensembl_db_name[ortho_species])

            ortho_exons   = []
            
            ortho_exons += get_novel_exons (cursor, ortho_gene_id, 'sw_exon')
            ortho_exons += get_novel_exons (cursor, ortho_gene_id, 'usearch_exon')

            if not ortho_exons: continue # nothing new here, move on

            ortho_exons +=  get_known_exons (cursor, ortho_gene_id, ortho_species)
            ortho_exons +=  get_predicted_exons (cursor, ortho_gene_id, ortho_species)

            maps = make_maps (cursor, ensembl_db_name, cfg, acg, ortho_species, human_exons, ortho_exons) 

            if not maps:
                print "\t", ortho_species, "no maps"
                continue

            switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
            store (cursor, maps, ensembl_db_name)
                
    cursor.close()
    db.close()

    return True


#########################################
def main():
    
    no_threads = 1
    special    = 'one'

    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <set name> <number of threads> <method>"
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
        print "using all protein coding genes that have an sw# or usearch patch"
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        gene_list = human_genes_w_novel_exon_orthologues (cursor, ensembl_db_name)

    cursor.close()
    db.close()

    parallelize (no_threads, maps_for_gene_list, gene_list[0:15000], [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
    main()

'''

    for gene_id in [412667]: #  wls
    for gene_id in [378768]: #  p53
     #for gene_id in [378766]: #  dynein
 

        # COMMENT THIS OUT PERHAPS?
        # get rid of the old maps
        # map_cleanup (cursor, ensembl_db_name, human_exons)
        

'''
