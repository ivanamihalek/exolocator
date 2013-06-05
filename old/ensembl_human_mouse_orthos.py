#!/usr/bin/python


import StringIO
import MySQLdb, commands, re, sys
from hashlib import sha1
from random  import random
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  *
from   el_utils.utils   import  *
from   el_utils.threads import  parallelize
from   el_utils.map     import  get_maps, Map

from   el_utils.special_gene_sets  import  get_theme_ids
from   el_utils.almt_cmd_generator import  AlignmentCommandGenerator
from   el_utils.config_reader      import  ConfigurationReader
from   el_utils.alignment          import  smith_waterman, exon_aware_smith_waterman
from   alignment import * # C implementation of smith waterman



#########################################
def main():

    local_db=False

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

    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
    gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        
    for gene_id in gene_list:
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        orthologues  = get_orthos (cursor, gene_id, 'orthologue') # get_orthos changes the db pointer
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        orthologues += get_orthos (cursor, gene_id, 'unresolved_ortho')
        human_stable = gene2stable(cursor, gene_id, ensembl_db_name['homo_sapiens'])
        human_descr = get_description (cursor, gene_id)
        for [ortho_id, species] in orthologues:
            if not species == 'mus_musculus':
                continue
            mouse_stable =  gene2stable(cursor, ortho_id, ensembl_db_name['mus_musculus'])
            if not mouse_stable: continue
            if not human_descr: human_descr=''
            print  "\t".join([human_stable, mouse_stable, human_descr])
 
    cursor.close()
    db.close()



#########################################
if __name__ == '__main__':
    main()
