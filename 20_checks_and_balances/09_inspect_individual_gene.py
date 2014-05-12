#!/usr/bin/python

import MySQLdb
import sys, commands
from random import choice
from   el_utils.mysql       import  connect_to_mysql, search_db
from   el_utils.ensembl     import  *
from   el_utils.exon        import  Exon
from   el_utils.threads     import  parallelize
from   el_utils.almt_cmd_generator  import AlignmentCommandGenerator
from   el_utils.special_gene_sets   import  get_theme_ids
from   el_utils.config_reader       import ConfigurationReader

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#########################################
def main():

    if (len(sys.argv) < 3):
        print "Usage: %s <species> <stable gene id>" % sys.argv[0]
        exit(1)

    [species, stable_id] = sys.argv[1:3]


    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
  

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)


    print species, stable_id, ensembl_db_name[species]
  
    switch_to_db (cursor, ensembl_db_name[species])
    gene_id = stable2gene(cursor, stable_id)

    print get_description(cursor,gene_id)
    print "gene id:", gene_id

    # find all exons we are tracking in the database
    human_exons     = gene2exon_list(cursor, gene_id)
    canonical_human_exons = []
    for human_exon in human_exons:
        if not human_exon.is_canonical or  not human_exon.is_coding:
            continue
        canonical_human_exons.append(human_exon)

    # the exons are not guaranteed to be in order
    canonical_human_exons.sort(key=lambda exon: exon.start_in_gene)

    print "exons:"
    for exon in canonical_human_exons:
        exon_seqs = get_exon_seqs (cursor, exon.exon_id, 1)
        [exon_pep_seq, trsl_from, trsl_to, exon_left_flank,
         exon_right_flank, exon_dna_seq] = exon_seqs [1:]
        print "exon:", exon.exon_id, "covering exon:", exon.covering_exon,  "pepseq:", exon_pep_seq
        if  not exon.covering_exon == -1:
            [exon_pep_seq_2, trsl_from, trsl_to, exon_left_flank,
             exon_right_flank, exon_dna_seq] =  get_exon_seqs (cursor, exon.covering_exon, 1)[1:]
            print "\t", exon.covering_exon, " seq:", exon_pep_seq_2

    #canonical_translation = get_canonical_transl (acg, cursor, gene_id, species)
    #print canonical_translation

    cursor.close()
    db    .close()

    

#########################################
if __name__ == '__main__':
    main()

