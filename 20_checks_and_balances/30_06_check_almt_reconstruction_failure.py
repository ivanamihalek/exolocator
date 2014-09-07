#!/usr/bin/python -u
# make the best alignment we can using the maps
# we currently have at hand

import pdb
#pdb.set_trace()

import MySQLdb, commands, re, os
import random, string
import inspect, pdb
from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.utils   import  *
from el_utils.el_specific   import  *
from el_utils.map     import  Map, get_maps, map2exon
from el_utils.tree    import  species_sort
from el_utils.ncbi    import  *
from el_utils.almt_cmd_generator import AlignmentCommandGenerator
from el_utils.config_reader      import ConfigurationReader
from el_utils.translation        import phase2offset, translation_bounds, crop_dna, translate
from el_utils.special_gene_sets  import *
from el_utils.processes import parallelize
from el_utils.exon_boundary_hacks import *
from bitstring import Bits
from random  import choice
# BioPython
from Bio          import  SeqIO
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

import pdb

verbose = True

#########################################
def make_exon_alignments(cursor, ensembl_db_name, canonical_human_exons,
                         mitochondrial, min_similarity, flank_length):
    alnmt_pep = {}
    alnmt_dna = {}
    first_human_exon = True
    for human_exon in canonical_human_exons:
        # make_exon_alignment defined in el_utils/el_specific.py
        # we need the info about the first human exon
        # to get a bit more lenient whent the first exon consists of M only
        # (I'll have to take a look at one point what's up with that)
        [alnmt_pep[human_exon], alnmt_dna[human_exon]]  =   make_exon_alignment(cursor, ensembl_db_name,  human_exon.exon_id, 
                                                                                human_exon.is_known,  mitochondrial, min_similarity, 
                                                                                flank_length, first_human_exon)   
        first_human_exon = False

    return [alnmt_pep, alnmt_dna] 


#########################################
def main():
    local_db   = False
    
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    print "using all protein coding genes"
    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
    min_similarity = cfg.get_value('min_accptbl_exon_sim') 
    flank_length = 10
    gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        
    new_afas = 0
    old_afas = 0
    ancient_afas = 0

    failed_afas = []
    for gene_id in gene_list:
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)
        if  check_afa_age (cfg, stable_id, max_days=30) == "new": 
            new_afas += 1
            continue                               
        elif  check_afa_age (cfg, stable_id, max_days=300) == "new": 
            old_afas += 1
            failed_afas.append(gene_id)
            continue                               
        elif  check_afa_age (cfg, stable_id, max_days=1000) == "new": 
           ancient_afas += 1
           failed_afas.append(gene_id)
           continue                               
            
    no_exons  = 0
    cases_with_no_orthos = 0
    no_exon_ids = []
    for gene_id in failed_afas:

        if ( (failed_afas.index(gene_id))%10 == 0 ): 
            print failed_afas.index(gene_id), "out of ", len(failed_afas)

        canonical_human_exons = get_canonical_coding_exons (cursor, gene_id, ensembl_db_name['homo_sapiens'])

        if not canonical_human_exons: 
            no_exon_ids.append(gene_id)
            no_exons += 1
            continue

        # reconstruct  per-exon alignments with orthologues
        mitochondrial = is_mitochondrial(cursor, gene_id)
        [alnmt_pep, alnmt_dna] = make_exon_alignments(cursor, ensembl_db_name, canonical_human_exons,
                                                      mitochondrial, min_similarity, flank_length)

        no_orthos = True
        for human_exon, almt in alnmt_pep.iteritems():
            if ( type(almt) is str or len(almt.keys()) >= 1): 
                no_orthos = False
                break

        if no_orthos:
            cases_with_no_orthos += 1
            continue
         
    
    print
    print "total genes", len(gene_list)
    print "new  afas", new_afas
    print "old  afas", old_afas
    print "ancient afas", ancient_afas

    print
    print "failure cases"
    print "\t no exons", no_exons
    print "\t no orthologues ", cases_with_no_orthos
    
    print
    for gene_id in no_exon_ids:
        print gene_id
        for exon in gene2exon_list(cursor, gene_id):
            print "\t", exon.is_canonical, exon.is_coding

    cursor.close()
    db.close()



#########################################
if __name__ == '__main__':
    main()
