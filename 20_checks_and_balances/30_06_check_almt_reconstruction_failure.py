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
from   random  import choice
# BioPython
from Bio          import  SeqIO
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

import pdb

verbose = True

from reconstruct_ortho_alnmts import *


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
            
    print "total genes", len(gene_list)
    print "new  afas", new_afas
    print "old  afas", old_afas
    print "ancient afas", ancient_afas

    no_exons  = 0
    no_orthos = 0
    for gene_id in failed_afas:
        canonical_human_exons = filter (lambda x:  x.is_canonical and x.is_coding, gene2exon_list(cursor, gene_id))
        if not canonical_human_exons: 
            no_exons += 1
            continue
        # the exons are not guaranteed to be in order
        canonical_human_exons.sort(key=lambda exon: exon.start_in_gene)
        # reconstruct  per-exon alignments with orthologues
        [alnmt_pep, alnmt_dna] = make_exon_alignments(cursor, ensembl_db_name, canonical_human_exons,
                                                      mitochondrial, min_similarity, flank_length)

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # we want to be able to retrieve the info starting from whichever end, so we construct the following maps:
        # to find all exons from an ortohologue, that map to a given human exon:
        #      human_exon_to_ortho_exon:  human_exon_to_ortho_exon[concat_seq_name][human_exon].append(exon_seq_name)
        # given a sequence name, retrieve all exons that belong to it
        #     sequence_to_exon: sequence_name_to_exon_names[concat_seq_name].append(exon_seq_name)
        # in the other direction - to find all human exons that a given exon  from orthologous sequence maps to
        #      ortho_exon_to_human_exon:  ortho_exon_to_human_exon[exon_seq_name].append(human_exon)
        # finally, collect in one place info about maps between human and orthologue that are not one-to-one
        #      overlapping_maps: overlapping_maps[concat_seq_name].append([human_exons, ortho_exons])
        [human_exon_to_ortho_exon, sequence_name_to_exon_names, 
         ortho_exon_to_human_exon, overlapping_maps] = make_atlas(cursor, ensembl_db_name, canonical_human_exons, 
                                                                  alnmt_pep, trivial_name)
        # the alignment always has human sequence, but if it is the only one
        # (see for example RPL41, ENSG00000229117, a 25 residue peptide,  for which NCBI REfseq
        # reports a single  confirmed  homologue in mouse, but Ensembl reports no orthologues at all)
        if ( len(sequence_name_to_exon_names) <= 1):
            no_orthos += 1
            continue
         
    
    print
    print "failure cases"
    print "\t no exons", no_exons
    print "\t no orthologues ", no_orthos


    cursor.close()
    db.close()



#########################################
if __name__ == '__main__':
    main()
