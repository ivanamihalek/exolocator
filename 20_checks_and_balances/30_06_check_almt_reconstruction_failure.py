#!/usr/bin/python -u
# make the best alignment we can using the maps
# we currently have at hand

import pdb
#pdb.set_trace()

import MySQLdb, commands, re, os, time
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
    for gene_id in gene_list:
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)
        if  check_afa_age (cfg, stable_id, max_days=30) == "new": 
            new_afas += 1
            continue                               
        elif  check_afa_age (cfg, stable_id, max_days=300) == "new": 
            old_afas += 1
            continue                               
    print "total genes", len(gene_list)
    print "old  afas", old_afas
    print "new  afas", new_afas

    cursor.close()
    db.close()
