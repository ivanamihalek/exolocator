#!/usr/bin/python -u
# make the best alignment we can using the maps
# we currently have at hand

import pdb
#pdb.set_trace()

import MySQLdb, commands, re, os
import inspect
from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.utils   import  *
from el_utils.el_specific   import  *
from el_utils.map     import  Map, get_maps, map2exon
from el_utils.tree    import  species_sort
from el_utils.ncbi    import  taxid2trivial
from el_utils.config_reader      import ConfigurationReader
from bitstring import Bits
from alignment import * # C implementation of smith waterman
from   random  import choice
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

import pdb

verbose = True

#########################################
def find_human_orthologue(cursor,  ensembl_db_name, species, gene_id):

    genome_db_id  = species2genome_db_id (cursor, species)

    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    qry = "select gene_id from orthologue "
    qry += " where cognate_gene_id=%d " % gene_id
    qry += " and cognate_genome_db_id=%d " % genome_db_id

    rows = search_db (cursor, qry)
    
    if not rows:
        search_db (cursor, qry, verbose=True)
        exit(1)

    try:
        retval = int(rows[0][0])
    except:
        retval = None

    return retval


#########################################
def find_human_exon_map (cursor, ensembl_db_name, human_gene_id, species, exon_id, exon_known):

    genome_db_id  = species2genome_db_id (cursor, species)
    
    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    qry = "select exon_id, exon_known from exon_map " 
    qry += " where  cognate_genome_db_id = %d " % int(genome_db_id)
    qry += " and  cognate_exon_id    = %d "     % int(exon_id)
    qry += " and  cognate_exon_known = %d "     % int(exon_known)
    
    rows = search_db (cursor, qry)
    
    if not rows:
        search_db (cursor, qry, verbose=True)
        exit(1)

    try:
        retval = [int(rows[0][0]), int(rows[0][1]) ]
    except:
        retval = None

    return retval


#########################################
def  reconstruct_alignment (cursor,  cfg, ensembl_db_name, species, exon_id, exon_known, sorted_species, output_fnm_root):


    switch_to_db (cursor, ensembl_db_name[species])
    # which gene does it belong to
    gene_id = exon_id2gene_id (cursor, ensembl_db_name[species],  exon_id, exon_known)
    human_gene_id = find_human_orthologue(cursor,  ensembl_db_name, species, gene_id)

    # is it mitochondrial
    mitochondrial = is_mitochondrial(cursor, human_gene_id)

    # find the human_exon that this exon maps to
    [human_exon_id, human_exon_known] = find_human_exon_map (cursor, ensembl_db_name, human_gene_id, 
                                                             species, exon_id, exon_known)

    min_similarity = cfg.get_value('min_accptbl_exon_sim') # minimum acceptable similarity between mapped exons
    # reconstruct the alignment
    is_known       =  1
    flank_length   = 20
    [alnmt_pep, alnmt_dna] = \
        make_exon_alignment(cursor, ensembl_db_name, human_exon_id, is_known,  
                            mitochondrial, min_similarity, flank_length)   
    
    # sort the names according to the distance from the query species
    sorted_names = sort_names (sorted_species, alnmt_pep)
    # spit out
    output_fasta( output_fnm_root+".pep.afa", sorted_names, alnmt_pep)
    output_fasta( output_fnm_root+".dna.afa", sorted_names, alnmt_dna)

    return True


#########################################
def main():

    if len(sys.argv) < 5:
        print "Usage: %s <species>  <exon_id> <exon_known> <output_name_root>" % sys.argv[0]
        exit(1)

    species         = sys.argv[1]
    exon_id         = int(sys.argv[2])
    exon_known      = int(sys.argv[3])
    output_fnm_root = sys.argv[4]

    local_db = False
    
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    sorted_species  = species_sort(cursor, all_species, species)

    reconstruct_alignment (cursor,  cfg, ensembl_db_name, species, exon_id, exon_known, sorted_species, output_fnm_root) 
    cursor.close()
    db.close()

    
    return True


#########################################
if __name__ == '__main__':
    main()
