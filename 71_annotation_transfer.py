#!/usr/bin/python

import MySQLdb, subprocess, re, commands
import copy, pdb

from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.utils   import  *
from el_utils.map     import  Map, get_maps
from el_utils.tree    import  species_sort

from el_utils.special_gene_sets  import get_theme_ids, get_complement_ids
from el_utils.almt_cmd_generator import AlignmentCommandGenerator
from el_utils.config_reader      import ConfigurationReader
from el_utils.threads import  *
from subprocess import Popen, PIPE, STDOUT

#########################################
verbose = True
ignorance_indicators = ['novel', 'uncharacterized']

########################################
def find_orthologues(cursor, ensembl_db_name, gene_id, query_species, species):

    ortho_gene_ids = {}

    switch_to_db (cursor, ensembl_db_name[species])
    stable_id = gene2stable(cursor, gene_id=gene_id)
    switch_to_db (cursor, ensembl_db_name['compara'])
    member_id = stable2member(cursor, stable_id)
    
    
    for ortho_type in ['ortholog_one2one','ortholog_one2many','ortholog_many2many']:
        ortho_gene_ids[ortho_type]  = get_orthologues_from_species(cursor, ortho_type, member_id, species)
        
    return ortho_gene_ids
    
########################################
def find_annotation (cursor, ensembl_db_name, species_list, gene_id):

    source_species = species_list[0]
    annotation     = 'none'
    query_species  = species_list[0]

    for species in species_list:

        if species==query_species:
            orthologous_gene_ids = {'self':[gene_id]}
            continue
        else:
            orthologous_gene_ids = find_orthologues(cursor, ensembl_db_name, gene_id, query_species, species)
            if not orthologous_gene_ids: continue

        switch_to_db (cursor, ensembl_db_name[species])
        for ortho_type, gene_ids  in orthologous_gene_ids.iteritems():
            for orthologous_gene_id in gene_ids:
                # does the orthologue have the description?
                description = get_description(cursor, orthologous_gene_id)
                if description and not filter (lambda x: x in description.lower(), ignorance_indicators):
                    source_species = species
                    annotation     = description  
                    if not species == 'oryctolagus_cuniculus':
                        print 'annotation found in', species
                        print 'orthology type', ortho_type
                        print annotation
        if not annotation == 'none': 
            #break
            exit(1)
    return [source_species, annotation]
  
#########################################
def annotate(gene_list, db_info):
    # 
    [local_db, all_species, ensembl_db_name, species] = db_info
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
        acg = AlignmentCommandGenerator()
    else:
        db  = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    if verbose: print "thread %s annotating %s "% (get_thread_name(), species)

    if not species == 'oryctolagus_cuniculus':
        print 'The preferred list of species is hardcoded for the rabbit. Consider modifying.'
        exit(1)
    preferred_species = [species, 'mus_musculus', 'rattus_norvegicus', 'homo_sapiens']
    nearest_species_list = species_sort(cursor, all_species, species)
    species_list = preferred_species + filter (lambda x: x not in preferred_species, nearest_species_list)

    for gene_id in gene_list:
        switch_to_db (cursor, ensembl_db_name[species])
        # Get stable id and description of this gene
        stable      = gene2stable    (cursor, gene_id)

        ####################
        [annot_source, annotation] = find_annotation (cursor, ensembl_db_name, species_list, gene_id)
        
        # find splices
        # output

    cursor.close()
    db.close()

#########################################
#########################################
def main():
    
    species    = 'oryctolagus_cuniculus'
    no_threads = 1
 

    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <species> <number of threads>" % sys.argv[0]
        exit(1)
    elif len(sys.argv)==3:

        species = sys.argv[1].lower()

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
    ensembl_db_name['compara']     = get_compara_name(cursor)

    print '======================================='
    print "running %s for %s " % (sys.argv[0], species)

    switch_to_db (cursor,  ensembl_db_name[species])
    gene_list = get_gene_ids (cursor, biotype='protein_coding')

    cursor.close()
    db.close()

    parallelize (no_threads, annotate, gene_list, [local_db, all_species, ensembl_db_name, species])
    
    return True

#########################################
if __name__ == '__main__':
	main()

