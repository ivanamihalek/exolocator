#!/usr/bin/python

import MySQLdb, subprocess, re, commands
import copy, pdb

from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.el_specific import  *
from el_utils.utils   import  *
from el_utils.map     import  Map, get_maps
from el_utils.tree    import  species_sort

from el_utils.special_gene_sets  import get_theme_ids, get_complement_ids
from el_utils.almt_cmd_generator import AlignmentCommandGenerator
from el_utils.config_reader      import ConfigurationReader
from el_utils.threads import  *
from subprocess import Popen, PIPE, STDOUT

#########################################
verbose = False
ignorance_indicators = ['novel', 'uncharacterized', 'cDNA']

########################################
def find_orthologues(cursor, ensembl_db_name, gene_id, query_species, species):

    ortho_gene_ids = {}

    switch_to_db (cursor, ensembl_db_name[query_species])
    query_stable_id = gene2stable(cursor, gene_id=gene_id)
    switch_to_db (cursor, ensembl_db_name['compara'])
    query_member_id = stable2member(cursor, query_stable_id)
    
    for ortho_type in ['ortholog_one2one','ortholog_one2many','ortholog_many2many','possible_ortholog', 'apparent_ortholog_one2one']:
        ortho_gene_ids[ortho_type]  = get_orthologues_from_species(cursor, ensembl_db_name, 
                                                                   ortho_type, query_member_id, species)
        if ortho_gene_ids[ortho_type]: break # if we got one2one, we're happy an we move on
    return ortho_gene_ids
    
########################################
def find_annotation (cursor, ensembl_db_name, species_list, gene_id):

    query_species  = species_list[0]

    for species in species_list:

        source_species = 'none'
        source_ids     = 'none'
        annotation     = 'none'
        ortho_type     = 'none'

        if species==query_species:
            # if the original annotation exists, we're done
            # (in the function call, the query species is the first)
            orthologous_gene_ids = {'self':[gene_id]}
            ortho_type = 'self'
        else:
            # find orthologues in this species (unfortunately, there might be more then one)
            orthologous_gene_ids = find_orthologues(cursor, ensembl_db_name, gene_id, query_species, species)
            if not orthologous_gene_ids: continue


        switch_to_db (cursor, ensembl_db_name[species])
        for ortho_type, gene_ids in orthologous_gene_ids.iteritems():
            # can I have several orthology types for the same species? 
            # it shouldn't be so ...
            if not gene_ids: continue
            source_species  = species
            description = ''
            source_ids  = ''
            #print species, ortho_type, gene_ids
            for orthologous_gene_id in gene_ids:
                # does the orthologue have the description?
                this_gene_description = get_description(cursor, orthologous_gene_id)
                if this_gene_description and not filter (lambda x: x in this_gene_description.lower(), ignorance_indicators):
                    if description: description += '; '
                    description += this_gene_description  
                    if source_ids: source_ids += ';'
                    source_ids  += gene2stable(cursor,orthologous_gene_id) 
            if description: 
                annotation = description
                break
        if not annotation=='none': break

    return [source_species, ortho_type, annotation, source_ids]
  

#########################################
def annotate(gene_list, db_info):
    # 
    [local_db, all_species, ensembl_db_name, species] = db_info
    db  = connect_to_mysql()
    cfg = ConfigurationReader()
    acg = AlignmentCommandGenerator()
    cursor = db.cursor()

    if verbose: print "thread %s annotating %s "% (get_thread_name(), species)

    if not species == 'oryctolagus_cuniculus':
        print 'The preferred list of species is hardcoded for the rabbit. Consider modifying.'
        exit(1)

    preferred_species = [species, 'mus_musculus', 'rattus_norvegicus', 'homo_sapiens']
    nearest_species_list = species_sort(cursor, all_species, species)
    species_list = preferred_species + filter (lambda x: x not in preferred_species, nearest_species_list)

    inf = erropen("temp_out.fasta", "w")

    for gene_id in gene_list:
    #for gene_id in [90020]:
        switch_to_db (cursor, ensembl_db_name[species])
        ####################
        # get stable id and description of this gene
        stable_id      = gene2stable    (cursor, gene_id)
        if not gene_list.index(gene_id)%100: print gene_list.index(gene_id), "out of", len(gene_list)
        if verbose: print "============================================="
        if verbose: print gene_id, stable_id
        ####################
        # find the annotation from the preferred source organism
        [annot_source, orthology_type, annotation, ortho_stable_ids] = find_annotation (cursor, ensembl_db_name, species_list, gene_id)
        if verbose: print annot_source, "**", orthology_type, '**', annotation

        ###################
        # find splices (for now find the canonical splice)
        switch_to_db (cursor, ensembl_db_name[species])
        canonical_splice = get_canonical_transl (acg, cursor, gene_id, species)

        # output
        if orthology_type == 'self' or annotation== 'none':
            header = ">{0} {1}".format(stable_id, annotation)
        else:
            header = ">{0} {1} [by sim to {2}, {3}]".format(stable_id, annotation, annot_source, ortho_stable_ids)

        print >>inf, header
        print >>inf, canonical_splice
 
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

    db  = connect_to_mysql()
    cfg = ConfigurationReader()
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    ensembl_db_name['compara']     = get_compara_name(cursor)

    print
    print "running %s for %s " % (sys.argv[0], species)
    print 

    switch_to_db (cursor,  ensembl_db_name[species])
    gene_list = get_gene_ids (cursor, biotype='protein_coding')

    cursor.close()
    db.close()

    parallelize (no_threads, annotate, gene_list, [local_db, all_species, ensembl_db_name, species])
    
    return True

#########################################
if __name__ == '__main__':
	main()

