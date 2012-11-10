#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  get_species, get_gene_ids, gene2stable, stable2gene
from   el_utils.ensembl import  get_compara_name, get_description, gene2exon_list
from   el_utils.ensembl import  get_exon_pepseq, genome_db_id2species
from   el_utils.utils   import  erropen
from   el_utils.map     import  Map
#########################################
def main():

      
    db     = connect_to_mysql()
    cursor = db.cursor()
    acg    = AlignmentCommandGenerator()

    [all_species, ensembl_db_name] = get_species (cursor)
    # find db ids adn common names for each species db


    # build a tree, and sort all species according to how they are from human
    # sorted_species = species_sort(all_species)

    species='homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)


    for gene_id in gene_ids:
        print gene_id, gene2stable(cursor, gene_id), get_description (cursor, gene_id)
        human_exons = gene2exon_list(cursor, gene_id)
        # for each species from my ordered list 
        #     for each human canonical exon
        #         see if there is a map on this species
        #         if yes, append the exon to this species' list of exons
        #         if not, append an empty string
        #     if at least one exon nonempty, append this species to the list of usable species
        # for each exon 
        #     make fasta with exons from all species
        #     align using mafft
        #     read in the alignment
        # for each species
        #     stitch all the pieces from all aliggnments 
        #     decorate with Z's to indicate the edges of each exon
        #     if the  species was not in this alignment, add  a stretch of gaps
        # output to a file, 'best_afas/stable_id.afa
        # output notes file
    pass

#########################################
if __name__ == '__main__':
    main()
