#!/usr/bin/python

import MySQLdb, commands, re
from   el_utils.mysql   import  connect_to_mysql, switch_to_db
from   el_utils.mysql   import  search_db, store_or_update
from   el_utils.ensembl import  get_species, get_gene_ids, gene2stable, stable2gene
from   el_utils.ensembl import  get_compara_name, get_description, gene2exon_list
from   el_utils.ensembl import  get_exon_pepseq, genome_db_id2species, species2taxid
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  Map,  get_maps
from   el_utils.tree    import  species_sort
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader

#########################################
def main():

    db     = connect_to_mysql()
    cursor = db.cursor()
    cfg    = ConfigurationReader()
    acg    = AlignmentCommandGenerator()
    # find db ids adn common names for each species db
    [all_species, ensembl_db_name] = get_species (cursor)
    

    species  = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    # for each human gene
    for gene_id in gene_ids:
        print gene_id, gene2stable(cursor, gene_id), get_description (cursor, gene_id)
        per_species_map = {}
        # find all exons we are tracking in the database
        human_exons = gene2exon_list(cursor, gene_id)
        human_exons.sort(key=lambda exon: exon.start_in_gene)
        for human_exon in human_exons:
            # find all orthologous exons the human exon  maps to
            maps  = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
            for map in maps:
                if ( not per_species_map.has_key(map.species_2) ): 
                     per_species_map[map.species_2] = {}
                per_species_map[map.species_2][human_exon] = map

        for species in all_species:
            print "\t ", species
            if (  not per_species_map.has_key(species) ): 
                print "\t      no gene found "
                continue
            human_exon_ct = 0
            #if 0:
            for human_exon in human_exons:
                human_exon_ct += 1
                if per_species_map[species].has_key(human_exon):
                    map =  per_species_map[species][human_exon]
                    print "\t      %3d  %8.2f" % ( human_exon_ct, map.similarity)
                else:
                    print "\t      %3d  none" % ( human_exon_ct)

        exit (1)




#########################################
if __name__ == '__main__':
    main()
