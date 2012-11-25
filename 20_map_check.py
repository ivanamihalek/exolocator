#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  get_maps, Map
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader


from bitstring import Bits


#########################################
def main():

    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    # find db ids adn common names for each species db
    [all_species, ensembl_db_name] = get_species (cursor)
    species                        = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_list                      = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
    with_map = 0
    tot = 0
    #or gene_id in [412667]: #  wls
    for gene_id in gene_list: 
        
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        #print  gene2stable(cursor, gene_id), get_description (cursor, gene_id)

        # find all exons we are tracking in the database
        human_exons = gene2exon_list(cursor, gene_id)
        human_exons.sort(key=lambda exon: exon.start_in_gene)
        has_a_map = False
        for human_exon in human_exons:
            if ( not human_exon.is_canonical or  not human_exon.is_coding): continue
            #print  
            #print "\t human",   human_exon.exon_id,  human_exon.is_known
            #print "\t", get_exon_pepseq(cursor, human_exon.exon_id, human_exon.is_known, 
            #                                         ensembl_db_name['homo_sapiens'])
            #print "checking maps ..."
            maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
            if maps:
                has_a_map = True
            #else:
            #print"no maps"
            if 0:
                for map in maps:
                    species = map.species_2
                    unaligned_sequence = get_exon_pepseq(cursor, map.exon_id_2, map.exon_known_2, 
                                                         ensembl_db_name[map.species_2])
                    print "\t", species, map.exon_id_2, map.exon_known_2
                    print "\t", unaligned_sequence
        if has_a_map: with_map +=1
        tot += 1
        print gene_id, tot, with_map

    cursor.close()
    db.close()

    print tot, with_map
            

#########################################
if __name__ == '__main__':
    main()
