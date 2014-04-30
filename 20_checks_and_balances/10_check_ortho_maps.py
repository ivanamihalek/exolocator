#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import sys, MySQLdb, commands, re
from bitstring import Bits
from   random           import  choice
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  Map, get_maps, map2exon 
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.special_gene_sets  import  get_theme_ids


verbose = True

#########################################
def main():


    no_threads = 1
    special    = None

    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <set name> <number of threads> " % sys.argv[0]
        exit(1)
    elif len(sys.argv)==3:

        special = sys.argv[1]
        special = special.lower()
        if special == 'none': special = None

        no_threads = int(sys.argv[2])

    local_db = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    # find db ids adn common names for each species db
    [all_species, ensembl_db_name] = get_species (cursor)
    species                        = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])

    if special:
        print "using", special, "set"
        gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special )
    else:
        print "using all protein coding genes"
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        
    incomplete = 0
    genes_checked = 0
    for gene_id in gene_list: 
    #for sampling_count in range(1000):
 
        #gene_id = choice(gene_list)
        genes_checked += 1
        with_map = 0
        tot      = 0
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        print  gene2stable(cursor, gene_id), get_description (cursor, gene_id)

        # find all exons we are tracking in the database
        human_exons = gene2exon_list(cursor, gene_id)
        human_exons.sort(key=lambda exon: exon.start_in_gene)
        has_a_map = False
        for human_exon in human_exons:
            if not human_exon.exon_id == 8339254: continue
            if (not human_exon.is_canonical or  not human_exon.is_coding): continue
            if verbose:
                print  
                print "\t human",   human_exon.exon_id,  human_exon.is_known
                print "\t ", get_exon_pepseq(cursor, human_exon, ensembl_db_name['homo_sapiens'])
                print "\t checking maps ..."
            maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
            tot += 1
            if maps:
                has_a_map = True
                with_map += 1
                #print "ok"
            else:
                print"no maps"
                print human_exon
                exit(1)
                pass
            if verbose:
                for map in maps:
                    species            = map.species_2
                    #if not species == 'procavia_capensis': continue
                    exon               = map2exon(cursor, ensembl_db_name, map)
                    unaligned_sequence = get_exon_pepseq(cursor, exon, ensembl_db_name[species])
                    if ( map.similarity):
                        print "\t", species,  map.source, map.exon_id_2, map.exon_known_2
                        print "\tmaps to ",  map.exon_id_1, map.exon_known_1
                        print "\tsim",  map.similarity,
                        print "\tsource",  map.source
                        print "\t", unaligned_sequence
                        if not map.bitmap:
                            print "\t bitmap not assigned"
                        else:
                            bs = Bits(bytes=map.bitmap)
                            reconst_pepseq = ''
                            if (not bs.count(1) == len(unaligned_sequence)): 
                                print "\talnd seq mismatch"
                            
                            else:
                                usi = iter(unaligned_sequence)
                                for c in bs.bin:
                                    if c == '0': reconst_pepseq += '-'
                                    else:        reconst_pepseq += next(usi)
                                print "\tbinary   : ", bs.bin
                                print "\talnd seq: ", reconst_pepseq
                        print
        if not tot== with_map:
            print "####  gene id: %d   total exons: %d     with map:  %d   ( = %d%%) " % \
                (gene_id,  tot,  with_map, int(float(with_map)/tot*100) )
            incomplete += 1

    print "genes checked: %d,  incomplete: %d"  %  (genes_checked, incomplete)
    cursor.close()
    db.close()

    print tot, with_map
            

#########################################
if __name__ == '__main__':
    main()
