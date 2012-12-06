#!/usr/bin/python

import MySQLdb, commands, re
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  Map,  get_maps
from   el_utils.tree    import  species_sort
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   bitstring        import Bits

#########################################
def make_exon_alignment(cursor, ensembl_db_name, human_exon):

    sequence = {}

    # find all other exons that map to the human exon
    maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)

    for map in maps:
        # get the raw (unaligned) sequence for the exon that maps onto human
        exon_seqs = get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
        if (not exon_seqs ):
            print map
            exit (1)
        exon_seqs = exon_seqs[1:] # the first entry is database id
        [pep_seq, left_flank, right_flank, dna_seq] = exon_seqs

        # inflate the compressed sequence
        unaligned_sequence = pep_seq
        if not unaligned_sequence:continue

        bs = Bits(bytes=map.bitmap)
        # check bitmap has correct number of 1s
        if ( not bs.count(1) == len(unaligned_sequence)):
            print "bitmap check fails (?)"
            continue

        # rebuild aligned sequence
        usi = iter(unaligned_sequence)
        reconstructed_sequence = "".join(('-' if c=='0' else next(usi) for c in bs.bin))

        # come up with a unique name for this sequence
        species       = map.species_2                
        sequence_name = species + "_" + str(map.exon_id_2)

        sequence[sequence_name] = reconstructed_sequence
    
    return sequence


#########################################
def print_alignment (sorted_species, alignment):
            
    for species in sorted_species:
        for seq_name, sequence in alignment.iteritems():
            if (species in seq_name):
                print ">"+seq_name
                print sequence
                cognate_found = True
            

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

    # note this function: it will return the list of species in Ensembl
    # according to their taxonomical distance from the last argument in the function
    sorted_species = species_sort(cursor, all_species, 'homo_sapiens')
    
    species  = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    # for each human gene
    #for gene_id in gene_ids:
    for gene_id in [412667]: #  wls
        print gene_id, gene2stable(cursor, gene_id), get_description (cursor, gene_id)
        alignment = {}
        # find all human exons we are tracking in the database
        human_exons = gene2exon_list(cursor, gene_id)
        human_exons = filter (lambda exon: exon.covering_exon < 0 and exon.is_canonical, human_exons)
        human_exons.sort(key=lambda exon: exon.start_in_gene)
        
        # note: the map [between the human and other species' exons]
        # is already stored in the database - check exon_map table
        # the function make_exon_alignment() is using precisely that info to
        # re-assemble the alignment
        #
        # and the other way round - the results should be stored in the exon_map table, 
        # wiht the source annotated as "sw"
        for human_exon in human_exons:
            alignment[human_exon] =  make_exon_alignment(cursor, ensembl_db_name, human_exon)
 
        # now fix the exon, and check in species direction - which species are missing the cognate exon?
        print " %30s  %3d " % (species, len(human_exons))
        for species in sorted_species[1:]: # the  first is human
            exon_ct = 0
            for human_exon in human_exons:
                for seq_name in alignment[human_exon].keys():
                    if species in seq_name:
                        exon_ct += 1
                        break
            print " %30s  %3d " % (species, exon_ct)
            

        exit(1)


#########################################
if __name__ == '__main__':
    main()


