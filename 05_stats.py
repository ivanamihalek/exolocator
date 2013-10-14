#!/usr/bin/python

import MySQLdb
from random import choice
from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.map     import  Map, get_maps
from el_utils.special_gene_sets  import  get_theme_ids
from el_utils.config_reader      import ConfigurationReader


#########################################
def ortho_stats (cursor, ensembl_db_name, species, gene_ids):
    
    if species == 'homo_sapiens':
        tables = ['orthologue', 'unresolved_ortho', 'paralogue']
    else:
        tables = ['paralogue']

    # agjafgj 
    for ortho_table in tables:
        # how many orthologue pairs, in principle
        # orthologue tables are in homo_sapiens database
        # each species has its own paralogue table
        qry = "select count(1) from %s.%s" % (ensembl_db_name[species], ortho_table)
        rows = search_db (cursor, qry, verbose=True)
        if (not rows):
            rows = search_db (cursor, qry, verbose=True)
        print ortho_table, " table size: ", rows[0][0]

    # how many per gene
    if False:
        for ortho_table in tables:
        
            print "histogram for ", ortho_table

            histogram = {}
            for gene_id in gene_ids:

                qry  = "select count(1) from %s.%s" % (ensembl_db_name[species], ortho_table)
                qry += " where gene_id = %d " % gene_id
                rows = search_db (cursor, qry)
                if (not rows):
                    rows = search_db (cursor, qry, verbose=True)
                if (not histogram.has_key(rows[0][0])):
                    histogram[rows[0][0]] = 0
                histogram[rows[0][0]] += 1
 
            for number_of_orthos in sorted (histogram.keys()):
                 print " %4d  %4d " % (number_of_orthos, histogram[number_of_orthos])


#########################################
def make_map_table (cursor, ensembl_db_name, all_species, human_exons):
    
    # make 'table' of maps, which is either pointer to the map if it exists, or None
    map_table  = {}
    for species in all_species:
        map_table[species] = {}
        for he in human_exons:
            map_table[species][he] = None

    maps_for_exon = {}
    for he in human_exons:
        maps_for_exon[he] =  get_maps(cursor, ensembl_db_name, he.exon_id, he.is_known) # exon data
        for m in maps_for_exon[he]:
            #if not m.source == 'ensembl': continue
            #if m.similarity < 0.33333: continue
            if not m.species_2 in all_species: continue
            map_table[m.species_2][he] = m
           #if m.source =='sw_sharp': print m.source
    # get rid of species that do not have the gene at all
    for species in all_species:
        one_exon_found = False
        for he in human_exons:
            if map_table[species][he]:
                one_exon_found = True
                break
        if not one_exon_found:
            del map_table[species]

    return map_table


#########################################
def exon_stats (cursor, ensembl_db_name, all_species, human_gene_list):

    tot =             0
    exons_not_found = 0
    map_fail =  0

    table_size = 0
    holes      = 0
    from_sw_sharp = 0
    from_usearch  = 0

    for human_gene_id in human_gene_list:
    #for ct in range (500):
        #human_gene_id =  choice(human_gene_list)
        tot += 1

	human_exons = [e for e in gene2exon_list(cursor, human_gene_id, db_name=ensembl_db_name['homo_sapiens']) 
                       if e.covering_exon < 0 and e.is_canonical and e.is_known]
        
        if not human_exons: 
            exons_not_found += 1
            continue

        map_table = make_map_table (cursor, ensembl_db_name, all_species, human_exons)
        if not map_table:
            map_fail += 1

        for species in map_table.keys():
            for he in human_exons:
                table_size += 1
                map = map_table[species][he]
                if not map or map.similarity < 0.2:
                    holes += 1
                elif map.source =='sw_sharp':
                    holes += 1
                    from_sw_sharp += 1
                elif map.source =='usearch':
                    holes += 1
                    from_usearch  += 1
        if table_size and holes: print " % 4d  %4d  %5.2f   %5.2e   %5.2e " %  \
                (table_size, holes,  float(holes)/table_size, float(from_usearch)/holes, float(from_sw_sharp)/holes)

 
        print " %20s   %5d " % ("tot", tot)

    print " %20s   %5d " % ("exons_not_found", exons_not_found)
    print " %20s   %5d " % ("map_fail", map_fail)
    if table_size and holes: print " % 4d  %4d  %5.2f   %5.2e   %5.2e " %  \
            (table_size, holes,  float(holes)/table_size, float(from_usearch)/holes, float(from_sw_sharp)/holes)


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

    mammals = ['ailuropoda_melanoleuca',   'bos_taurus',  'callithrix_jacchus',  'canis_familiaris',  
               'cavia_porcellus',  'choloepus_hoffmanni',  'dasypus_novemcinctus',  'dipodomys_ordii',  
               'echinops_telfairi',  'equus_caballus',  'erinaceus_europaeus',  'felis_catus',   'gorilla_gorilla',  
               'ictidomys_tridecemlineatus',   'loxodonta_africana',  'macaca_mulatta',  'macropus_eugenii',    
               'microcebus_murinus',  'monodelphis_domestica',  'mus_musculus',  'mustela_putorius_furo',  
               'myotis_lucifugus',  'nomascus_leucogenys',  'ochotona_princeps',   'ornithorhynchus_anatinus',  
               'oryctolagus_cuniculus', 'otolemur_garnettii', 'pan_troglodytes', 'pongo_abelii',  
               'procavia_capensis', 'pteropus_vampyrus', 'rattus_norvegicus', 'sarcophilus_harrisii',  
               'sorex_araneus', 'sus_scrofa', 'tarsius_syrichta',  'tupaia_belangeri',  'tursiops_truncatus',  
               'vicugna_pacos']

    for species in all_species:
        switch_to_db (cursor,  ensembl_db_name[species])
        print
        print "** species: ", species
        known_genes = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        print "\t known genes: ", len(known_genes)
        predicted_genes = get_gene_ids (cursor, biotype='protein_coding', is_known=0)
        print "\t predicted genes: ", len(predicted_genes)
        # sanity check:
        all_genes = get_gene_ids (cursor, biotype='protein_coding')
        print "\t predicted genes: ", len(all_genes)

    if 0:
        for special in  ['wnt_pathway']:

            species   = 'homo_sapiens'
            switch_to_db (cursor,  ensembl_db_name[species])

            # how many genes are we talking about, in the first place?
            if special:
                print "using", special, "set"
                human_gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special )
            else:
                print "using all protein coding genes"
                switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
                human_gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
            print
            print "protein coding genes in human:  %10d " %  len(human_gene_list),
            print " (known, not predicited)"

            # how many  have orthologues reported?
            for species in all_species:
                ortho_stats (cursor,  ensembl_db_name, species, human_gene_list)    

            # how many  have exons reported? (this becomes meaningful only later in the pypeline)
            #exon_stats (cursor, ensembl_db_name, mammals, human_gene_list)    


    # how often does it happen that one  exon does not have
    # a map while the others do

    # how many of those can actually be found, and how many are gaps in the seqeunces
    # (are the gaps in the sequence commesurat withe the coverage?)
    
    # ow many are at the scaffold boundary?

    # how many of those cane be patched by brute force?
    # which sequences are more patchable?

    # what are the lessons learned (1) about biology, (2) about tchnology/tools that 
    # we need to find the missing exons?

    #
       

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

'''
       'telomere_maintenance',
                     'nonhom_end_joining', 'egfr_signaling', 'cell_cycle_checkpoints',
                     'genecards_top500', 'wnt_pathway',  'enzymes'

       prev_end = 0
        for human_exon in human_exons:

            first_exon = (human_exons.index(human_exon) == 0)
            if  first_exon: origin = human_exon.start_in_gene
            print " &&&   %10d    %10d   %10d" % (human_exon.start_in_gene-origin, 
                                                     human_exon.end_in_gene-human_exon.start_in_gene+1,
                                                     human_exon.start_in_gene - prev_end)
            prev_end = human_exon.end_in_gene

'''
