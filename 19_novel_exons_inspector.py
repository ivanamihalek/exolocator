#!/usr/bin/python

import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl import *
from el_utils.special_gene_sets  import get_theme_ids, get_complement_ids
from el_utils.config_reader      import ConfigurationReader
from el_utils.map       import *
from el_utils.utils       import *

#########################################
def  map_cleanup (cursor, ensembl_db_name, human_exons):
    
    switch_to_db(cursor,ensembl_db_name['homo_sapiens']) 
    for exon in human_exons:
        qry  = "delete from exon_map where exon_id = %d " % exon.exon_id
        qry += " and exon_known = %d " % exon.is_known
        qry += " and cognate_exon_known > 1 " 
        qry += " and similarity is NULL" 
        rows = search_db (cursor, qry, verbose=False)


    return True


#########################################
def main():
    
    special    = 'one'
    no_threads = 1
    local_db = False
    

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    print "using", special, "set"
    gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special)

    # loop over all genes
    for human_gene_id in gene_list:
        
 	human_stable      = gene2stable    (cursor, human_gene_id)
        human_description = get_description(cursor, human_gene_id)
	print human_gene_id, human_stable, human_description
   
  	human_exons = [e for e in gene2exon_list(cursor, human_gene_id) 
                       if e.covering_exon < 0 and e.is_canonical and e.is_known]
        if not human_exons: 
            print "\t\t no exons found"
            continue

        map_cleanup (cursor, ensembl_db_name, human_exons)

	human_exons.sort(key=lambda exon: exon.start_in_gene)
        # loop over all exons in this gene
        maps_for_exon = {}
        for he in human_exons:
            he.stable_id = exon2stable (cursor, he.exon_id, ensembl_db_name['homo_sapiens'])
            he.pepseq = get_exon_pepseq (cursor, he,  ensembl_db_name['homo_sapiens'])
            # maps cleanup: get rid of maps that have "none" as similarity

            maps_for_exon[he] =  get_maps(cursor, ensembl_db_name, he.exon_id, he.is_known) # exon data
            if not maps_for_exon[he]: continue

            maps_for_exon[he] = filter (lambda m: m.source == 'sw_sharp' or m.source == 'usearch', 
                                        maps_for_exon[he])

            if not maps_for_exon[he]: continue

            print
            print "======================================"
            print he.exon_id, he.stable_id 

            for m in maps_for_exon[he]:

                ret = get_exon_seqs (cursor, m.exon_id_2, m.exon_known_2, ensembl_db_name[m.species_2])
                print 
                print m.source, m.similarity, m.species_2, m.exon_id_1, m.exon_id_2
                
                if not ret:  
                    print "exon seqs not found (?)"
                else:
                    [exon_seq_id, protein_seq, pepseq_transl_start, 
                     pepseq_transl_end, left_flank, right_flank, dna_seq] = ret
                    [seq1, seq2] = unfold_cigar_line (he.pepseq, protein_seq, m.cigar_line)
                    print seq1
                    print seq2

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
