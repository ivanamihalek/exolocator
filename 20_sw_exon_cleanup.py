#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl import  *
from el_utils.config_reader      import ConfigurationReader
from el_utils.special_gene_sets  import human_genes_w_sw_sharp_annotation, get_theme_ids
# BioPython
from Bio.Seq      import Seq

#########################################
def get_ok_human_exons (cursor, ensembl_db_name, gene_id):
    canonical_human_exons = []
    # find all exons we are tracking in the database
    human_exons     = gene2exon_list(cursor, gene_id)
    for human_exon in human_exons:
        if not human_exon.is_canonical or  not human_exon.is_coding or  not human_exon.is_known:
            continue
        
        if not get_exon_seqs (cursor, human_exon.exon_id, 1, ensembl_db_name['homo_sapiens']):
            continue
        canonical_human_exons.append(human_exon)
    return canonical_human_exons

#########################################
def main():
    
    no_threads = 1
    local_db   = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [species_list, ensembl_db_name] = get_species   (cursor)
    gene_list                       = get_theme_ids (cursor, ensembl_db_name, cfg, 'genecards_top500')

    for gene_id in gene_list:

        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)
        print "#############################################"
        print gene_id, stable_id, get_description (cursor, gene_id)

        human_exons = get_ok_human_exons (cursor,ensembl_db_name,  gene_id)
        tot    = 0
        tot_ok = 0

        for human_exon in human_exons:
            [exon_seq_id, human_protein_seq, pepseq_transl_start, pepseq_transl_end, 
             left_flank, right_flank, dna_seq] = get_exon_seqs (cursor, human_exon.exon_id,  1, ensembl_db_name['homo_sapiens'])
            human_exon_phase = get_exon_phase (cursor, human_exon.exon_id,  1)

            for species in species_list:
                switch_to_db(cursor, ensembl_db_name[species])
                qry      = "select * from sw_exon where maps_to_human_exon_id = %d " % human_exon.exon_id
                sw_exons = search_db(cursor, qry)

                if not sw_exons:
                    #print "no sw exons for ", species
                    continue

                print "\n", species

                ct = 0
                ok = 0
                for sw_exon in sw_exons:
                    ct += 1

                    has_stop = False
                    has_NNN  = False
                    [sw_exon_id, gene_id, start_in_gene, end_in_gene,  maps_to_human_exon_id, exon_seq_id,
                     template_exon_seq_id, template_species,  strand, phase, has_NNN, has_stop, has_3p_ss, has_5p_ss] = sw_exon
                    [exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, 
                     left_flank, right_flank, dna_seq] = get_exon_seq_by_db_id (cursor, exon_seq_id, ensembl_db_name[species])


                    if not has_stop and not has_NNN:
                        ok  += 1

                    ########################
                    # see if the  protein seq matches the quoted boundaries
                    # coding dna sequence:
                    cds = dna_seq[pepseq_transl_start:pepseq_transl_end]
                    translated_cds = Seq(cds).translate().tostring()
                    if not  translated_cds == protein_seq:
                        print Seq(cds).translate().tostring()
                        print protein_seq
                        exit(1)

                    template_db_id = species2genome_db_id (cursor, template_species)

                    [templ_exon_seq_id, templ_protein_seq, templ_pepseq_transl_start, 
                     templ_pepseq_transl_end,  templ_left_flank, templ_right_flank, templ_dna_seq] \
                     = get_exon_seq_by_db_id (cursor, template_exon_seq_id, ensembl_db_name[template_species])

                    ########################
                    # 
                    print "\t     human", human_protein_seq, human_exon_phase
                    print "\t deposited", sw_exon_id, protein_seq
                    print "\t translted", sw_exon_id, translated_cds
                    print "\t  template", templ_protein_seq
                    print "\t  template", template_exon_seq_id, template_species, template_db_id

                tot_ok += ok
                tot += ct

                print species, "total: ", ct, "  ok: ", ok
        exit(1)
    
    print tot, tot_ok

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
