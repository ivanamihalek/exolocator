#!/usr/bin/python

import MySQLdb
import commands
from random import choice

from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  *
from   el_utils.exon    import  Exon
from   el_utils.threads import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.special_gene_sets  import  get_theme_ids
from   el_utils.config_reader      import ConfigurationReader

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna



#########################################
def check_ccds (cursor, transcript_stable_id):

    ccds = ""

    qry  = "select dna_align_feature.hit_name "
    qry += "from dna_align_feature, transcript, transcript_supporting_feature "
    qry += "   where dna_align_feature.dna_align_feature_id =  transcript_supporting_feature.feature_id "
    qry += "   and transcript_supporting_feature.feature_type ='dna_align_feature' "
    qry += "   and transcript_supporting_feature.transcript_id =transcript.transcript_id "
    qry += "   and transcript.stable_id = '%s' " % transcript_stable_id

    rows = search_db(cursor, qry)

    if not rows:
        return ccds


    for row in rows:
        if 'CCDS' in row[0]:
            ccds = row[0]

    return  ccds

#########################################
def transcript_id2exon_ids (cursor, transcript_id):

    exon_ids = []
    qry = "select exon_id from exon_transcript "
    qry += " where transcript_id = %d " %  transcript_id
    rows   = search_db (cursor, qry)
    if (not rows):
        return []
    for row in rows:
        exon_ids.append(row[0])

    return exon_ids

#########################################
def check_alt_splices (cursor, species, ensembl_db_name):

    
    print "############################"
    print 'checking alt splicing in ', species

    qry = "use " + ensembl_db_name[species]
    search_db(cursor, qry)
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    no_cover_and_no_seq = 0  
    no_cover_and_no_dna_seq = 0  
    no_seq_info_in_database = 0
    all_ok   = 0   
    cover_already_present  = 0
    cover_not_in_exon_set  = 0     
    cov_exon_not_in_ccds   = 0 
    genes_w_ccds           = 0
    tot_exons = 0
        
    #for gene_id in gene_ids[:100]:
    #for gene_id in [413198]:
    for count in range(1000):
        gene_id = choice (gene_ids)

        #print  gene_id,  gene2stable(cursor, gene_id), get_description (cursor, gene_id)
        transcript_ids = get_transcript_ids(cursor, gene_id)

        tr_w_ccds = []
        for [tr_id, tr_stable] in transcript_ids:
            ccds = check_ccds (cursor, tr_stable)
            if not ccds: continue
            tr_w_ccds.append([tr_id, tr_stable])

        if not tr_w_ccds: continue

        genes_w_ccds += 1

        # get all exons for this gene
        all_exons    = gene2exon_list (cursor, gene_id)
        
        exons_w_ccds = set([]) # get the unique_ids

        # find exons which are on the ccds list
        for [tr_id, tr_stable] in tr_w_ccds:
            exon_ids =  transcript_id2exon_ids (cursor, tr_id)
            exons_w_ccds.update( set(exon_ids))
           
        # for these exons check sequence
        is_known=1
        for exon_id in exons_w_ccds:
            tot_exons += 1
            exon = get_exon      (cursor, exon_id, is_known)
            seq  = get_exon_seqs (cursor, exon_id, is_known)
            if not seq:
                #print  gene_id,  gene2stable(cursor, gene_id), get_description (cursor, gene_id)
                #print exon_id, " no seq", "covered: ", exon.covering_exon
                #exit (1)
                no_seq_info_in_database += 1
                continue
            [exon_seq_id, protein_seq, pepseq_transl_start, 
             pepseq_transl_end, left_flank, right_flank, dna_seq] = seq

            #print exon_id, exon_seq_id, 
            #print "  %7d   %7d  %7d  " % (exon.start_in_gene, exon.end_in_gene,  exon.covering_exon),
            
            if exon.covering_exon < 0:
                #print protein_seq
                if not protein_seq:
                    no_cover_and_no_seq += 1
                elif not dna_seq:
                    no_cover_and_no_dna_seq += 1
                else:
                    all_ok += 1
            else:
                if exon.covering_exon_known and exon.covering_exon in exons_w_ccds:
                    #print " <<<< "
                    cover_already_present += 1
                else:
                    all_exon_ids =  map(lambda exon: exon.exon_id, all_exons)
                    if not exon.covering_exon in all_exon_ids:
                        cover_not_in_exon_set += 1
                        print "cover_not_in_exon_set: "
                        print  gene_id,  gene2stable(cursor, gene_id), get_description (cursor, gene_id)
                        print  exon_id, exon_seq_id, 
                        print "  %7d   %7d  %7d  " % (exon.start_in_gene, exon.end_in_gene,  exon.covering_exon)
                        #print covering_exon
                        print " ************"
                        for e in all_exons:
                            print e
                        exit (1)

                    elif not exon.covering_exon in exons_w_ccds:
                        cov_exon_not_in_ccds += 1
                        #print "covering exon is not in ccds set "
                       


    print "genes_w_ccds", genes_w_ccds
    print "tot_exons", tot_exons
    print "no_seq_info_in_database ",   no_seq_info_in_database
    print "all_ok",  all_ok 
    print "cover_already_present ", cover_already_present
    print "no_cover_and_no_seq   ", no_cover_and_no_seq 
    print "no_cover_and_no_dna_seq   ", no_cover_and_no_dna_seq 
    print "cov_exon_not_in_ccds",   cov_exon_not_in_ccds
    print "cover_not_in_exon_set ", cover_not_in_exon_set

    return True



#########################################
def main():

    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    # human and mouse are the only two species that have CCDs info
    for species in ['homo_sapiens']:
        check_alt_splices (cursor, species, ensembl_db_name)



    cursor.close()
    db    .close()

#########################################
if __name__ == '__main__':
    main()

