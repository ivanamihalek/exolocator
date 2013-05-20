#!/usr/bin/python


import MySQLdb
from tempfile import NamedTemporaryFile
from operator import itemgetter
from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.config_reader      import ConfigurationReader
from el_utils.special_gene_sets  import human_genes_w_sw_sharp_annotation, get_theme_ids
from el_utils.almt_cmd_generator import AlignmentCommandGenerator
from el_utils.threads import  parallelize
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
def check_left_flank (acg, left_flank, dna_seq, template_dna_seq):

    correction    = None
    phase         = 0
    left_flank_ok = False
    max_score     = -10
    # the things to check:
    # 1) we should not get out of range - in particular, are we hitting NNNN?
    # 2) how big is the difference in length between the template an the proposed exon
    # 3) splicing model for fish 
    # 4) are we lloking at U12 splice site by any chance

    if len(dna_seq) < len(template_dna_seq):
        shift_range_left  = 10 + abs(len(dna_seq)-len(template_dna_seq))
        shift_range_right = 10 
    else:
        shift_range_left  = 10 
        shift_range_right = 10 + abs(len(dna_seq)-len(template_dna_seq))
        
    splice_site = left_flank +dna_seq
    splice_zero = len(left_flank)

    shifts = []
    list_of_shifted_seqs = ""
    for shift in range(-shift_range_left, shift_range_right):
        start = splice_zero+shift-20
        end   = splice_zero+shift+3
        if start < 0: continue
        if end   > len(splice_site): continue
        if 'N' in splice_site[start:end]: continue
        list_of_shifted_seqs += splice_site[start:end] +"\n"
        shifts.append(shift)

    if not list_of_shifted_seqs:
       return [left_flank_ok, correction, phase, max_score]
 

    tmp_in_file = NamedTemporaryFile(delete=True)
    tmp_in_file.write(list_of_shifted_seqs)
    tmp_in_file.flush()


    cmd     = acg.generate_maxentscan_cmd (3, tmp_in_file.name)
    scores  = commands.getoutput (cmd).split("\n")
    tmp_in_file.close()



    # find the position of the max score        
    index_of_max_score =  min( enumerate(scores), key=itemgetter(1))[0]  
    max_score = float(scores[index_of_max_score])
    if index_of_max_score > len(shifts):
        print cmd
        print list_of_shifted_seqs
        exit(1)
    max_score_shift = shifts[index_of_max_score]
    if max_score < 6.0: return [left_flank_ok, correction, phase, max_score]
    left_flank_ok = True

    ########### important!
    # if no score is positive - abandon, this is probably not major class, 
    # or some other "cryptic" splice site
    
    # we need to move our presumtpive end ov exon by
    # this much in order to get the most likely splicing signal
    if max_score_shift == 0: 
        correction = None
        phase      = 0
    elif max_score_shift == -1:
        correction = None
        phase      = 2
    elif max_score_shift == -2:
        correction = None
        phase      = 1
    else: # we seem to have missed by a couple of residues
        correction  = max_score_shift
        sign        = correction/abs(correction)
        correction  = abs(correction)
        phase       = correction%3
        correction  = correction/3*3
        correction *= sign
        if sign < 0: phase = (3-phase)%3

    return [left_flank_ok, correction, phase, max_score]

#########################################
def check_right_flank(acg, right_flank, dna_seq, template_dna_seq):

    correction  = None
    phase       = 0
    right_flank_ok = False
    max_score   = -1
    tmp_in_file = NamedTemporaryFile(delete=True)

    shift_range = 10 + abs(len(dna_seq)-len(template_dna_seq))
    if len(dna_seq) < len(template_dna_seq):
        shift_range_left  = 10 
        shift_range_right = 10 + abs(len(dna_seq)-len(template_dna_seq))
    else:
        shift_range_left  = 10 + abs(len(dna_seq)-len(template_dna_seq))
        shift_range_right = 10 
        
    splice_site = dna_seq + right_flank
    splice_zero = len(dna_seq)

    shifts = []
    list_of_shifted_seqs = ""
    for shift in range(-shift_range_left, shift_range_right):
        start = splice_zero+shift-3
        end   = splice_zero+shift+9
        if start < 0: continue
        if end   > len(splice_site): continue
        if 'N' in splice_site[start:end]: continue
        list_of_shifted_seqs += splice_site[start:end]+"\n"
        shifts.append(shift)

    if not list_of_shifted_seqs:
       return [right_flank_ok, correction, phase, max_score]

    tmp_in_file = NamedTemporaryFile(delete=False)
    tmp_in_file.write(list_of_shifted_seqs)
    tmp_in_file.flush()

    cmd     = acg.generate_maxentscan_cmd (5, tmp_in_file.name)
    scores  = commands.getoutput (cmd).split("\n")
    tmp_in_file.close()

    if not scores:
        return [right_flank_ok, correction, phase, max_score]

    # find the position of the max score        
    index_of_max_score =  min( enumerate(scores), key=itemgetter(1))[0]  


    max_score       = float(scores[index_of_max_score])
    if index_of_max_score > len(shifts):
        print cmd
        print list_of_shifted_seqs
        exit(1)
    max_score_shift = shifts[index_of_max_score]
    if max_score < 6.0: return [right_flank_ok, correction, phase, max_score]
    right_flank_ok = True
   

    if max_score_shift   == 0: 
        correction = None
        phase      = 0
    elif max_score_shift == 1:
        correction = None
        phase      = 1    # this is endphase
    elif max_score_shift == 2:
        correction = None 
        phase      = 2
    else: # we seem to have missed by a couple of residues
        correction = max_score_shift
        sign       = correction/abs(correction)
        correction = abs(correction)
        phase      = correction%3
        correction = correction/3*3
        correction *= sign
        if sign < 0: phase = (3-phase)%3

    return [right_flank_ok, correction, phase, max_score]


#########################################
def exon_cleanup(gene_list, db_info):


    [local_db, ensembl_db_name] = db_info
    
    
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
        acg = AlignmentCommandGenerator()
    else:
        db  = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    # find db ids and common names for each species db
    all_species, ensembl_db_name = get_species (cursor)


    for human_gene_id in gene_list:
    #for human_gene_id in [397176]:
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id   = gene2stable(cursor, human_gene_id)
        description = get_description (cursor, human_gene_id)
        #print "#############################################"
        #print human_gene_id, stable_id, get_description (cursor, human_gene_id)

        human_exons = get_ok_human_exons (cursor,ensembl_db_name,  human_gene_id)
        tot         = 0
        tot_ok      = 0

        for human_exon in human_exons:
            [exon_seq_id, human_protein_seq, pepseq_transl_start, pepseq_transl_end, 
             left_flank, right_flank, dna_seq] = get_exon_seqs (cursor, human_exon.exon_id,  1, ensembl_db_name['homo_sapiens'])
            human_exon_phase = get_exon_phase (cursor, human_exon.exon_id,  1)

            for species in all_species:
            #for species in ['pelodiscus_sinensis']:
                 for table in ['sw_exon','usearch_exon']:
                     switch_to_db(cursor, ensembl_db_name[species])
                     qry      = "select * from %s where maps_to_human_exon_id = %d " % (table, human_exon.exon_id)
                     sw_exons = search_db(cursor, qry)

                     if not sw_exons:
                         #print "no", table,  "for", species
                         continue
                     ct = 0
                     ok = 0
                     for sw_exon in sw_exons:
                         ct += 1

                         has_stop = False
                         has_NNN  = False

                         [sw_exon_id, gene_id, start_in_gene, end_in_gene,  maps_to_human_exon_id, exon_seq_id,
                          template_exon_seq_id, template_species,  strand, phase, end_phase, has_NNN, has_stop, has_3p_ss, has_5p_ss] = sw_exon
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
                         # see if the left splice site is ok
                         [left_flank_ok, correction, phase, max_score] = \
                             check_left_flank (acg, left_flank, dna_seq, templ_dna_seq)

                         ########################
                         # 
                         # see if the right splice site is ok
                         [right_flank_ok, end_correction, end_phase, end_max_score] = \
                             check_right_flank(acg, right_flank, dna_seq, templ_dna_seq)

                         cds  = ""
                         pepseq_corrected = ""
                         if left_flank_ok and correction:
                             if correction > 0:
                                 cds = dna_seq[pepseq_transl_start+correction:pepseq_transl_end]
                             else:
                                 # correction is negative, therefore left_flank[correction:] is the tail of left_flank
                                 cds = left_flank[correction:]+ dna_seq[pepseq_transl_start:pepseq_transl_end]
                         
                         if right_flank_ok and end_correction:
                             if not cds: cds = dna_seq[pepseq_transl_start:pepseq_transl_end]
                             if end_correction < 0:
                                 cds = cds[:end_correction]
                             else:
                                 # correction is negative, therefore left_flank[correction:] is the tail of left_flank
                                 cds += right_flank[:end_correction]
                                 
                         if cds: 
                             pepseq_corrected = Seq(cds).translate().tostring()
                             if '*' in pepseq_corrected:
                                 has_stop = True

                         print "#############################################"
                         print human_gene_id, stable_id, description
                         print species, table

                         print "\t     human", human_protein_seq, human_exon_phase
                         print "\t deposited", protein_seq
                         print "\t translted", translated_cds
                         print "\t  template", templ_protein_seq
                         print "\t  template", template_exon_seq_id, template_species, template_db_id
                         print "\t  template left flank", templ_left_flank, templ_dna_seq[0:3]
                         print "\t           left flank", left_flank, dna_seq[0:3]
                         print "\t          ",  left_flank_ok, correction, phase, max_score
                         print "\t  template right flank", templ_dna_seq[-3:],templ_right_flank
                         print "\t           right flank", dna_seq[-3:], right_flank
                         print "\t          ", right_flank_ok, end_correction, end_phase, end_max_score

                         if pepseq_corrected: print "\t corrected", pepseq_corrected

                         #########################################################
                         # update the *_exon and exon_seq tables accordingly
                         switch_to_db(cursor, ensembl_db_name[species])
                         qry = "update %s set " % table
                         if left_flank_ok:
                             qry += " phase = %d,  " % phase
                             qry += " has_3p_ss = '%s', "  % ("me_score="+str(max_score))
                         if right_flank_ok:
                             qry += " end_phase = %d,  " % end_phase
                             qry += " has_5p_ss = '%s' "  % ("me_score="+str(end_max_score))
                         qry += " where exon_id=%d" %  sw_exon_id
                         
                         print  ensembl_db_name[species]
                         print qry

                         qry = "update %s set " % table
                         if left_flank_ok:
                             qry += " phase = %d,  " % phase
                             qry += " has_3p_ss = '%s', "  % ("me_score="+str(max_score))
                         if right_flank_ok:
                             qry += " end_phase = %d,  " % end_phase
                             qry += " has_5p_ss = '%s' "  % ("me_score="+str(end_max_score))
                         qry += " where exon_id=%d" %  sw_exon_id
                          

                         exit(1)

                     tot_ok += ok
                     tot += ct

            #print species, "total: ", ct, "  ok: ", ok
    
    print tot, tot_ok

    cursor.close()
    db.close()

#########################################
def main():
    
    no_threads = 1
    special    = 'test'

    local_db   = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)


    if special:
        print "using", special, "set"
        gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special )
    else:
        print "using all protein coding genes"
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    cursor.close()
    db.close()

    parallelize (no_threads, exon_cleanup, gene_list, [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
    main()
