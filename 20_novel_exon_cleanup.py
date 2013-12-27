#!/usr/bin/python


import MySQLdb
from tempfile import NamedTemporaryFile
from operator import itemgetter
from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.utils   import  *
from el_utils.config_reader      import ConfigurationReader
from el_utils.special_gene_sets  import *
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

    canonical_human_exons.sort(key=lambda exon: exon.start_in_gene)
    return canonical_human_exons

#########################################
def check_translation_start (mitochondrial, left_flank, dna_seq, template_dna_seq, templ_protein_seq):

    correction    = None
    phase         = 0
    left_flank_ok = False


    if not templ_protein_seq[0]=='M':
       return [left_flank_ok, correction, phase]    

    if mitochondrial:
        pepseq = Seq(dna_seq[:3]).translate(table="Vertebrate Mitochondrial").tostring()
    else:
        pepseq = Seq(dna_seq[:3]).translate().tostring()
    
    if pepseq[0] == 'M':
        left_flank_ok = True
        return [left_flank_ok, correction, phase]    

    if len(dna_seq) < 3*len(templ_protein_seq):
        shift_range_left  = 12 + abs(len(dna_seq)- 3*len(templ_protein_seq))/3*3
        shift_range_right = 12
    else:
        shift_range_left  = 12 
        shift_range_right = 12 + abs(len(dna_seq)- 3*len(templ_protein_seq))/3*3
        

    concat    = left_flank +dna_seq
    old_joint = len(left_flank)
   
    while shift_range_left > len(left_flank):
        shift_range_left -= 3
    while shift_range_right > len(dna_seq):
         shift_range_right -= 3

    for correction in range(0, -shift_range_left+3, -3):
        codon = concat[old_joint-correction-3:old_joint-correction]
        if mitochondrial:
            pepseq = Seq(codon).translate(table="Vertebrate Mitochondrial").tostring()
        else:
            pepseq = Seq(codon).translate().tostring()
        if pepseq and pepseq[0] == 'M':
            left_flank_ok = True
            return [left_flank_ok, correction, phase]    

    for correction in range(0, shift_range_right-3, +3):
        codon = concat[old_joint+correction:old_joint+correction+3]
        if mitochondrial:
            pepseq = Seq(codon).translate(table="Vertebrate Mitochondrial").tostring()
        else:
            pepseq = Seq(codon).translate().tostring()
        if pepseq and pepseq[0] == 'M':
            left_flank_ok = True
            return [left_flank_ok, correction, phase] 
  


    return [left_flank_ok, correction, phase]  

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


    cmd    = acg.generate_maxentscan_cmd (3, tmp_in_file.name)
    blah   = commands.getoutput (cmd).split("\n")
    scores = []
    for sc in blah:
        try:
            num = float(sc)
        except:
            num = None
        if not num == None:    scores.append(num)
    tmp_in_file.close()

    # find the position of the max score   
    index_of_max_score = None
    for i in range(len(scores)):
        sc = scores[i]
        if max_score < sc: 
            max_score =sc
            index_of_max_score = i
    if not index_of_max_score or index_of_max_score > len(shifts):
        return [left_flank_ok, correction, phase, max_score]


    if max_score < 3.0: return [left_flank_ok, correction, phase, max_score]

    left_flank_ok = True
    max_score_shift = shifts[index_of_max_score]

    # for novel exons, translation starts and ends at the exon boundaries

    ########### important!
    # if no score is positive - abandon, this is probably not major class, 
    # or some other "cryptic" splice site
    
    # we need to move our presumtpive end ov exon by
    # this much in order to get the most likely splicing signal
    correction = max_score_shift
    phase      = correction%3
    if phase < 0: phase = abs(phase)

    left_flank_ok = ( max_score >=  6.0)

    return [left_flank_ok, correction, phase, max_score]

#########################################
def check_right_flank(acg, right_flank, dna_seq, template_dna_seq):

    correction  = None
    phase       = 0
    right_flank_ok = False
    max_score   = -10
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
    blah   = commands.getoutput (cmd).split("\n")
    scores = []
    for sc in blah:
        try:
            num = float(sc)
        except:
            num = None
        if not num == None:    scores.append(num)



    tmp_in_file.close()

    if not scores:
        return [right_flank_ok, correction, phase, max_score]

    # find the position of the max score   
    index_of_max_score = None
    for i in range(len(scores)):
        sc = float(scores[i])
        if max_score < sc: 
            max_score =sc
            index_of_max_score = i

    if not index_of_max_score or index_of_max_score > len(shifts):
        # allscores are negative
        return [right_flank_ok, correction, phase, max_score]

    max_score_shift = shifts[index_of_max_score]
    if max_score < 3.0: return [right_flank_ok, correction, phase, max_score]
    right_flank_ok = True
   

    correction = max_score_shift
    phase      = correction%3
    if phase < 0: phase = (3-abs(phase))%3


    right_flank_ok = (max_score >= 6.0)

    return [right_flank_ok, correction, phase, max_score]

#########################################
def check_coordinates_in_the_gene (cursor, cfg, acg, ensembl_db_name, species, sw_exon, exon_dna_seq):

    [sw_exon_id, gene_id, start_in_gene, end_in_gene,  maps_to_human_exon_id, exon_seq_id,
     template_exon_seq_id, template_species,  strand, phase, end_phase, has_NNN, has_stop, 
     has_3p_ss, has_5p_ss] = sw_exon


    gene_coords =  get_gene_coordinates (cursor, gene_id, ensembl_db_name[species])  
    [gene_seq_region_id, gene_start, gene_end, gene_strand] = gene_coords


    if ( gene_strand >  0 ):
        region_start = gene_start + start_in_gene
        region_end   = gene_start + end_in_gene

    else:
        region_end   = gene_end - start_in_gene
        region_start = gene_end - end_in_gene
        

    qry  = "select name, file_name "
    qry += " from seq_region where seq_region_id = %d " % int(gene_seq_region_id)
    rows = search_db(cursor, qry)
    [name, file_names] = rows[0]
    filename = get_best_filename(file_names)
    fasta    = get_fasta (acg, species, name, filename, gene_strand, region_start, region_end)
    qry_seq = "".join(fasta.splitlines()[1:])

    delete = True

    resultstr = usearch (cfg, acg, exon_dna_seq, qry_seq, delete)
    match = parse_usearch_output (resultstr)

    if not match:# should do something in this case
        print "no match - the stored in-gene coordinates ",
        print " (%s, %s) probably wrong" % (start_in_gene, end_in_gene)
        print sw_exon
        #exit(1)
        return []
    else:
        [match_start, match_end, template_start, template_end, 
         aligned_target_seq, aligned_template_seq] = match
        
    if ( gene_strand >  0 ):
        start_in_gene_corrected = region_start + match_start - gene_start + 1
        end_in_gene_corrected   = region_start + match_end - gene_start   + 1
    else:
        start_in_gene_corrected = match_start - (region_end - gene_end) + 1
        end_in_gene_corrected   = match_end   - (region_end - gene_end) + 1
        
 

    return [start_in_gene, end_in_gene]


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

    tot         = 0
    tot_ok      = 0
    for human_gene_id in gene_list:

        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id   = gene2stable(cursor, human_gene_id)
        description = get_description (cursor, human_gene_id)


        mitochondrial = is_mitochondrial(cursor, human_gene_id)

        #print "#############################################"
        #print human_gene_id, stable_id, get_description (cursor, human_gene_id)

        human_exons = get_ok_human_exons (cursor,ensembl_db_name,  human_gene_id)
        
        for human_exon in human_exons:
            [exon_seq_id, human_protein_seq, pepseq_transl_start, pepseq_transl_end, 
             left_flank, right_flank, dna_seq] = get_exon_seqs (cursor, human_exon.exon_id,  1,
                                                                ensembl_db_name['homo_sapiens'])
            human_exon_phase = get_exon_phase (cursor, human_exon.exon_id,  1)

            first_exon = (human_exons.index(human_exon) == 0)

            for species in mammals: # maxentscan does not work for fish 
                 
                 for table in ['sw_exon','usearch_exon']:
                     switch_to_db(cursor, ensembl_db_name[species])
                     qry      = "select * from %s where maps_to_human_exon_id = %d " % (table, human_exon.exon_id)
                     sw_exons = search_db(cursor, qry)

                     if not sw_exons:
                         #print  "human_exon: ", human_exon.exon_id, "no", table,  "for", species
                         continue
                     ct = 0
                     ok = 0
                     for sw_exon in sw_exons:
                         ct += 1

                         has_stop = False
                         has_NNN  = False

                         [sw_exon_id, gene_id, start_in_gene, end_in_gene,  maps_to_human_exon_id, exon_seq_id,
                          template_exon_seq_id, template_species,  strand, phase, end_phase, has_NNN, has_stop, 
                          has_3p_ss, has_5p_ss] = sw_exon

                         tot +=1
             
                         exon_seqs =  get_exon_seq_by_db_id (cursor, exon_seq_id, ensembl_db_name[species])
                         if not exon_seqs: 
                             print "exon seq not stored "
                             continue

                         [exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, 
                          left_flank, right_flank, dna_seq] = exon_seqs

                         len_ok   =  (pepseq_transl_end-pepseq_transl_start) == len (dna_seq)
                         if not len_ok:
                             # if it is not the case, then make it be so
                             left_flank += dna_seq[:pepseq_transl_start]
                             right_flank = dna_seq[pepseq_transl_end:] + right_flank
                             dna_seq   = dna_seq[pepseq_transl_start:pepseq_transl_end]
                             pepseq_transl_start = 0
                             pepseq_transl_end   = len(dna_seq)

                         phase_ok =  (len (dna_seq)%3 == 0)
                         if not phase_ok:
                             phase = len(dna_seq)%3
                             cds = dna_seq[phase:]
                             pepseq_corrected = Seq(cds).translate().tostring()
                             if pepseq_corrected == protein_seq:
                                 left_flank += dna_seq[:phase]
                                 dna_seq     = dna_seq[phase:]
                             else:
                                 cds = dna_seq[:-phase]
                                 pepseq_corrected = Seq(cds).translate().tostring()

                                 if pepseq_corrected == protein_seq:
                                     right_flank += dna_seq[-phase:] + right_flank
                                     dna_seq     = dna_seq[:-phase]
                                 else: 
                                     print "no match ..."
                                     continue # don't want to shut-off the pipeline here
                             
                             pepseq_transl_start = 0
                             pepseq_transl_end   = len(dna_seq)
                       

                         # retrieve the template
                         template_db_id = species2genome_db_id (cursor, template_species)

                         [templ_exon_seq_id, templ_protein_seq, templ_pepseq_transl_start, 
                          templ_pepseq_transl_end,  templ_left_flank, templ_right_flank, templ_dna_seq] \
                          = get_exon_seq_by_db_id (cursor, template_exon_seq_id, ensembl_db_name[template_species])

                         correction = 0
                         phase      = 0
                         end_phase  = 0

                         # if this is the first exon, check if we are starting from methionine
                         if first_exon:
                             [left_flank_ok, correction, phase] = \
                                 check_translation_start (mitochondrial, left_flank, dna_seq, templ_dna_seq, templ_protein_seq)
                         # see if the left splice site is ok
                         else:
                             [left_flank_ok, correction, phase, max_score] = \
                                 check_left_flank (acg, left_flank, dna_seq, templ_dna_seq)

                         ########################
                         # 
                         # see if the right splice site is ok
                         [right_flank_ok, end_correction, end_phase, end_max_score] = \
                             check_right_flank(acg, right_flank, dna_seq, templ_dna_seq)

                         #if not left_flank_ok and not right_flank_ok: continue

                         pepseq_corrected = ""
                         new_left_flank   = ""
                         new_right_flank  = ""
                         new_dna_seq      = ""
                         if left_flank_ok: 
                             offset = (3-phase)%3
                             if correction:
                                 if correction > 0:
                                     new_dna_seq     =              dna_seq[correction:]
                                     new_left_flank  = left_flank + dna_seq[:correction]
                                 else:
                                     # correction is negative, therefore left_flank[correction:] is the tail of left_flank
                                     new_dna_seq     = left_flank[correction:]+ dna_seq
                                     new_left_flank  = left_flank[:correction]
                             else:
                                 new_left_flank = left_flank

                             pepseq_transl_start = offset

                         else:
                             new_left_flank = left_flank

                         if right_flank_ok:
                             if not new_dna_seq: new_dna_seq = dna_seq
                             if end_correction:
                                 if end_correction < 0:
                                     new_right_flank = new_dna_seq[end_correction:] + right_flank
                                     new_dna_seq     = new_dna_seq[:end_correction]
                                 else:
                                     # correction is negative, therefore right_flank[correction:] is the tail of right_flank
                                     new_right_flank  = right_flank[end_correction:]
                                     new_dna_seq     += right_flank[:end_correction]
                             else:
                                 new_right_flank = right_flank
                             pepseq_transl_end  = len(new_dna_seq) 
                             pepseq_transl_end -= end_phase
                         else:
                             new_right_flank = right_flank

                         # if only one flank is ok, use that side to decide if there is a phase on the other
                         if left_flank_ok and not right_flank_ok:
                             end_phase = (pepseq_transl_end-pepseq_transl_start)%3
                             pepseq_transl_end -= end_phase
                         
                         if right_flank_ok and not left_flank_ok:
                             phase = (pepseq_transl_end-pepseq_transl_start)%3
                             pepseq_transl_start += phase
 
                         # check that the lengths match
                         has_stop = None
                         if new_dna_seq: 
                             len_old = len(    left_flank +     dna_seq +     right_flank)
                             len_new = len(new_left_flank + new_dna_seq + new_right_flank)
                             if not len_old==len_new:
                                 
                                 print len_old, len_new
                                 print correction, end_correction
                                 print map (len, [left_flank, dna_seq, right_flank])
                                 print map (len, [new_left_flank, new_dna_seq, new_right_flank])
                                 continue
                             cds = new_dna_seq[pepseq_transl_start:pepseq_transl_end]
                             if mitochondrial:
                                 pepseq_corrected = Seq(cds).translate(table="Vertebrate Mitochondrial").tostring()
                             else:
                                 pepseq_corrected = Seq(cds).translate().tostring()
                             if '*' in pepseq_corrected:
                                 has_stop = 1
                             else:
                                 has_stop = 0

                         if has_stop and not '*' in protein_seq: continue # abort, abort

                         if True:
                             print "#############################################"
                             print human_gene_id, stable_id, description
                             print species, table

                             print "\t  template", template_exon_seq_id, template_species, template_db_id
                             print "\t  template left flank", templ_left_flank, templ_dna_seq[0:3]
                             print "\t           left flank", left_flank, dna_seq[0:3]
                             print "\t          ",  left_flank_ok, correction, phase,
                             if not first_exon:
                                 print max_score
                             else:
                                 print
                             print "\t  template right flank", templ_dna_seq[-3:],templ_right_flank
                             print "\t           right flank", dna_seq[-3:], right_flank
                             print "\t          ", right_flank_ok, end_correction, end_phase, end_max_score

                             print "\t     human", human_protein_seq, human_exon.exon_id, human_exon_phase 
                             print "\t  template", templ_protein_seq
                             print "\t deposited", protein_seq
                             if pepseq_corrected: print "\t corrected", pepseq_corrected
                             
                         if new_dna_seq:
                             if (pepseq_transl_end-pepseq_transl_start)%3:
                                 print "length not divisible by 3 "
                                 print pepseq_transl_start, pepseq_transl_end
                                 print phase, end_phase
                                 print len(new_dna_seq)
                                 print "%%%%% "
                                 continue
                         else:
                             new_dna_seq = dna_seq 

                         #########################################################
                         # 18_find_exons is sometimes messing up the coordinates 
                         # I do not know why
                         ret = check_coordinates_in_the_gene (cursor, cfg, acg, ensembl_db_name, 
                                                              species, sw_exon, new_dna_seq)
                         if not ret: 
                             print "\t coordinate check failed"
                             continue
                         [start_in_gene_corrected, end_in_gene_corrected] = ret
                             

                         #########################################################
                         # update the *_exon and exon_seq tables accordingly
                         switch_to_db(cursor, ensembl_db_name[species])
                         
                     
                         qry = "update %s set " % table
                         
                         set_fields = ""
                         if not start_in_gene_corrected == start_in_gene:
                            if set_fields: set_fields += ", "
                            set_fields += " start_in_gene = %d  " % start_in_gene_corrected

                         if not end_in_gene_corrected == end_in_gene:
                            if set_fields: set_fields += ", "
                            set_fields += " end_in_gene = %d  " % end_in_gene_corrected

                         if not has_stop is None:
                             if set_fields: set_fields += ", "
                             set_fields += " has_stop  = %d" % has_stop
                            
                         if left_flank_ok:
                             if set_fields: set_fields += ", "
                             set_fields += " phase = %d,  " % phase
                             if first_exon:
                                 set_fields += " has_3p_ss = '%s' "  % ("first exon; starts with M")
                             else:
                                 set_fields += " has_3p_ss = '%s' "  % ("me_score="+str(max_score))

                         if right_flank_ok:
                             if set_fields: set_fields += ", "
                             set_fields += " end_phase = %d,  " % end_phase
                             set_fields += " has_5p_ss = '%s' "  % ("me_score="+str(end_max_score))
                         
                         
                         qry += set_fields + " where exon_id=%d" %  sw_exon_id
                         
                         if set_fields:
                             search_db(cursor, qry)
 
                         # update exon sequence
                         if pepseq_corrected: 
                             # we might have changed our mind as to what is the cDNA seq, and what is flanking
                             qry  = "update exon_seq set "
                             qry += " protein_seq = '%s', " % pepseq_corrected
                             qry += " dna_seq = '%s',     " % new_dna_seq
                             qry += " left_flank = '%s',  " % new_left_flank
                             qry += " right_flank = '%s',       " % new_right_flank
                             qry += " pepseq_transl_start = %d, " % pepseq_transl_start
                             qry += " pepseq_transl_end   = %d  " % pepseq_transl_end
                             table_id = 2 if table=='sw_exon' else 3
                             qry += " where exon_id=%d and is_known=%d" % (sw_exon_id, table_id)
                             search_db(cursor, qry)
                            
                         # gene2exon --> have to go back to 07_gene2exon for that
                         tot_ok += 1

    
    cursor.close()
    db.close()

#########################################
def main():
    
    no_threads = 1
    special    = 'test'

    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <set name> <number of threads>" % sys.argv[0]
        exit(1)
    elif len(sys.argv)==3:

        special = sys.argv[1]
        special = special.lower()
        if special == 'none': special = None

        no_threads = int(sys.argv[2])

 
    local_db   = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)


    print '======================================='
    print sys.argv[0]
    if special:
        print "using", special, "set"
        if special == 'complement':
            gene_list = get_complement_ids(cursor, ensembl_db_name, cfg)
        else:
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
