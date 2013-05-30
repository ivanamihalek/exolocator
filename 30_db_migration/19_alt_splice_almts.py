#!/usr/bin/python

import MySQLdb
import commands
from random import choice

from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.utils   import  *
from   el_utils.ensembl import  *
from   el_utils.exon    import  Exon
from   el_utils.threads import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.special_gene_sets  import  get_theme_ids
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.translation   import translate

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


verbose = True

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
def get_gene_sequence (cursor, acg, species, gene_id):

    sequence = ""
    ret = get_primary_seq_info (cursor, gene_id, species)
    if not ret: return sequence

    [seq_name, file_names, seq_region_start, seq_region_end, 
     seq_region_strand, is_mitochondrial] = ret
   
    sequence = extract_gene_seq(acg, species, seq_name, file_names, seq_region_strand,  
                                seq_region_start, seq_region_end)

    return [seq_region_strand, is_mitochondrial, sequence]

#########################################
def  cleanup_endphase (cursor, exon):


    qry = "select phase, end_phase from exon where exon_id =  %d " % exon.exon_id
    rows  = search_db (cursor, qry)
    if not rows:
        exon.phase     = 0
        exon.end_phase = 0
    else:
        exon.phase     = rows[0][0]
        exon.end_phase = rows[0][1]


    #exon.end_phase = (exon.end_in_gene -(exon.start_in_gene+ (3-exon.phase)%3 )) % 3

#########################################
def alt_splice_almt (cursor, cfg, acg, species, ensembl_db_name):

    flank_length = 10

    print "############################"
    print 'checking alt splicing in ', species

    qry = "use " + ensembl_db_name[species]
    search_db(cursor, qry)
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    if species == 'homo_sapiens':
        spec_short = 'HSA'
    else:
        spec_short = 'MMU'

    outdir   = "{0}/alt/{1}".format(cfg.dir_path['afs_dumps'], spec_short)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ########################################
    ########################################
    ########################################

    #gene_ids.reverse()

    for gene_id in gene_ids:
    #for gene_id in [429349]:
    #for count in range(1000):
        #gene_id = choice (gene_ids)

        stable_gene_id = gene2stable(cursor, gene_id)
        if verbose: print  gene_id, stable_gene_id, get_description (cursor, gene_id)
        transcript_ids = get_transcript_ids(cursor, gene_id)

        tr_w_ccds = []
        for [tr_id, tr_stable] in transcript_ids:
            ccds = check_ccds (cursor, tr_stable)
            if not ccds: continue
            tr_w_ccds.append([tr_id, tr_stable])

        if not tr_w_ccds: continue

        # get all exons for this gene
        all_exons    = gene2exon_list (cursor, gene_id)
        
        exons_w_ccds = set([]) # get the unique_ids

        # find exons which are on the ccds list
        for [tr_id, tr_stable] in tr_w_ccds:
            exon_ids =  transcript_id2exon_ids (cursor, tr_id)
            exons_w_ccds.update( set(exon_ids))
           
        # for these exons check sequence
        is_known = 1
        bad_exon = set([])
        for exon_id in exons_w_ccds:
            exon = get_exon      (cursor, exon_id, is_known)
            seq  = get_exon_seqs (cursor, exon_id, is_known)
            if not seq:
                bad_exon.add(exon_id)
                continue
            [exon_seq_id, protein_seq, pepseq_transl_start, 
             pepseq_transl_end, left_flank, right_flank, dna_seq] = seq

            if exon.covering_exon < 0:
                if not dna_seq:
                     bad_exon.add(exon_id)
            else:
                if exon.covering_exon_known and exon.covering_exon in exons_w_ccds:
                    pass
                else:
                    all_exon_ids =  map(lambda exon: exon.exon_id, all_exons)
                    if not exon.covering_exon in all_exon_ids:
                        bad_exon.add(exon_id)
                        
        # which transcripts seem to be completely ok?
        if verbose: print  "reconstructing alt splice almts for "
        if verbose: print  gene_id,  gene2stable(cursor, gene_id), get_description (cursor, gene_id)
        if verbose: print "there are ", len(tr_w_ccds), " transscripts with ccds"

        # get the gene_sequence
        [seq_region_strand, mitochondrial, gene_sequence] = get_gene_sequence (cursor, acg, species,  gene_id)
        output_seq    = {}
        global_boundaries = []
        local_boundaries  = {}



        # sort exons by the order in which they appear in the gene
        all_exons.sort(key=lambda exon: exon.start_in_gene)


        # a bit of a cleanup
        for exon in all_exons:
            cleanup_endphase (cursor, exon)


        # check if any of the translations is complete:
        no_ok_transcripts = True
        for [tr_id, tr_stable] in tr_w_ccds:
            tr_exon_ids =  transcript_id2exon_ids (cursor, tr_id)
            if bad_exon & set(tr_exon_ids): continue

            if verbose: print tr_stable, " ok "
            no_ok_transcripts = False

        if no_ok_transcripts:
            if verbose: print " no ok transcripts found"
            continue

        # main loop
        cary = "" # for patching up codons split by intron
        for [tr_id, tr_stable] in tr_w_ccds:
            tr_exon_ids =  transcript_id2exon_ids (cursor, tr_id)
            if bad_exon & set(tr_exon_ids): continue

            # translation is from where to where?
            ret = get_translation_coords (cursor, tr_id)
            if not ret:
                continue
            [seq_start, start_exon_id, seq_end, end_exon_id] = ret
            for exon in all_exons:
                if exon.exon_id == start_exon_id: start_exon=exon
                if exon.exon_id == end_exon_id:   end_exon=exon
           
            transl_start_in_gene = start_exon.start_in_gene + seq_start
            transl_end_in_gene   =   end_exon.start_in_gene + seq_end
                

            local_boundaries[tr_stable] = []
            output_seq[tr_stable] = "-"*len(gene_sequence)
            output_seq[tr_stable+"_pep"] = "-"*len(gene_sequence)
            transl_end = ""

            for exon in all_exons:
                if not exon.exon_id in tr_exon_ids: continue
                
                start       = exon.start_in_gene
                start_flank = exon.start_in_gene - flank_length
                if start_flank  < 0: 
                    start_flank  = 0
                else:
                    if not start_flank-1 in global_boundaries: global_boundaries.append(start_flank-1)
                    local_boundaries[tr_stable].append(start_flank)

                end       = exon.end_in_gene
                end_flank = exon.end_in_gene + flank_length
                if end_flank > len(gene_sequence): 
                    end_flank = len(gene_sequence)
                else:
                    if not end_flank in global_boundaries: global_boundaries.append(end_flank)
                    local_boundaries[tr_stable].append(end_flank)
                
                tmp_dna  = output_seq[tr_stable][:start_flank]  + gene_sequence[start_flank:start].lower()
                tmp_dna += gene_sequence[start:end]
                tmp_dna += gene_sequence[end:end_flank].lower() + output_seq[tr_stable][end_flank:]

                output_seq[tr_stable] = tmp_dna


                #################################################
                # now try and handle translation to protein
                prev_transl_end = transl_end

                # where does translation start:
                if exon.end_in_gene < transl_start_in_gene:
                    transl_start = -1
                elif exon.exon_id == start_exon_id:
                    # if this is the first exon, the transl start given above
                    transl_start =  exon.start_in_gene+seq_start-1
                else:
                    # otherwise it is the exon start - except that if this is not the
                    # first exon and the codon is split, we want to start with the
                    # translation of the stitched up exon
                    transl_start = exon.start_in_gene
                    start_flank  = exon.phase

                # where does translation end: 
                if exon.start_in_gene >  transl_end_in_gene:
                    transl_end  = -1
                elif exon.exon_id == end_exon_id:
                    # if this is the first exon, the transl start given above
                    transl_end = exon.start_in_gene+seq_end
                else:
                    # otherwise it is the exon start - except that if this is not the
                    # first exon and the codon is split, we want to start with the
                    # translation of the stitched up exon
                    transl_end  = exon.end_in_gene - exon.end_phase+1
                    end_flank   = exon.end_phase


                if transl_start < 0 or transl_end < 0 :
                     continue

                if exon.phase > 0 and  prev_transl_end:
                    cary = gene_sequence[prev_transl_end:prev_transl_end+exon.phase]
                else:
                    cary = ""

                [phase, pepseq] = translate (cary + gene_sequence[transl_start:transl_end], 
                                             0,  mitochondrial, strip_stop = False)
                prev_transl_end = transl_end

                pepseq_padded = ""
                for aa in pepseq:
                    pepseq_padded += "-"+aa+"-"

                pepseq_name = tr_stable+"_pep"

                tmp_pep  = output_seq[pepseq_name][:transl_start-len(cary)] 
                tmp_pep += pepseq_padded
                tmp_pep += output_seq[pepseq_name][transl_end:]

                output_seq[pepseq_name] = tmp_pep



        global_boundaries.sort()
        for [tr_id, tr_stable] in tr_w_ccds:
            seq =  output_seq[tr_stable]
            tmp_seq   = ""
            prev_bdry = 0
            for bdry in global_boundaries:
                tmp_seq += seq[prev_bdry:bdry] 
                if bdry >= len(seq): continue
                if bdry in local_boundaries[tr_stable]:
                    marker = "-Z-"
                else:
                    marker = "---" 
                tmp_seq += marker 

                prev_bdry = bdry

            output_seq[tr_stable] = tmp_seq

            pepseq_name = tr_stable+"_pep"
            seq =  output_seq[pepseq_name]
            tmp_seq   = ""
            prev_bdry = 0
            for bdry in global_boundaries:
                tmp_seq += seq[prev_bdry:bdry] 
                if bdry >= len(seq): continue
                if bdry in local_boundaries[tr_stable]: # note here
                    marker = "-Z-"
                else:
                    marker = "---" 
                tmp_seq += marker 

                prev_bdry = bdry

            output_seq[pepseq_name] = tmp_seq


        output_seq = strip_gaps(output_seq)


        # define the order in which we  want the sequences output
        name_order = []
        for [tr_id, tr_stable] in tr_w_ccds:
           pepseq_name = tr_stable+"_pep"
           name_order.append (pepseq_name)
           name_order.append (tr_stable)


        afa_fnm  = "{0}/{1}.afa".format(outdir, stable_gene_id)
        ret = output_fasta (afa_fnm, name_order, output_seq)

        print afa_fnm

    return True



#########################################
def main():

    local_db = False

    if local_db:
        db  = connect_to_mysql()
        acg = AlignmentCommandGenerator()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

  
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    # human and mouse are the only two species that have CCDs info
    for species in [ 'homo_sapiens', 'mus_musculus']:
        alt_splice_almt (cursor, cfg, acg, species, ensembl_db_name)



    cursor.close()
    db    .close()

#########################################
if __name__ == '__main__':
    main()

