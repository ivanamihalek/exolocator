#!/usr/bin/python


import StringIO
import MySQLdb, commands, re, sys
from hashlib import sha1
from random  import random, choice
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  *
from   el_utils.utils   import  *
from   el_utils.threads import  parallelize
from   el_utils.map     import  get_maps, Map
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.alignment          import smith_waterman, exon_aware_smith_waterman
from   alignment import * # C implementation of smith waterman


#########################################
def decorate_and_concatenate (exons):
    decorated_seq = ""
    count = 1
    for  exon in exons:
        pepseq = exon.pepseq
        padded_count = "{0:03d}".format(count)
        decorated_seq += 'B'+padded_count+pepseq+'Z'
        count += 1

    return decorated_seq

#########################################
def moveB (seq):
    seqlist = list(seq)
    begin_pattern = re.compile("B\d{3}\-+")

    for begin_label in begin_pattern.finditer(seq):
        start =  begin_label.start()
        end   =  begin_label.end()
        label =  seq[start:end]
        for i in range(4):
            seqlist[start+i] = '-'
        for i in range(4):
            seqlist[end-4+i] = label[i]

    return "".join(seqlist)

#########################################
def find_exon_positions(seq):

    exon_position = {}
    
    exon_pattern = re.compile("B.*?Z")
    for match in exon_pattern.finditer(seq):
        start       = match.start()
        end         = match.end()
        exon_seq_no = seq[start+1:start+4]
        exon_position[exon_seq_no] = [start+4, end-1]  #B+3 digits on one end,  Z on the other

    return exon_position

#########################################
def  pad_the_alnmt (exon_seq_human, human_start, exon_seq_other, other_start):
    
    seq_human = ""
    seq_other = ""

    padding = ""
    if ( human_start > other_start):
        for i in range (human_start-other_start):
            padding += "-"
    seq_human = padding + exon_seq_human


    padding = ""
    if ( other_start > human_start):
        for i in range (other_start-human_start):
            padding += "-"
    seq_other = padding + exon_seq_other

    if ( len(seq_human) >  len(seq_other)):
        padding = ""
        for i in range  (len(seq_human)-len(seq_other)):
            padding += "-"
        seq_other += padding

    if ( len(seq_other) >  len(seq_human)):
        padding = ""
        for i in range  (len(seq_other)-len(seq_human)):
            padding += "-"
        seq_human += padding


    return [seq_human, seq_other] 

    

#########################################
def cigar_line (seq_human, seq_other):

    cigar_line     = []

    alignment_line = []

    if ( not len(seq_human) ==  len(seq_other) ):
        print "alignment_line:  the seqeunces must be aligned"
        return ""
    else:
        length = len(seq_human)

    if not length:
        print "zero length sequence (?)"
        return ""

    for i in range(length):
        if not seq_human[i] == "-" and  not seq_other[i] == "-":
            alignment_line.append ("M")

        elif seq_human[i] == "-" and  seq_other[i] == "-":
            pass
            #alignment_line.append ("-")

        elif (seq_human[i]  == "-" ):
            alignment_line.append ("A")

        elif (seq_other[i]  == "-" ):
            alignment_line.append ("B")

            
    prev_char = alignment_line[0]
    count     = 1
    for i in range(1,len(alignment_line)):
        if ( alignment_line[i] == prev_char):
            count += 1
        else:
            cigar_line.append( "{0}{1}".format(count, prev_char))
            prev_char = alignment_line[i]
            count     = 1
                               
    cigar_line.append("{0}{1}".format(count, prev_char))

    return  "".join(cigar_line)


#########################################
def unfold_cigar_line (seq_A, seq_B, cigar_line):

    seq_A_aligned = ""
    seq_B_aligned = ""


    char_pattern = re.compile("\D")
    a_ct     = 0
    b_ct     = 0
    prev_end = 0

    for match in char_pattern.finditer(cigar_line):
        this_start       = match.start()
        no_repeats = int(cigar_line[prev_end:this_start])
        alignment_instruction = cigar_line[this_start]
        prev_end = match.end()

        if alignment_instruction == 'M':
            seq_A_aligned += seq_A[a_ct:a_ct+no_repeats]
            a_ct         += no_repeats
            seq_B_aligned += seq_B[b_ct:b_ct+no_repeats]
            b_ct  += no_repeats

        elif alignment_instruction == 'A':
            seq_A_aligned += '-'*no_repeats 
            seq_B_aligned += seq_B[b_ct:b_ct+no_repeats]
            b_ct  += no_repeats
            
        elif alignment_instruction == 'B':
            seq_A_aligned += seq_A[a_ct:a_ct+no_repeats]
            a_ct          += no_repeats
            seq_B_aligned += '-'*no_repeats 

    return [seq_A_aligned, seq_B_aligned]


#########################################
def  fract_identity (cigar_line):

    fraction = 0

    char_pattern = re.compile("\D")
    total_length = 0
    common       = 0
    prev_end     = 0

    for match in char_pattern.finditer(cigar_line):
        this_start       = match.start()
        no_repeats = int(cigar_line[prev_end:this_start])
        alignment_instruction = cigar_line[this_start]
        prev_end = match.end()

        total_length += no_repeats
        if alignment_instruction == 'M':
            common += no_repeats
     
    if total_length:
        fraction = common/float(total_length)
        
    return  fraction

#########################################
def overlap (start, end, other_start, other_end):
    if ( other_end < start): 
        return False
    elif ( end < other_start):
        return False
    else:
        return True


#########################################
def maps_evaluate (template_exons, para_exons, aligned_seq, exon_positions):

    maps = []
   
    if len(aligned_seq.keys()) > 2:
        print "right now the mapping implemented for two species only"
        return []


    for template_exon_ct in range(len(template_exons)):

        padded_count_template = "{0:03d}".format(template_exon_ct+1)
        if ( not  exon_positions['template'].has_key(padded_count_template) ):
            continue
        [template_start, template_end] = exon_positions['template'][padded_count_template]


        for para_exon_ct in range(len(para_exons)):

            padded_count_para = "{0:03d}".format(para_exon_ct+1)
            if ( not  exon_positions['paralogue'].has_key(padded_count_para) ):
                continue
            [other_start, other_end] = exon_positions['paralogue'][padded_count_para]

            if ( overlap (template_start, template_end, other_start, other_end) ):
                
                map = Map()
                map.species_1     = 'template'
                map.species_2     = 'paralogue'
                
                map.exon_id_1     = template_exons[template_exon_ct].exon_id
                map.exon_id_2     = para_exons[para_exon_ct].exon_id

                map.exon_known_1  = template_exons[template_exon_ct].is_known
                map.exon_known_2  = para_exons[para_exon_ct].is_known

                exon_seq_template = aligned_seq['template'][template_start:template_end]
                exon_seq_other    = aligned_seq['paralogue'][other_start:other_end]
                [seq_template, seq_other] = pad_the_alnmt (exon_seq_template,template_start,
                                                           exon_seq_other, other_start)

                ciggy = cigar_line (seq_template, seq_other)
                [seq1, seq2] = unfold_cigar_line (seq_template.replace('-',''), seq_other.replace('-',''), ciggy)

                map.cigar_line = ciggy
                map.similarity = pairwise_tanimoto (seq1, seq2)
                
                maps.append(map)


    return maps

#########################################
def find_relevant_exons (cursor, all_exons):

    relevant_exons = []
    protein_seq    = []

    # 1) choose exons that I need
    for exon in all_exons:
        if (not exon.is_coding or  exon.covering_exon > 0):
            continue
        relevant_exons.append(exon)

    # 2) sort them by their start position in the gene
    to_remove = []
    relevant_exons.sort(key=lambda exon: exon.start_in_gene)
    for i in range(len(relevant_exons)):
        exon   = relevant_exons[i]
        pepseq = get_exon_pepseq (cursor, exon.exon_id, exon.is_known)
        if not pepseq:
            to_remove.append(i)
            continue
        pepseq = pepseq.replace ('X', '')
        if  not pepseq:
            to_remove.append(i)
        else:
            exon.pepseq = pepseq

    for i in range (len(to_remove)-1, -1, -1):
        del relevant_exons[to_remove[i]]
 

    return relevant_exons


#########################################
def make_para_maps (cursor, ensembl_db_name, cfg, acg, template_exons, para_exons):

    maps = []

    relevant_template_exons = find_relevant_exons (cursor, template_exons)
    #print "relevant template: ", map(lambda exon: exon.exon_id, relevant_template_exons)
    relevant_para_exons     = find_relevant_exons (cursor, para_exons)
    #print "relevant para:     ", map(lambda exon: exon.exon_id, relevant_para_exons)

    template_seq = decorate_and_concatenate (relevant_template_exons)
    para_seq     = decorate_and_concatenate (relevant_para_exons)
    
    if (not template_seq or not para_seq):
        return maps
    
    aligned_seq = {}
    [aligned_seq['template'], aligned_seq['paralogue']] \
            = smith_waterman_context (template_seq, para_seq, -5, -3)

    if (not aligned_seq['template'] or 
        not aligned_seq['paralogue']):
        return []

    # find the positions of the exons in the alignment
    exon_positions = {}
    for name, seq in aligned_seq.iteritems():
        # move B to beginning of each exon sequence
        seq = moveB(seq)
        # find beginning and end of each exon in the alignment
        exon_positions[name] = find_exon_positions(seq)

    # fill in the actual map values
    maps = maps_evaluate (relevant_template_exons, relevant_para_exons, aligned_seq, exon_positions)

    return maps
    
#########################################
def store (cursor, maps):

    for map in maps:
        fixed_fields  = {}
        update_fields = {}
        fixed_fields ['exon_id']              = map.exon_id_1
        fixed_fields ['exon_known']           = map.exon_known_1
        #fixed_fields ['cognate_genome_db_id'] = species2genome_db_id(cursor, map.species_2)
        fixed_fields ['cognate_exon_id']      = map.exon_id_2
        fixed_fields ['cognate_exon_known']   = map.exon_known_2
        update_fields['cigar_line']           = map.cigar_line
        update_fields['similarity']           = map.similarity
        update_fields['source']               = 'ensembl'
        #####
        store_or_update (cursor, 'para_exon_map', fixed_fields, update_fields)

    return True

#########################################
def  map_cleanup (cursor, ensembl_db_name, human_exons):
    
    switch_to_db(cursor,ensembl_db_name['homo_sapiens']) 
    for exon in human_exons:
        qry  = "delete from exon_map where exon_id = %d " % exon.exon_id
        qry += " and exon_known = %d " % exon.is_known
        rows = search_db (cursor, qry, verbose=False)


    return True


#########################################
def gene_has_a_para_map (cursor, species, ensembl_db_name, template_exons):

    has_a_map = False
    for template_exon in template_exons:
        maps = get_maps(cursor, ensembl_db_name, template_exon.exon_id,
                        template_exon.is_known, species=species, table='para_exon_map')
        if maps:
            has_a_map = True
            break

    return has_a_map


#########################################
def make_para_exon_maps(species_list, db_info):
    
    verbose = True

    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        cfg      = ConfigurationReader()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    missing_exon_info = 0
    missing_seq_info  = 0
    ct                = 0
    no_maps           = 0
    

    for species in species_list:
        print
        print "############################"
        print  species
        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)
        para_table  = 'paralogue'

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        
        
        for gene_id in gene_ids:
            ct += 1
            if not ct%100: print "\t", species, ct, " out of ", len(gene_ids)
            if verbose: 
                print
                print gene_id, gene2stable(cursor, gene_id), get_description (cursor, gene_id)

            # get the paralogues - only the representative for  the family will have this 
            paralogues = get_paras (cursor, gene_id)  
            if not paralogues:
                if verbose:  print "\t not a template or no paralogues"
                continue

            if verbose:  print "paralogues: ", paralogues

            # get _all_ exons
            template_exons = gene2exon_list(cursor, gene_id)
            if (not template_exons):
                if verbose: print 'no exons for ', gene_id
                continue
            if verbose:
                print
                print "template_exons"
                for ex in template_exons:
                    print ex.exon_id
                print

            # check maps
            # if gene_has_a_para_map (cursor, species, ensembl_db_name, template_exons):
            #     if verbose: print "\t has a map"
            #     continue

            if verbose: print "\t no map found - making new one"

            # COMMENT THIS OUT PERHAPS? yeah, unless you know exactly what you're doing
            # get rid of the old maps
            # map_cleanup (cursor, ensembl_db_name, human_exons)


            # 
            for para_gene_id  in paralogues:
                description =  get_description (cursor, para_gene_id)
                para_exons = gene2exon_list(cursor, para_gene_id)
                if not para_exons:
                    missing_exon_info += 1
                    if verbose: print "\t",description, "no exon info"
                    continue
                if verbose: print "\t", description, "making maps ..."
                maps = make_para_maps (cursor, ensembl_db_name,  cfg, acg,  template_exons, para_exons)   
                if not maps:
                    missing_seq_info += 1
                    if verbose: print "\t", description, "no maps"
                    continue

                if verbose: print "\t", description, "maps ok"
                no_maps += len(maps)
 

                store (cursor, maps)
                exit(1)


    cursor.close()
    db.close()

    return True

#########################################
def main():
    
    no_threads = 1

    local_db = False

    if local_db:
        db = connect_to_mysql()
    else:
        db = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    parallelize (no_threads, make_para_exon_maps, all_species, [local_db, ensembl_db_name])

    return True

#########################################
if __name__ == '__main__':
    main()

'''

    for gene_id in [412667]: #  wls
    for gene_id in [378768]: #  p53
     #for gene_id in [378766]: #  dynein
 


'''
