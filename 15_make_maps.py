#!/usr/bin/python


import MySQLdb, commands, re, sys
from hashlib import sha1
from random  import random
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
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
def alignment_line (seq_human, seq_other):

    alignment_line = []

    if ( not len(seq_human) ==  len(seq_other) ):
        print "alignment_line: the sequences must be aligned"
        exit(1)
    else:
        length = len(seq_human)

    for i in range(length):
        if not seq_human[i] == "-" and  not seq_other[i] == "-":
            alignment_line.append ("|")

        elif seq_human[i] == "-" and  seq_other[i] == "-":
            #pass
            alignment_line.append ("-")

        elif (seq_human[i]  == "-" ):
            alignment_line.append ("A")

        elif (seq_other[i]  == "-" ):
            alignment_line.append ("B")
    return  "".join(alignment_line)


#########################################
def cigar_line (seq_human, seq_other):

    cigar_line     = []

    alignment_line = []

    if ( not len(seq_human) ==  len(seq_other) ):
        print "alignment_line:  the seqeunces must be aligned"
        exit(1)
    else:
        length = len(seq_human)

    if not length:
        print "zero length sequence (?)"
        exit (1)

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
def maps_evaluate (human_exons, ortho_exons, aligned_seq, exon_positions):

    maps = []

    if len(aligned_seq.keys()) > 2:
        print "right now the mapping implemented for two species only"
        exit (1)

    for species in aligned_seq.keys():
        if species == 'homo_sapiens': continue
        other_species = species
        break

    for human_exon_ct in range(len(human_exons)):

        padded_count_human = "{0:03d}".format(human_exon_ct+1)
        if ( not  exon_positions['homo_sapiens'].has_key(padded_count_human) ):
            continue
        [human_start, human_end] = exon_positions['homo_sapiens'][padded_count_human]


        for ortho_exon_ct in range(len(ortho_exons)):

            padded_count_ortho = "{0:03d}".format(ortho_exon_ct+1)
            if ( not  exon_positions[other_species].has_key(padded_count_ortho) ):
                continue
            [other_start, other_end] = exon_positions[other_species][padded_count_ortho]

            if ( human_start <= other_start <= human_end or
                 human_start <= other_end <= human_end):
                
                map = Map()
                map.species_1 = 'homo_sapiens'
                map.species_2 = other_species
                
                map.exon_id_1 = human_exons[human_exon_ct].exon_id
                map.exon_id_2 = ortho_exons[ortho_exon_ct].exon_id

                map.exon_known_1 = human_exons[human_exon_ct].is_known
                map.exon_known_2 = ortho_exons[ortho_exon_ct].is_known

                exon_seq_human = aligned_seq['homo_sapiens'][human_start:human_end]
                exon_seq_other = aligned_seq[other_species][other_start:other_end]
                [seq_human, seq_other] = pad_the_alnmt (exon_seq_human,human_start,
                                                        exon_seq_other, other_start)

                ciggy = cigar_line (seq_human, seq_other)
                [seq1, seq2] = unfold_cigar_line (seq_human.replace('-',''), seq_other.replace('-',''), ciggy)

                map.cigar_line = ciggy
                map.similarity = fract_identity (ciggy)
                
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
def mafft_align (cfg, acg, seq1, seq2):

    aligned_seqs = []

    # generate rand string for the name
    name = sha1(str(random())).hexdigest()[:15]
    # output
    fastafile = "{0}/{1}.fa".format  (cfg.dir_path['scratch'], name)
    output_fasta (fastafile, ['seq1', 'seq2'], {'seq1':seq1, 'seq2':seq2}) 
    # align
    mafftcmd = acg.generate_mafft_command (fastafile)
    almt     = commands.getoutput(mafftcmd)

    seq = ""
    for line in almt.split('\n'):
        if '>' in line: 
            if seq:
                aligned_seqs.append(seq)
                seq = ""
        else:
            seq += line
    aligned_seqs.append(seq)
    
    commands.getoutput("rm "+fastafile)
 

    return aligned_seqs

#########################################
def  pairwise_fract_identity (seqs):
    
    fract_identity = 0.0
    [seq1, seq2]   = seqs
    if ( not len(seq1)):
        return fract_identity

    for i in range(len(seq1)):
        if (seq1[i] == '-'): continue
        if seq1[i] == seq2[i]: fract_identity += 1.0
    
    fract_identity /= float(len(seq1))
    return fract_identity

#########################################
def make_maps (cursor, ensembl_db_name, cfg, acg, ortho_species, human_exons, ortho_exons):

    maps = []

    #print "############################## human"
    switch_to_db(cursor,  ensembl_db_name['homo_sapiens'])
    relevant_human_exons = find_relevant_exons (cursor, human_exons)
    #print "##############################", ortho_species
    switch_to_db(cursor,  ensembl_db_name[ortho_species])
    relevant_ortho_exons = find_relevant_exons (cursor, ortho_exons)

    human_seq = decorate_and_concatenate (relevant_human_exons)
    ortho_seq = decorate_and_concatenate (relevant_ortho_exons)
    
    if (not human_seq or not ortho_seq):
        return maps
    
    aligned_seq = {}
    if 0:
        [aligned_seq['homo_sapiens'], aligned_seq[ortho_species]] \
            = exon_aware_smith_waterman (human_seq, ortho_seq)
    else: # C implementation
        [aligned_seq['homo_sapiens'], aligned_seq[ortho_species]] \
            = smith_waterman_context (human_seq, ortho_seq)

    if (not aligned_seq['homo_sapiens'] or 
        not aligned_seq[ortho_species]):
        return maps # this is empty

    # find the positions of the exons in the alignment
    exon_positions = {}
    for species, seq in aligned_seq.iteritems():
        # move B to beginning of each exon sequence
        seq = moveB(seq)
        #beginning and end of each exon in the alignment
        exon_positions[species] = find_exon_positions(seq)

    # fill in the actual map values
    maps = maps_evaluate (relevant_human_exons, relevant_ortho_exons, aligned_seq, exon_positions)

    return maps
    
#########################################
def store (cursor, maps, ensembl_db_name):

    for map in maps:
        fixed_fields  = {}
        update_fields = {}
        fixed_fields ['exon_id']              = map.exon_id_1
        fixed_fields ['exon_known']           = map.exon_known_1
        fixed_fields ['cognate_genome_db_id'] = species2genome_db_id(cursor, map.species_2)
        fixed_fields ['cognate_exon_id']      = map.exon_id_2
        fixed_fields ['cognate_exon_known']   = map.exon_known_2
        update_fields['cigar_line']           = map.cigar_line
        update_fields['similarity']           = map.similarity
        update_fields['source']               = 'ensembl'
        #####
        switch_to_db(cursor,ensembl_db_name['homo_sapiens']) 
        store_or_update (cursor, 'exon_map', fixed_fields, update_fields)

    return True


#########################################
def gene_has_a_map (cursor, ensembl_db_name, human_exons):

    has_a_map = False
    for human_exon in human_exons:
        if ( not human_exon.is_canonical or  not human_exon.is_coding): continue
        maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
        if maps:
            has_a_map = True
            break

    return has_a_map


#########################################
def maps_for_gene_list(gene_list, db_info):
    

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
    

    for gene_id in gene_list:

        ct += 1
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        #print gene_id, gene2stable(cursor, gene_id), get_description (cursor, gene_id)
        
        # get _all_ exons
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        human_exons = gene2exon_list(cursor, gene_id)
        if (not human_exons):
            print 'no exons for ', gene_id
            sys.exit(1)
        # check maps
        if gene_has_a_map (cursor, ensembl_db_name, human_exons):
            continue

        # one2one   orthologues
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        known_orthologues      = get_orthos (cursor, gene_id, 'orthologue')
        # not-clear orthologues
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        unresolved_orthologues = get_orthos (cursor, gene_id, 'unresolved_ortho')
        # 
        for [ortho_gene_id, ortho_species] in known_orthologues+unresolved_orthologues:
            ortho_exons = gene2exon_list(cursor, ortho_gene_id, db_name=ensembl_db_name[ortho_species] )
            if not ortho_exons:
                missing_exon_info += 1
                print "\t", ortho_species, "no exon info"
                continue
            #print "\t", ortho_species, "making maps ..."
            maps = make_maps (cursor, ensembl_db_name,  cfg, acg, ortho_species, human_exons, ortho_exons)   
            if not maps:
                missing_seq_info += 1
                #print "\t", ortho_species, "no maps"
                continue
            no_maps += len(maps)
            
            store (cursor, maps, ensembl_db_name)

        if (not ct%10):
            print "processed ", ct, "genes,  out of ", len(gene_list), "  ",
            print no_maps, " maps;   no_exon_info: ", missing_exon_info , "no_seq_info:", missing_seq_info 
 
    cursor.close()
    db.close()

    return True

#########################################
def main():
    
    no_threads = 7

    local_db = False

    if local_db:
        db = connect_to_mysql()
    else:
        db = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)


    species                        = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_list                      = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
    cursor.close()
    db.close()

    parallelize (no_threads, maps_for_gene_list, gene_list[15000:], [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
    main()

'''

    for gene_id in [412667]: #  wls
    for gene_id in [378768]: #  p53
     #for gene_id in [378766]: #  dynein
 


'''
