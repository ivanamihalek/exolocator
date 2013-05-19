#!/usr/bin/python


import StringIO
import MySQLdb, commands, re, sys
from hashlib import sha1
from random  import random
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  *
from   el_utils.utils   import  *
from   el_utils.threads import  parallelize
from   el_utils.map     import  get_maps, Map
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.alignment          import smith_waterman, exon_aware_smith_waterman
from   el_utils.special_gene_sets  import human_genes_w_sw_sharp_annotation, get_theme_ids
from   alignment import * # C implementation of smith waterman


#########################################
def decorate_and_concatenate (exons):
    decorated_seq = ""
    count = 1
    for  exon in exons:
        pepseq = exon.pepseq
        padded_count   = "{0:03d}".format(count)
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

    if ( len(seq_human) > len(seq_other)):
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
def overlap (start, end, other_start, other_end):
    if ( other_end < start): 
        return False
    elif ( end < other_start):
        return False
    else:
        return True


#########################################
def maps_evaluate (cursor, ensembl_db_name, human_exons, ortho_exons, aligned_seq, exon_positions):

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

        human_exon = human_exons[human_exon_ct]
        old_maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)


        for ortho_exon_ct in range(len(ortho_exons)):

            padded_count_ortho = "{0:03d}".format(ortho_exon_ct+1)
            if ( not  exon_positions[other_species].has_key(padded_count_ortho) ):
                continue
            [other_start, other_end] = exon_positions[other_species][padded_count_ortho]

            if ( overlap (human_start, human_end, other_start, other_end) ):
                
                map = Map()
                map.species_1    = 'homo_sapiens'
                map.species_2    = other_species
                
                map.exon_id_1    = human_exons[human_exon_ct].exon_id
                map.exon_id_2    = ortho_exons[ortho_exon_ct].exon_id

                map.exon_known_1 = human_exons[human_exon_ct].is_known
                map.exon_known_2 = ortho_exons[ortho_exon_ct].is_known

                if ortho_exons[ortho_exon_ct].analysis_id == -1:
                    map.source   = 'sw_sharp' 
                elif ortho_exons[ortho_exon_ct].analysis_id == -2:
                    map.source   = 'usearch' 
                else:
                    map.source   = 'ensembl'

                exon_seq_human   = aligned_seq['homo_sapiens'][human_start:human_end].replace('#','-')
                exon_seq_other   =  aligned_seq[other_species][other_start:other_end].replace('#','-')
                [seq_human, seq_other] = pad_the_alnmt (exon_seq_human,human_start,
                                                        exon_seq_other, other_start)

                ciggy = cigar_line (seq_human, seq_other)
                [seq1, seq2] = unfold_cigar_line (seq_human.replace('-',''), seq_other.replace('-',''), ciggy)

                map.cigar_line = ciggy
                map.similarity = pairwise_fract_similarity (seq1, seq2)
                
                # bit of paranoia, but not misplaced:  do we already have a map for this exon by any chance?
                better_map = False
                for old_map in old_maps:
                    if old_map.species_2  == other_species:
                        if  old_map.similarity >  map.similarity:
                            better_map = True
                            break
                if not better_map: maps.append(map)

    return maps

########################################
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
        # local version of the function, to search possibly for sw_exon
        pepseq = get_exon_pepseq (cursor, exon)
        if not pepseq:
            to_remove.append(i)
            continue
        pepseq = pepseq.replace ('X', '')
        if  len (pepseq) < 3:
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
    if 0:# python implementation
        [aligned_seq['homo_sapiens'], aligned_seq[ortho_species]] \
            = exon_aware_smith_waterman (human_seq, ortho_seq)
    else: # C implementation
        [aligned_seq['homo_sapiens'], aligned_seq[ortho_species]] \
            = smith_waterman_context (human_seq, ortho_seq)

    if (not aligned_seq['homo_sapiens'] or 
        not aligned_seq[ortho_species]):
        return []

    # find the positions of the exons in the alignment
    exon_positions = {}
    for species, seq in aligned_seq.iteritems():
        # move B to beginning of each exon sequence
        seq = moveB(seq)
        #beginning and end of each exon in the alignment
        exon_positions[species] = find_exon_positions(seq)

    # fill in the actual map values
    maps = maps_evaluate (cursor, ensembl_db_name, relevant_human_exons, relevant_ortho_exons, aligned_seq, exon_positions)

    return maps
    
#########################################
def store (cursor, maps, ensembl_db_name, source = None):

    for map in maps:
        fixed_fields  = {}
        update_fields = {}
        fixed_fields ['exon_id']              = map.exon_id_1
        fixed_fields ['exon_known']           = map.exon_known_1
        fixed_fields ['cognate_genome_db_id'] = species2genome_db_id(cursor, map.species_2)
        fixed_fields ['cognate_exon_id']      = map.exon_id_2
        fixed_fields ['cognate_exon_known']   = map.exon_known_2
        fixed_fields ['source']               = map.source 
        update_fields['cigar_line']           = map.cigar_line
        update_fields['similarity']           = map.similarity
        #####
        switch_to_db(cursor,ensembl_db_name['homo_sapiens']) 
        store_or_update (cursor, 'exon_map', fixed_fields, update_fields)

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
def exonify(cursor, ensembl_db_name, row_from_sw_exon_table, method):

    exon = Exon()
    
    [sw_exon_id, gene_id, start_in_gene, end_in_gene, human_exon_id,
     exon_seq_id, template_exon_id, template_species,
     strand, phase, has_NNN, has_stop, has_3p_ss, has_5p_ss] = row_from_sw_exon_table

    human_coding = is_coding (cursor, human_exon_id, ensembl_db_name['homo_sapiens'])
    
    exon.gene_id             = gene_id
    exon.exon_id             = sw_exon_id
    exon.start_in_gene       = start_in_gene
    exon.end_in_gene         = end_in_gene
    exon.canon_transl_start  = -1
    exon.canon_transl_end    = -1
    exon.exon_seq_id         = exon_seq_id
    exon.strand              = strand
    exon.phase               = phase
    exon.is_coding           = 1 if (human_coding and not has_stop) else 0
    exon.is_canonical        = 0
    exon.is_constitutive     = 0
    exon.covering_exon       = -1
    exon.covering_exon_known = -1 

    if method == 'sw_sharp':
        exon.is_known        =  2 # arbitrary index, used to indicate sw_sharp findings here
        exon.analysis_id     = -1
    elif  method == 'usearch':
        exon.is_known        =  3 # arbitrary index, used to indicate usearch  findings 
        exon.analysis_id     = -2
    else:
        exon.is_known        =  0 
        exon.analysis_id     =  0
        

    return exon


#########################################
def maps_for_gene_list(gene_list, db_info):
    

    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        cfg    = ConfigurationReader()
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

    #######################################
    for gene_id in gene_list:

        ct += 1
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        if not ct%10: print ct, "out of ", len(gene_list) 
        
        # get _all_ exons
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        print ct, len(gene_list),  gene_id,  gene2stable(cursor, gene_id), get_description (cursor, gene_id)

        human_exons = gene2exon_list(cursor, gene_id)
        if (not human_exons):
            print 'no exons for ', gene_id
            continue
            #sys.exit(1)

        # get rid of the old maps # can't do that here bcs this script is only updating sw exons
        #map_cleanup (cursor, ensembl_db_name, human_exons)
        ##########
        for ortho_species, db_name in ensembl_db_name.iteritems():
            if ortho_species == 'homo_sapiens': continue

            ortho_exons   = []
            ortho_gene_id = None
            for human_exon in human_exons:
                switch_to_db (cursor, db_name)

                # novel exon
                for table in ['sw_exon', 'usearch_exon']:
                    method = 'sw_sharp' if table=='sw_exon' else 'usearch'
                    qry  = "select * from %s " % table
                    qry += " where maps_to_human_exon_id = %d " % human_exon.exon_id
                    rows = search_db (cursor, qry)
                    if not rows: continue

                    for row in rows:
                        # turn them to actual exon objects
                        exon = exonify (cursor, ensembl_db_name, row, method)
                        if not ortho_gene_id:
                            ortho_gene_id = exon.gene_id
                        elif not ortho_gene_id == exon.gene_id:
                            print "funny error in gene assignment ..."
                            continue
                        #if not exon.is_coding: continue ? whats this?
                        ortho_exons.append(exon)

            if not ortho_exons: continue
      
            # other exons from this species
            ortho_exons += gene2exon_list(cursor, ortho_gene_id, db_name=db_name)

            maps = make_maps (cursor, ensembl_db_name, cfg, acg, ortho_species, human_exons, ortho_exons) 
            if not maps:
                print "\t", ortho_species, "no maps"
                continue
            
            switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
            store (cursor, maps, ensembl_db_name)
                
    cursor.close()
    db.close()

    return True


#########################################
def main():
    
    no_threads = 10
    special    = 'telomere_maintenance'

    local_db   = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    if special:
        print "using", special, "set"
        gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special )
    else:
        print "using all protein coding genes that have an sw# patch"
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        gene_list = human_genes_w_sw_sharp_annotation (cursor, ensembl_db_name)

    cursor.close()
    db.close()

    parallelize (no_threads, maps_for_gene_list, gene_list[0:15000], [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
    main()

'''

    for gene_id in [412667]: #  wls
    for gene_id in [378768]: #  p53
     #for gene_id in [378766]: #  dynein
 

        # COMMENT THIS OUT PERHAPS?
        # get rid of the old maps
        # map_cleanup (cursor, ensembl_db_name, human_exons)
        

'''
