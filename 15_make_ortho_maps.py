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
from   el_utils.map     import  *

from   el_utils.special_gene_sets  import  get_theme_ids
from   el_utils.almt_cmd_generator import  AlignmentCommandGenerator
from   el_utils.config_reader      import  ConfigurationReader
from   el_utils.alignment          import  smith_waterman, exon_aware_smith_waterman
from   alignment import * # C implementation of smith waterman

#########################################
verbose = True


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

    if ( len(seq_human) > len(seq_other)):
        padding = ""
        for i in range  (len(seq_human)-len(seq_other)):
            padding += "-"
        seq_other += padding

    if ( len(seq_other) > len(seq_human)):
        padding = ""
        for i in range  (len(seq_other)-len(seq_human)):
            padding += "-"
        seq_human += padding

    seq_human_no_common_gaps = ""
    seq_other_no_common_gaps = ""
    
    for i in range (len(seq_human)):
        if seq_human[i] == '-' and seq_other[i] == '-': continue
        seq_human_no_common_gaps += seq_human[i]
        seq_other_no_common_gaps += seq_other[i]

    return [seq_human_no_common_gaps, seq_other_no_common_gaps] 

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

    if not cigar_line: 
        return [seq_A_aligned, seq_B_aligned]

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
        pepseq = get_exon_pepseq (cursor, exon)
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


    for gene_id in gene_list:

        ct += 1
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        if verbose: print gene_id, gene2stable(cursor, gene_id), get_description(cursor, gene_id)
        
        # get _all_ exons
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        human_exons = gene2exon_list(cursor, gene_id)
        if (not human_exons):
            print 'no exons for ', gene_id
            continue

        # get rid of the old maps
        map_cleanup (cursor, ensembl_db_name, human_exons)

        # human as its own orthologue - let's be systematic
        maps = self_maps (cursor, ensembl_db_name, human_exons)
        store (cursor, maps, ensembl_db_name)

        # one2one   orthologues
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        known_orthologues      = get_orthos (cursor, gene_id, 'orthologue')

        # not-clear orthologues
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        unresolved_orthologues = get_orthos (cursor, gene_id, 'unresolved_ortho')

        orthologues =  known_orthologues+unresolved_orthologues
        for [ortho_gene_id, ortho_species] in orthologues:
            #if not ortho_species == 'dipodomys_ordii': continue
            ortho_exons = gene2exon_list(cursor, ortho_gene_id, db_name=ensembl_db_name[ortho_species] )
            # in the first round of building the database there will be no 'novel' exons
            # but we put this here just in case we want to rerun the script anyway
            # (the script 21_make_novel_exon_maps.py rebuilds the maps only for the species for which
            # there are some novel exons - this one builds the maps from scratch)
            ortho_exons += get_novel_exons (cursor, ortho_gene_id, 'sw_exon')
            ortho_exons += get_novel_exons (cursor, ortho_gene_id, 'usearch_exon')

            usearch_ct =  len(get_novel_exons (cursor, ortho_gene_id, 'usearch_exon'))
            if usearch_ct: print ortho_species+", usearch: ", usearch_ct

            if not ortho_exons:
                missing_exon_info += 1
                print "\t", ortho_species, "no exon info"
                continue

            maps = make_maps (cursor, ensembl_db_name,  cfg, acg, ortho_species, human_exons, ortho_exons)   
            if not maps:
                missing_seq_info += 1
                print "\t", ortho_species, "no maps"
                continue

            no_maps += len(maps)
            store (cursor, maps, ensembl_db_name)

        if  not ct%100:
            datastring = StringIO.StringIO()
            print >> datastring, "processed ", ct, "genes,  out of ", len(gene_list), "  ",
            print >> datastring, no_maps, " maps;  no_exon_info: ", missing_exon_info, "no_seq_info:", missing_seq_info 
            print datastring.getvalue()
    cursor.close()
    db.close()

    return True


#########################################
def main():
    
    no_threads = 10
    special    = None

    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <set name> <number of threads> " % sys.argv[0]
        exit(1)
    elif len(sys.argv)==3:

        special = sys.argv[1]
        special = special.lower()
        if special == 'none': special = None

        no_threads = int(sys.argv[2])

    local_db = False

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
        gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special )
    else:
        print "using all protein coding genes"
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        

    cursor.close()
    db.close()

    parallelize (no_threads, maps_for_gene_list, gene_list, [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
    main()

#########################################
#########################################
#########################################


'''

    for gene_id in [412667]: #  wls
    for gene_id in [378768]: #  p53
     #for gene_id in [378766]: #  dynein
 

        # COMMENT THIS OUT PERHAPS?
        # get rid of the old maps
        # map_cleanup (cursor, ensembl_db_name, human_exons)



#########################################
def alignment_line (seq_human, seq_other):

    alignment_line = []

    if ( not len(seq_human) ==  len(seq_other) ):
        print "alignment_line: the sequences must be aligned"
        return ""
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



        

'''
