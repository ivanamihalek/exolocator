#!/usr/bin/python


import StringIO
import MySQLdb, commands, re, sys
from hashlib import sha1
from random  import random
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.utils   import  *
from   el_utils.threads import  parallelize
from   el_utils.map     import  get_maps, Map
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.alignment          import smith_waterman, exon_aware_smith_waterman
from   el_utils.special_gene_sets  import *
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
def overlap (start, end, other_start, other_end):
    if ( other_end < start): 
        return False
    elif ( end < other_start):
        return False
    else:
        return True


#########################################
def maps_evaluate (cursor, cfg, ensembl_db_name, human_exons, ortho_exons, aligned_seq, exon_positions):

    maps = []
   
    if len(aligned_seq.keys()) > 2:
        print "right now the mapping implemented for two species only"
        exit (1)

    for species in aligned_seq.keys():
        if species == 'homo_sapiens': continue
        other_species = species
        break

    min_similarity = cfg.get_value('min_accptbl_exon_sim')

    for human_exon_ct in range(len(human_exons)):

        padded_count_human = "{0:03d}".format(human_exon_ct+1)
        if ( not  exon_positions['homo_sapiens'].has_key(padded_count_human) ):
            continue
        [human_start, human_end] = exon_positions['homo_sapiens'][padded_count_human]

        human_exon = human_exons[human_exon_ct]
        old_maps   = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)


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

                if ortho_exons[ortho_exon_ct].is_known == 2:
                    map.source   = 'sw_sharp' 
                elif ortho_exons[ortho_exon_ct].is_known == 3:
                    map.source   = 'usearch' 
                else:
                    map.source   = 'ensembl'

                exon_seq_human   = aligned_seq['homo_sapiens'][human_start:human_end].replace('#','-')
                exon_seq_other   = aligned_seq[other_species][other_start:other_end].replace('#','-')
                [seq_human, seq_other] = pad_the_alnmt (exon_seq_human,human_start,
                                                        exon_seq_other, other_start)
                seq = {'human':seq_human, 'other':seq_other}
                seq = strip_gaps(seq)
                if not seq:  
                    c=inspect.currentframe()
                    print " in %s:%d" % ( c.f_code.co_filename, c.f_lineno)
                    exit(1)

                map.similarity = pairwise_tanimoto(seq['human'], seq['other'])

                if False  and not map.source == 'ensembl':
                    print
                    print other_species, map.source
                    print seq['human']
                    print seq['other']
                    print map.similarity
                    print

                if map.similarity < min_similarity: continue

                ciggy = cigar_line (seq['human'], seq['other'])
                map.cigar_line = ciggy
                                                   

                maps.append(map)

    return maps

########################################
def find_relevant_exons (cursor, all_exons, human):

    relevant_exons = []
    protein_seq    = []

    # 1) choose exons that I need
    for exon in all_exons:
        if (exon.covering_exon > 0):
            continue
        if human and not exon.is_coding:
            continue
        if exon.is_known>1 and not exon.exon_seq_id:
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
def make_maps (cursor, ensembl_db_name, cfg, acg, ortho_species, human_exons, ortho_exons):

    maps = []

    #print "############################## human"
    switch_to_db(cursor,  ensembl_db_name['homo_sapiens'])
    relevant_human_exons = find_relevant_exons (cursor, human_exons, human=True)


    #print "##############################", ortho_species
    switch_to_db(cursor,  ensembl_db_name[ortho_species])
    relevant_ortho_exons = find_relevant_exons (cursor, ortho_exons, human=False)

    if not relevant_ortho_exons:
        print "bleep 2"
        return maps

    human_seq = decorate_and_concatenate (relevant_human_exons)
    ortho_seq = decorate_and_concatenate (relevant_ortho_exons)
   
    if (not human_seq or not ortho_seq):
        print "bleep 3"
        return maps
    
    aligned_seq = {}
    if 0:# python implementation
        [aligned_seq['homo_sapiens'], aligned_seq[ortho_species]] \
            = exon_aware_smith_waterman (human_seq, ortho_seq)
    else: # C implementation
        ret = smith_waterman_context (human_seq, ortho_seq, -5, -3)
        if len(ret) == 2:
            [aligned_seq['homo_sapiens'], aligned_seq[ortho_species]] = ret
        else:
            print ret
            return []

    if (not aligned_seq['homo_sapiens'] or 
        not aligned_seq[ortho_species]):
        print "bleep 4"
        return []

    # find the positions of the exons in the alignment
    exon_positions = {}
    for species, seq in aligned_seq.iteritems():
        # move B to beginning of each exon sequence
        seq = moveB(seq)
        #beginning and end of each exon in the alignment
        exon_positions[species] = find_exon_positions(seq)

    # fill in the actual map values
    maps = maps_evaluate (cursor, cfg, ensembl_db_name, relevant_human_exons, 
                          relevant_ortho_exons, aligned_seq, exon_positions)

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
        ################################
        switch_to_db(cursor,ensembl_db_name['homo_sapiens']) 
        store_or_update (cursor, 'exon_map', fixed_fields, update_fields)

    return True

#########################################
def  map_cleanup (cursor, ensembl_db_name, human_exons):
    
    switch_to_db(cursor,ensembl_db_name['homo_sapiens']) 
    for exon in human_exons:
        qry  = "delete from exon_map where exon_id = %d " % exon.exon_id
        qry += " and exon_known = %d " % exon.is_known
        qry += " and cognate_exon_known > 1 " 
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

    #######################################
    for gene_id in gene_list:

        ct += 1
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        if not ct%10: print ct, "out of ", len(gene_list) 
        
        # get _all_ exons
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        print  gene_id,  gene2stable(cursor, gene_id), get_description (cursor, gene_id)

        human_exons = gene2exon_list(cursor, gene_id)
        if (not human_exons):
            print 'no exons for ', gene_id
            continue
            #sys.exit(1)

        # get rid of the old maps # can't do that here bcs this script is only updating sw exons
        map_cleanup (cursor, ensembl_db_name, human_exons)

        orthologues  = get_orthos (cursor, gene_id, 'orthologue') # get_orthos changes the db pointer
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        orthologues += get_orthos (cursor, gene_id, 'unresolved_ortho')

        ##########
        for [ortho_gene_id, ortho_species] in orthologues:

            #if not ortho_species=='ochotona_princeps': continue
            #print ortho_species, ortho_gene_id, species2genome_db_id (cursor, ortho_species)

            switch_to_db (cursor, ensembl_db_name[ortho_species])

            ortho_exons   = []
            
            ortho_exons += get_novel_exons (cursor, ortho_gene_id, 'sw_exon')
            ortho_exons += get_novel_exons (cursor, ortho_gene_id, 'usearch_exon')

            if not ortho_exons: continue # nothing new here, move on

            ortho_exons +=  get_known_exons (cursor, ortho_gene_id, ortho_species)
            ortho_exons +=  get_predicted_exons (cursor, ortho_gene_id, ortho_species)

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
    
    no_threads = 1
    special    = 'one'


    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <set name> <number of threads> <method>"
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