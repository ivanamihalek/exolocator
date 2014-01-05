
import os
from ensembl     import *
from mysql       import *
from el_specific import *
from exon        import Exon
from   py_alignment          import  smith_waterman, exon_aware_smith_waterman
from   alignment import * # C implementation of smith waterman

#########################################################
class Map:    # this particular map is between exons

    #################################################
    # print
    def __str__ (self):

        printstr = ""

        for attr, value in self.__dict__.iteritems():
            if ( not value is None):
                printstr += attr+ "    "+ str(value)
            else:
                printstr += attr+ "    None"
                
            printstr += "\n"
        return printstr

    ###################################
    # when something is defined as an exon ....
    def __init__ (self):
        
        self.species_1          = 'homo_sapiens'
        self.species_2          = None
        self.exon_id_1          = None
        self.exon_id_2          = None
        self.exon_known_1       = None
        self.exon_known_2       = None
  
        self.cigar_line         = None
        self.similarity         = None
        self.bitmap             = None
        self.source             = None
        self.warning            = None

    ###################################
    def load_from_db (self, db_row, cursor, paralogue=False):
        
        if (paralogue):
            [exon_map_id, exon_id, cognate_exon_id,  
             exon_known, cognate_exon_known,  
             cigar_line, similarity, source, msa_bitmap] = db_row
            warning = None
        else:
            [exon_map_id, exon_id, cognate_exon_id,  
             exon_known, cognate_exon_known, cognate_genome_db_id, 
             cigar_line, similarity, source, msa_bitmap, warning] = db_row
            self.species_2      = genome_db_id2species (cursor,  cognate_genome_db_id)

        self.exon_id_1          = exon_id
        self.exon_id_2          = cognate_exon_id
        self.exon_known_1       = exon_known
        self.exon_known_2       = cognate_exon_known  
        self.cigar_line         = cigar_line
        self.similarity         = similarity
        self.source             = source
        self.bitmap             = msa_bitmap
        self.warning            = warning
       
#########################################
def get_maps(cursor, ensembl_db_name, exon_id, is_known, species = 'homo_sapiens', table='exon_map'):
    
    maps = []

    switch_to_db (cursor,  ensembl_db_name[species])
    qry  = "select * from %s where exon_id = %d " % (table, exon_id)
    qry += " and exon_known = %d " % is_known
    rows = search_db (cursor, qry)
    if not rows or "ERROR" in rows[0]:
        return []

    for row in rows:
        map = Map()
        paralogue = (table=='para_exon_map')
        map.load_from_db(row, cursor, paralogue)
        maps.append(map)

    return maps

#########################################
def map2exon(cursor, ensembl_db_name, map, paralogue=False):

    # this is fake exon info! to be passe to get_exon_pepseq
    exon = Exon ()
    exon.exon_id     = map.exon_id_2
    exon.is_known    = map.exon_known_2
    if map.source == 'sw_sharp':
        exon.analysis_id = -1 
        if not paralogue:  # move to the other species
            rows = switch_to_db (cursor, ensembl_db_name[map.species_2])
            if  not rows:
                exon.exon_seq_id = -1
                return exon
        else:
            qry  = "select exon_seq_id from sw_exon where exon_id = %d " % exon.exon_id 
            rows = search_db (cursor, qry)
            if not rows or not rows[0][0]:
                exon.exon_seq_id = -1
            else:
                exon.exon_seq_id = int(rows[0][0])

    elif map.source == 'usearch':
        exon.analysis_id = -2
        if not paralogue: 
            rows = switch_to_db (cursor, ensembl_db_name[map.species_2])
            if  not rows:
                exon.exon_seq_id = -1
                return exon
        else:
            qry  = "select exon_seq_id from usearch_exon where exon_id = %d " % exon.exon_id 
            rows = search_db (cursor, qry)
            if not rows or not rows[0][0]:
                exon.exon_seq_id = -1
            else:
                exon.exon_seq_id = int(rows[0][0])
    else:
        exon.analysis_id = 1
 
    return exon



#########################################
def self_maps (cursor, ensembl_db_name, human_exons):

    maps = []
    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    # this should fill in the seqs for the coding exons
    relevant_human_exons = find_relevant_exons (cursor, human_exons, human=True) 
    for he in relevant_human_exons:
        map = Map()
        map.species_1    = 'homo_sapiens'
        map.species_2    = 'homo_sapiens'                
        map.exon_id_1    = he.exon_id
        map.exon_id_2    = he.exon_id

        map.exon_known_1 = he.is_known
        map.exon_known_2 = he.is_known
        map.similarity   = 1.0

        pepseq = get_exon_pepseq (cursor, he, verbose=True)
        map.cigar_line   = cigar_line (pepseq, pepseq)
        map.source       = 'ensembl'
        maps.append(map)  
    return maps


###################################################################################################
###################################################################################################

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
def overlap (start, end, other_start, other_end):
    if ( other_end < start): 
        return False
    elif ( end < other_start):
        return False
    else:
        return True



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
def maps_evaluate (cfg, ensembl_db_name, human_exons, ortho_exons, aligned_seq, exon_positions):

    maps = []
   
    if len(aligned_seq.keys()) > 2:
        print "right now the mapping implemented for two species only"
        return []

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
                    return []

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
def make_maps (cursor, ensembl_db_name, cfg, acg, ortho_species, human_exons, ortho_exons):

    maps = []

    #print "############################## human"
    switch_to_db(cursor,  ensembl_db_name['homo_sapiens'])
    relevant_human_exons = find_relevant_exons (cursor, human_exons, human=True)
    #print "##############################", ortho_species
    switch_to_db(cursor,  ensembl_db_name[ortho_species])
    relevant_ortho_exons = find_relevant_exons (cursor, ortho_exons, human=False)

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
            = smith_waterman_context (human_seq, ortho_seq, -5, -3)

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
    maps = maps_evaluate (cfg, ensembl_db_name, relevant_human_exons, 
                          relevant_ortho_exons, aligned_seq, exon_positions)

    return maps
