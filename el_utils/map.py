
import os
from ensembl import genome_db_id2species
from mysql   import switch_to_db, search_db
from exon    import Exon

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
    relevant_human_exons = find_relevant_exons (cursor, human_exons) 
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
def maps_evaluate (cfg, human_exons, ortho_exons, aligned_seq, exon_positions):

    maps = []
   
    if len(aligned_seq.keys()) > 2:
        print "right now the mapping implemented for two species only"
        return maps

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

                exon_seq_human   = aligned_seq['homo_sapiens'][human_start:human_end].replace('#','-')
                exon_seq_other   = aligned_seq[other_species][other_start:other_end].replace('#','-')
                [seq_human, seq_other] = pad_the_alnmt (exon_seq_human,human_start,
                                                        exon_seq_other, other_start)

                seq = {'human':seq_human, 'other':seq_other}
                seq = strip_gaps(seq)
                if not seq:  
                    c=inspect.currentframe()
                    print " in %s:%d" % ( c.f_code.co_filename, c.f_lineno)
                    continue

                map.similarity = pairwise_tanimoto(seq['human'], seq['other'])
                if map.source == 'usearch':
                    print "\t", other_species,  map.similarity, min_similarity

                if map.similarity < min_similarity: continue

                ciggy = cigar_line (seq['human'], seq['other'])
                map.cigar_line = ciggy
                                                   
                maps.append(map)                

    return maps


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
    maps = maps_evaluate (cfg, relevant_human_exons, relevant_ortho_exons, aligned_seq, exon_positions)

    return maps
