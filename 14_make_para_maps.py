#!/usr/bin/python


import StringIO
import MySQLdb, commands, re, sys
from hashlib import sha1
from random  import random, choice
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.utils   import  *
from   el_utils.threads import  parallelize
from   el_utils.map     import  *
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   alignment import * # C implementation of smith waterman

#########################################
verbose = False


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
        
        ct = 0
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

            if verbose: print "\t no map found - making new one"

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

        print species, "done"

    cursor.close()
    db.close()

    return True

#########################################
def main():
    
    no_threads = 10

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
