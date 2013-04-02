#!/usr/bin/python

import MySQLdb, subprocess, re, commands
import copy
from tempfile     import NamedTemporaryFile
from bisect       import bisect

from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna
from Bio.Data     import CodonTable

from difflib      import SequenceMatcher
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  Map, get_maps
from   el_utils.tree    import  species_sort
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.threads import  parallelize
from subprocess import Popen, PIPE, STDOUT

#########################################
verbose = True



#########################################
MAX_SEARCH_LENGTH = 100000
FLANK_LENGTH      = 60



#########################################
class Seq_Region:
    def __init__(self, name, filename, strand, start, end):
        self.species = filename.split('.')[0].lower()
	self.name = name
	self.filename = filename
	self.strand = 1 if strand == 1 else 2
	self.start = start
	self.end = end
    def __str__(self):
        return "(" + ",".join((self.name,self.filename,str(self.strand),str(self.start),str(self.end))) + ")"


#########################################
def get_best_filename(names_string):
    names_list = names_string.split()
    chromo_names = [n for n in names_list if 'chromosom' in n]
    if chromo_names: return chromo_names[0]
    else: return names_list[0]

#########################################
def get_gene_start_end(cursor, ensembl_db_name, gene_id):
    switch_to_db (cursor, ensembl_db_name)
    qry  = "select seq_region_start, seq_region_end from gene where gene_id=%d" % gene_id
    rows = search_db(cursor, qry)
    if not rows: 
        return None
    return  map(int, rows[0])

#########################################
def get_seq_region(cursor, ensembl_db_name, exon_id, exon_known):

    switch_to_db (cursor, ensembl_db_name)
    qry  = "select name, file_name, seq_region_strand, seq_region_start, seq_region_end "
    if exon_known == 1:
        qry += "from exon "
    elif exon_known == 0:
        qry += "from prediction_exon "
    elif exon_known == 2:
        print "get_seq_region() from sw exon not implemented - need to finish it up from the gene region "
        return None
    else:
        print "error in get_seq_region() "
        return None
        
    qry += "join seq_region on exon.seq_region_id = seq_region.seq_region_id "
    qry += "where exon_id=%d" % exon_id
    rows = search_db(cursor, qry)
    if not rows:   return False
    if not len(rows[0]) ==5:
        print rows[0]
        return False
    name, file_names, strand, start, end = rows[0]
    return Seq_Region(name, get_best_filename(file_names), int(strand), int(start), int(end))

#########################################
def get_alt_seq_region(cursor,ensembl_db_name, exon_id, exon_known):

    switch_to_db (cursor, ensembl_db_name)
    qry  = "select name, file_name, seq_region_strand, exc_seq_region_start, exc_seq_region_end "
    if exon_known == 1:
        qry += "from exon "
    elif exon_known == 0:
        qry += "from prediction_exon "
    elif exon_known == 2:
        print "get_alt_seq_region() from sw exon not implemented - need to fishi it up from the gene region "
        return None
    else:
        print "error in get_alt_seq_region() "
        return None
    
    qry += "join assembly_exception as ae on ae.seq_region_id=exon.seq_region_id "
    qry += "join seq_region on seq_region.seq_region_id=ae.exc_seq_region_id "
    qry += "where exon_id=%d" % exon_id
    rows = search_db(cursor, qry)
    if not rows: return False
    name, file_names, strand, start, end = rows[0]
    return Seq_Region(name, get_best_filename(file_names), int(strand), int(start), int(end))

#########################################
def find_sw_exon_id (cursor, human_exon_id):

    query  = "select exon_id  from sw_exon "
    query += " where maps_to_human_exon_id = '%s' " % human_exon_id
    rows   = search_db(cursor, query)
    if not rows or  'Error' in rows[0]:
        return ''
    else:
        return int(rows[0][0])

#########################################
def find_exon_seq_id (cursor, sw_exon_id):

    query  = "select exon_seq_id  from exon_seq "
    query += " where exon_id = %d  and is_sw = 1 " % sw_exon_id
    rows   = search_db(cursor, query)
    if not rows or  'Error' in rows[0]:
        return ''
    else:
        return int (rows[0][0])

#########################################
def store_sw_exon (cursor, db_name, human_exon_id, gene_id, start_in_gene, 
		   end_in_gene,  strand,  dnaseq, left_flank, right_flank, pepseq,
                   has_NNN, has_stop, template_exon_id, template_species):

    switch_to_db(cursor, db_name)

    #################################
    fixed_fields  = {}
    update_fields = {}

    fixed_fields['maps_to_human_exon_id'] = human_exon_id

    update_fields['gene_id']          = gene_id
    update_fields['start_in_gene']    = start_in_gene
    update_fields['end_in_gene']      = end_in_gene
    update_fields['strand']           = strand
    update_fields['has_NNN']          = has_NNN
    update_fields['has_stop']         = has_stop

    update_fields['template_exon_id'] = template_exon_id
    update_fields['template_species'] = template_species

    store_or_update (cursor, 'sw_exon', fixed_fields, update_fields)

    #################################
    sw_exon_id = find_sw_exon_id (cursor, human_exon_id)
    if not sw_exon_id: return False

    #################################
    fixed_fields  = {}
    update_fields = {}

    fixed_fields ['exon_id']             = sw_exon_id
    fixed_fields ['is_sw']               = 1

    update_fields['is_known']            = 0
    update_fields['dna_seq']             = dnaseq
    update_fields['left_flank']          = left_flank
    update_fields['right_flank']         = right_flank
    update_fields['protein_seq']         = pepseq
    update_fields['pepseq_transl_start'] = 0
    update_fields['pepseq_transl_end']   = len(dnaseq)

    store_or_update (cursor, 'exon_seq', fixed_fields, update_fields)

    #################################
    exon_seq_id = find_exon_seq_id (cursor, sw_exon_id)
    if not exon_seq_id: return False

    #################################
    qry  = "update sw_exon set exon_seq_id = {0} where exon_id = {1}".format(exon_seq_id, sw_exon_id)
    cursor.execute(qry)

    return sw_exon_id


#########################################
def translate(searchseq, search_start, search_end, mitochondrial):
    dnaseq = searchseq[search_start:search_end]

    if mitochondrial:
        try:
            pepseq = str(Seq(dnaseq,generic_dna).translate(table="Vertebrate Mitochondrial"))
        except:
            pepseq = ""
    else:
        try:
            pepseq = str(Seq(dnaseq,generic_dna).translate())
        except:
            print dnaseq
            print Seq(dnaseq,generic_dna).translate()
            pepseq = ""
  
    return pepseq


#########################################
def find_the_most_similar_frame (cursor, mitochondrial, template_pepseq, searchseq,  
                                 match_start, match_end, patched_target_seq):


    patched_searchseq  = searchseq[:match_start]    

    patched_searchseq += patched_target_seq
    patched_searchseq += searchseq[match_end:]


    searchseq = patched_searchseq

    pepseq = translate(searchseq, match_start, match_end, mitochondrial)

    if pepseq==template_pepseq:
        return  [pepseq, match_start, match_end]

    s = SequenceMatcher(None, pepseq, template_pepseq)
    if s.ratio() > 0.5 and not '*' in pepseq:
        return  [pepseq, match_start, match_end]

    ratio = {}
    for b in range(-2, 3): # begin offset
        if match_start+b < 0: continue
        if match_start+b > len(searchseq): continue
        for e in range(-2, 3): # end offset
            if match_end + e < 0: continue
            if match_end + e > len(searchseq): continue
            
            pepseq = translate(searchseq, match_start+b, match_end+e,  mitochondrial)
            s = SequenceMatcher(None, pepseq, template_pepseq)
            key = (b,e)
            ratio[key] =  s.ratio()
            
    
    sorted_keys = sorted(ratio, key=lambda key: -ratio[key])
    
    key = sorted_keys[0]
    (b,e) = key
    if ratio[key] < 0.5:
        return ["", match_start, match_end]

    pepseq = translate(searchseq, match_start+b, match_end+e,  mitochondrial)

    if "*" in pepseq and len(pepseq) > len(template_pepseq): # if the stop is outside of the matching region with template, 
                                                             # return only the matching region
        print 'resolving a stop codon'
        s      = SequenceMatcher(None, pepseq, template_pepseq)
        blocks = s.get_matching_blocks()
        blocks = filter( lambda (i,j,n): n>0, blocks)
        if len(blocks)==1:
            (i, j, n) = blocks[0]
            if n == len(template_pepseq):
                b += 3*i
                e = match_start + b + 3*n - match_end  
                pepseq = pepseq[i:i+n]
        else:
            (i_b,j_b,n_b) = blocks[0]
            (i_e,j_e,n_e) = blocks[-1]
            matching_region_length = i_e+n_e-i_b
            if matching_region_length == len(template_pepseq):
                b += 3*i_b
                e = match_start + b + 3*(i_e+n_e-i_b)  - match_end
                pepseq = pepseq[i_b:i_e+n_e]

    if  "*" in pepseq:
        pepseq = "" # failure

    return  [pepseq, match_start+b, match_end+e]
        

#########################################
def  check_sw_exon_exists(cursor, ensembl_db_name, species, maps_to_human_exon_id):
    switch_to_db (cursor, ensembl_db_name[species])
    qry  = "select count(1) from sw_exon where maps_to_human_exon_id=%d" % maps_to_human_exon_id
    rows = search_db(cursor, qry)
    if not rows or 'ERROR' in rows[0] or rows[0][0] == 0:
        return False
    else:
        return True
    
#########################################
def get_fasta (acg, species, searchname, searchfile, searchstrand, searchstart, searchend):

    fasta    = None
    fastacmd = acg.generate_fastacmd_gene_command(species, searchname, searchfile, searchstrand, searchstart, searchend)
    p      = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
    fasta, errmsg = p.communicate()
    if errmsg and 'From location cannot be greater than' in errmsg:
        match = re.search('\d+', errmsg )
        max_searchend = match.group(0)
        if max_searchend < searchstart: return fasta

        fastacmd = acg.generate_fastacmd_gene_command(species, searchname, searchfile, searchstrand, searchstart, max_searchend)
        p      = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
        fasta, errmsg = p.communicate()


    return fasta

#########################################
def sw_search(acg, species,  prev_seq_region, next_seq_region, template_seq):

    resulststr  = None
    searchstr   = None
    searchstart = None
    ###########################################################
    # determine the search region, adn extract it using fastascmd

    if not next_seq_region and not prev_seq_region:
        print "no regions specified in sw_searc()"
        return [resulststr, searchstr, searchstart]
    
    if not prev_seq_region:
        searchname   = next_seq_region.name
        searchfile   = next_seq_region.filename
        searchstrand = next_seq_region.strand
        if searchstrand==1:
            searchstart = next_seq_region.start - MAX_SEARCH_LENGTH
            searchend   = next_seq_region.start
        else:
            searchstart = next_seq_region.end
            searchend   = next_seq_region.end   + MAX_SEARCH_LENGTH

    elif not next_seq_region:
        searchname   = prev_seq_region.name
        searchfile   = prev_seq_region.filename
        searchstrand = prev_seq_region.strand

        if searchstrand==1:
            searchstart = prev_seq_region.end
            searchend   = prev_seq_region.end   + MAX_SEARCH_LENGTH
        else:
            searchstart = prev_seq_region.start - MAX_SEARCH_LENGTH
            searchend   = prev_seq_region.start

    else:
        if prev_seq_region.name == next_seq_region.name:
            searchname =  prev_seq_region.name
        else:
            print "prev_seq_region.name != next_seq_region.name ",  
            print prev_seq_region.name,  next_seq_region.name
            print "(cannot handle such cases yet)"
            return [resulststr, searchstr, searchstart]

        if prev_seq_region.filename == next_seq_region.filename:
            searchfile =  prev_seq_region.filename
        else:
            print "prev_seq_region.filename != next_seq_region.filename ",  
            print prev_seq_region.filename,  next_seq_region.filename
            print "(cannot handle such cases yet)"
            return [resulststr, searchstr, searchstart]

        if prev_seq_region.strand == next_seq_region.strand:
            searchstrand = prev_seq_region.strand
        else:
            print "prev_seq_region.strand != next_seq_region.strand ",  
            print prev_seq_region.strand,  next_seq_region.strand
            print "(cannot handle such cases yet)"
            return [resulststr, searchstr, searchstart]

        searchstart = prev_seq_region.end
        searchend   = next_seq_region.start

        if searchstart > searchend: # we are on the other strand
            searchstart = next_seq_region.end
            searchend   = prev_seq_region.start
    
    if searchstart < 1: searchstart=1

    searchtmp = NamedTemporaryFile(delete=True)
    querytmp  = NamedTemporaryFile(delete=True)

    ###########################################################
    # extract search region  using fastacmd
    fasta = get_fasta (acg, species, searchname, searchfile, searchstrand, searchstart, searchend)
    if not fasta: return [resulststr, searchstr, searchstart]
    searchtmp.write(fasta)
    searchtmp.flush()

    ###########################################################
    # Write out query sequences(? how many of them?)
    querytmp.write(">{0}\n{1}\n".format("query", template_seq)) # dna_seq
    querytmp.flush()

    ###########################################################
    # do  SW# search
    swsharpcmd = acg.generate_SW_nt(querytmp.name, searchtmp.name)
    #print swsharpcmd
   
    resultstr  = commands.getoutput (swsharpcmd)
    searchtmp.close()
    querytmp.close()

    #print fasta
    #print resultstr
    #print 
    #print

    ###########################################################
    # give me that sequence,  now that you have it    
    searchseq = "".join(fasta.splitlines()[1:])

    return [resultstr, searchseq, searchstart]

    
#########################################
def get_template (cursor, ensembl_db_name, map_table, species, he):

    template_species = None
    template_seq     = None

    nearest_species = species_sort(cursor, map_table.keys(), species)[1:]
    # I have a problem with the lamprey - it is an outlayer to everything else
    if species=='petromyzon_marinus':
        nearest_species.reverse()

    exon = Exon()
    len_human_protein_seq = 1.0*len(he.pepseq)
    for nearest in nearest_species:
        if not map_table[nearest][he]: continue

        m = map_table[nearest][he]

        template_seqs = get_exon_seqs (cursor,  m.exon_id_2,  m.exon_known_2,  ensembl_db_name[nearest])
        if not template_seqs:
            template_species = None
        else:
            [exon_seq_id, protein_seq, pepseq_transl_start, 
             pepseq_transl_end, left_flank, right_flank, dna_seq] = template_seqs
            if len(protein_seq)/len_human_protein_seq< 0.3: continue
            template_species    = nearest
            template_exon_id    = m.exon_id_2
            template_exon_known = m.exon_known_2
            break


    return [template_species, template_exon_id, template_exon_known, dna_seq, protein_seq]



#########################################
# he stands for "human exon"
def get_neighboring_region  (cursor, ensembl_db_name,  map_table, species, he, nbr_he):
    
    nbr_region = None

    if not nbr_he:
        # the bound will be my present left boud 
        #  + default range of search on the left
        pass
    else:
        nbr_map        = map_table[species][nbr_he]
        if not nbr_map: return None
        nbr_exon_id    = nbr_map.exon_id_2
        nbr_exon_known = nbr_map.exon_known_2

        nbr_region = get_seq_region(cursor, ensembl_db_name[species], nbr_exon_id, nbr_exon_known)

    return nbr_region

#########################################
def patch_aligned_seq (aligned_target_seq, aligned_template_seq,  template_start, mitochondrial):
    
    patched_target_seq = aligned_target_seq
    patched_positions  = []
    patch_failure      = True

    # if there are no gaps, nothing to be done here
    if not '-' in patched_target_seq:
        patch_failure = False
        return [patch_failure, patched_target_seq, patched_positions]

    # only the simplest of simple cases, because 
    # 1) it is the easiest to script [otherwise I get into wuagmire of guessing whcih things might be patchable
    # 2) it hase the best chance of really being a stutter in sequencing
    # thus require no gaps in template, and only isolated gaps in target
    if '-' in aligned_template_seq:
        patch_failure = True
        return [patch_failure, patched_target_seq, patched_positions]
    if '--' in aligned_target_seq:
        patch_failure = True
        return [patch_failure, patched_target_seq, patched_positions]

    # patch if needed and possible (can it be that something is mitochondrial in one species but not in the other)
    if mitochondrial:
        codon_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
    else:
        codon_table = CodonTable.unambiguous_dna_by_name["Standard"]

    # find isolated gaps
    # do we have phase in the template seq by any chance?
    prev_pos        = (3-template_start)%3
    next_pos        = prev_pos+3
    template_length = len(aligned_template_seq)

    patched_target_seq = ""
    for pos in range(template_length)[next_pos::3]:

        # find what the corresponding codon is in template
        target_codon   = aligned_target_seq[prev_pos:pos]
        template_codon = aligned_template_seq[prev_pos:pos]
        
        if '-' in target_codon: 
            # if the two codons match, copy the third
            test_codon = ""
            for i  in range(3):
                if not target_codon[i]=='-':
                    test_codon       += target_codon[i]
                    patched_position  = prev_pos + i
                else:
                    test_codon += template_codon[i]
            # if we end up with the stop codon rigth there,  abandon the attempt
            if test_codon in codon_table.stop_codons:
                continue

            patched_target_seq += test_codon
            patched_positions.append(patched_position)
            
        else:
            patched_target_seq += target_codon
            
        prev_pos = pos

    # return warning if the patching was done
    patch_failure = False
    return [patch_failure, patched_target_seq, patched_positions]


#########################################
def parse_sw_output (resultstr):
		    
    best_match = None
    longest    = -1

    for r in (f.splitlines() for f in resultstr.split("#"*80+"\n")):
        
        if len(r) < 14: continue # Skip blank or malformed results

        # Parse result
        seqlen = min(int(re.split('\D+',r[1])[1]),int(re.split('\D+',r[3])[1]))
        identity, matchlen = map(int, re.split('\D+', r[7])[1:3])
        #similarity = int(re.split('\D+',r[8])[1])
        #gaps       = int(re.split('\D+',r[9])[1])
        #score      = float(r[10].split()[1])

        
        # Reject if identity too low or too short -- seqlen is the length of the query
        if identity < 0.4*seqlen or identity < 10: continue

        # FOUND AN EXON!
        # ... but lets keep what might be the best match
        if matchlen > longest:
            longest = matchlen
            [search_start, search_end, template_start, template_end] = map(int,re.split('\D+',r[6])[1:5])

            aligned_qry_seq = ""
            for row in r[13::3]: # every third row
                seq = re.split('\s+',row)[2]
                aligned_qry_seq += seq
            aligned_target_seq = ""
            for row in r[12::3]: # every third row
                seq = re.split('\s+',row)[2]
                aligned_target_seq += seq
            best_match = [search_start-1, search_end-1, template_start-1, template_end-1, 
                          aligned_target_seq, aligned_qry_seq]

    return best_match 

#########################################
def translation_check ( searchseq, match_start, match_end, mitochondrial, pepseq):
    if not pepseq == translate(searchseq, match_start, match_end, mitochondrial):
        print " ! "
        print pepseq
        print translate(searchseq, match_start, match_end, mitochondrial)
        return False
    return True


#########################################
def organize_and_store_sw_exon (cursor, ensembl_db_name,  species, gene_id, 
                                gene_coords, search_start, searchseq, he,
                                match_start, match_end, mitochondrial, pepseq,
                                template_species, template_exon_id, 
                                template_start, template_end, verbose=False):

    
    [gene_seq_id, gene_start, gene_end, gene_strand] = gene_coords

    dnaseq = searchseq [match_start:match_end]
    search_end = search_start + len(searchseq)-1

    # check one more time that the translation that we are storing matches:
    if not translation_check (searchseq, match_start, match_end, mitochondrial, pepseq): return None
                
    # flanking seqeunces
    flanklen    = FLANK_LENGTH
    left_flank  = searchseq[match_start-flanklen if match_start>flanklen else 0 : match_start]
    right_flank = searchseq[match_end : match_end+flanklen if match_end+flanklen<len(searchseq) else len(searchseq)]
    # some warning signs
    has_NNN  = 1 if "NNN" in searchseq else 0
    has_stop = 1 if "*"   in pepseq    else 0

    # figure out where the start in gene is
    if gene_strand==1:
        start_in_gene = search_start + match_start - gene_start
        end_in_gene   = search_start + match_end   - gene_start - 1
    else:
        start_in_gene = gene_end - (search_end - match_start)
        end_in_gene   = gene_end - (search_end - match_end + 1)
        

    gene_strand = -1 if gene_strand != 1 else 1
    sw_exon_id  = None
    template_stable_exon_id = 
    sw_exon_id  = store_sw_exon (cursor, ensembl_db_name[species], he.exon_id, gene_id, 
                                 start_in_gene, end_in_gene, gene_strand, dnaseq, left_flank, right_flank, pepseq,
                                 has_NNN, has_stop, template_stable_exon_id, template_species)
    
    if not sw_exon_id: return None

    # the reutn is used for diagnostic purposes
    return [left_flank, right_flank, has_stop, sw_exon_id]


#########################################
def search_and_store (cursor, ensembl_db_name, acg, human_exon, species, gene_id, gene_coords, prev_seq_region, next_seq_region, 
                      template_info, mitochondrial):

    matching_region = None

    [template_species, template_exon_id, template_dna, template_pepseq] = template_info

    # do sw search
    [resultstr, searchseq, search_start] = sw_search (acg, species, prev_seq_region, next_seq_region, template_dna)
    if not resultstr: return matching_region

    # parse the output
    match = parse_sw_output(resultstr)
    if not match: return matching_region
    [match_start, match_end, template_start, template_end, 
     aligned_target_seq, aligned_template_seq] = match

    # patch isolated gaps == keep track of the positions that are patched
    [patch_failure, patched_target_seq, patched_positions] = patch_aligned_seq (aligned_target_seq, 
                                                                                aligned_template_seq, 
                                                                                template_start, mitochondrial)
    if patch_failure: return matching_region

    # how different is translation for the translated template?
    # align in all three frames and pick one which is the most similar to the template 
    # check for the length of the peptide # should I perhaps check if one of the
    # solutions that are usb-optimal in length give lnger peptide?
    # something I chould defineitley come back to one fine day
    [pepseq, match_start, match_end] = find_the_most_similar_frame (cursor, mitochondrial,
                                                                    template_pepseq, searchseq, 
                                                                    match_start, match_end,
                                                                    patched_target_seq)
    if len(pepseq)<3: return matching_region

    matching_region = copy.copy(next_seq_region)
    if prev_seq_region.strand == 1:
        matching_region.start = prev_seq_region.end + match_start
        matching_region.end   = prev_seq_region.end + match_end
    else:
        matching_region.end   = next_seq_region.start   - match_start
        matching_region.start = next_seq_region.start   - match_end

    
    ret = []
    ret = organize_and_store_sw_exon(cursor, ensembl_db_name, species, gene_id,  
                               gene_coords, search_start, searchseq, human_exon,
                               match_start, match_end, mitochondrial, pepseq,
                               template_species, template_exon_id, 
                               template_start, template_end, verbose=True)


    if verbose:# Print out some debugging info
        if ret:
            [left_flank, right_flank, has_stop, sw_exon_id] = ret
        print "============================================"
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
	human_stable      = gene2stable     (cursor, human_exon.gene_id)
        human_description = get_description (cursor, human_exon.gene_id)
        print human_exon.gene_id, human_stable, human_description 
        print "found sequence for {0}:".format(species)
        #print dnaseq
        if  ret:
            print "left flank:    " + left_flank 
            print "right flank:   " + right_flank
        #print "translation:   " + pepseq + (" (Stop codon)" if has_stop == 1 else "")
        print "translation:   " + pepseq 
        print "template:      " + template_pepseq
        print "human version: " + human_exon.pepseq
        print "based on {0}, exon {1}, position {2} - {3}:".format(template_species,template_exon_id,
                                                                   template_start, template_end)
        print "storing to ", ensembl_db_name[species]
        print "stored as exon " if  ret else "not stored"
        print 

    return matching_region

#########################################
def left_region (seq_region, region_length):

    new_region = copy.copy(seq_region)
    new_region.start = seq_region.start-2*region_length
    new_region.end   = seq_region.start-region_length

    return new_region

#########################################
def right_region (seq_region, region_length):

    new_region = copy.copy(seq_region)
    new_region.start = seq_region.end   +   region_length
    new_region.end   = seq_region.start + 2*region_length

    return new_region


#########################################
def find_missing_exons(human_gene_list, db_info):

    # 
    [local_db, ensembl_db_name] = db_info
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
        acg = AlignmentCommandGenerator()
    else:
        db  = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    # find db ids and common names for each species db
    all_species, ensembl_db_name = get_species (cursor)

    # for each human gene
    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])

    gene_ct = 0
    found   = 0
    sought  = 0
    human_gene_list.reverse()
    #for human_gene_id in human_gene_list:
    #for human_gene_id in [370495]: # Known hit
    #for human_gene_id in [378768]: #  p53
    #for human_gene_id in [412667]: #  wls
    #for human_gene_id in [418249]: #  TRF2
    #for human_gene_id in [374433]: # Known hit
    for human_gene_id in [397321]: # nice example, finds 23 out of 61
    #for human_gene_id in [397176]:


	switch_to_db (cursor, ensembl_db_name['homo_sapiens'])

	# Get stable id and description of this gene -- DEBUG
	human_stable      = gene2stable    (cursor, human_gene_id)
        human_description = get_description(cursor, human_gene_id)
	if verbose:  print human_gene_id, human_stable, human_description

	# progress counter 
	gene_ct += 1
	if (not gene_ct%10): 
            print "processed ",   gene_ct, " out of ", len(human_gene_list), "genes"
            print "exons found: ", found, "out of ", sought, "sought"

	# find all human exons for this gene we are tracking in the database in order
	human_exons = [e for e in gene2exon_list(cursor, human_gene_id) 
                       if e.covering_exon < 0 and e.is_canonical and e.is_known]
	human_exons.sort(key=lambda exon: exon.start_in_gene)
        for he in human_exons:
            he.stable_id = exon2stable (cursor, he.exon_id)

	# make 'table' of maps, which is either pointer to the map if it exists, or None
	map_table  = {}
        for species in all_species:
            map_table[species] = {}
            for he in human_exons:
                map_table[species][he] = None

        for he in human_exons:
            maps_for_exon =  get_maps(cursor, ensembl_db_name, he.exon_id, he.is_known) # exon data
            for m in maps_for_exon:
               if m.similarity < 0.7: continue
               map_table[m.species_2][he] = m

        # get rid of species that do not have the gene at all
        for species in all_species:
            one_exon_found = False
            for he in human_exons:
                if map_table[species][he]:
                    one_exon_found = True
                    break
            if not one_exon_found:
                del map_table[species]
               
        # get rid of exons  that appear in no other species but human (?)
        bad_he = []
        for he in human_exons:
            one_species_found = False
            for species in  map_table.keys(): 
                if species =='homo_sapiens': continue
                if map_table[species][he]:
                    one_species_found = True
                    break
            if not one_species_found:
                bad_he.append(he)
        human_exons = filter (lambda he: not he in bad_he, human_exons)

        # fill in the peptide sequence field for each exon
        for he in human_exons:
            he.pepseq =   get_exon_pepseq (cursor, he, ensembl_db_name['homo_sapiens'])

        # keep track of nearest neighbors for each human exon
        previous = {}
        next     = {}
        prev     = None
        for he in human_exons:
            previous[he]        = prev
            if prev: next[prev] = he
            prev = he
        next[he] = None

        # fill,  starting from the species that are nearest to the human
        species_sorted_from_human = species_sort(cursor,map_table.keys(),species)[1:]
        for species in species_sorted_from_human:
        #for species in ['choloepus_hoffmanni']:
            # see which exons have which neighbors
            no_left  = []
            no_right = []
            has_both_neighbors = []
            one_existing_map   = None
            for he in human_exons:
                m =  map_table[species][he]
                if m: 
                    one_existing_map = m
                    continue
                prev = previous[he]
                nxt  = next[he]
                if prev and nxt and map_table[species][prev] and map_table[species][nxt]:
                    has_both_neighbors.append(he)
                elif not prev or not map_table[species][prev]:
                    no_left.append(he)
                elif not nxt  or not map_table[species][nxt]:
                    no_right.append(he)
            
            if not one_existing_map: continue # this shouldn't happen
            if not has_both_neighbors and not no_left and not no_right: continue

            # what is the gene that we are talking about?
            exon_id  = one_existing_map.exon_id_2
            is_known = one_existing_map.exon_known_2
            gene_id  = exon_id2gene_id (cursor, ensembl_db_name[species], exon_id, is_known)
            # is it mitochondrial?
            mitochondrial = is_mitochondrial(cursor, gene_id, ensembl_db_name[species])
            # where is the gene origin (position on the sequence)
            gene_coords =  get_gene_coordinates (cursor, gene_id, ensembl_db_name[species])
            if not gene_coords: continue
            [gene_seq_id, gene_start, gene_end, gene_strand] = gene_coords

            # fill in exons that have both neighbors:
            # human exon functions as a coordinate here
            for he in has_both_neighbors:

                m = map_table[species][he]
                # get template (known exon from the nearest species)
                [template_species, template_exon_id, template_exon_known,
                 template_dna, template_pepseq] = get_template (cursor, ensembl_db_name, 
                                                                map_table, species, he)
                if not template_species: continue
                # previous_ and next_seq_region are of the type Seq_Region defined on the top of the file

                # get previous region
                prev_seq_region = get_neighboring_region (cursor, ensembl_db_name, 
                                                          map_table, species, he, previous[he])
                # get following  region
                next_seq_region = get_neighboring_region (cursor, ensembl_db_name, map_table, 
                                                          species, he, next[he])

                sought += 1
                template_info = [template_species,template_exon_id, template_dna, template_pepseq]
                matching_region = search_and_store(cursor, ensembl_db_name, acg, he, species, gene_id, gene_coords, 
                                                    prev_seq_region, next_seq_region, template_info, mitochondrial)
                if matching_region: found += 1


            # work backwards
            # use the last known region on the left as the bound
            no_left.reverse()
            next_seq_region = None
            for he in no_left:
                m =  map_table[species][he]
                # get template (known exon from the nearest species)
                [template_species, template_exon_id, template_exon_known,
                 template_dna, template_pepseq] = get_template (cursor, ensembl_db_name, 
                                                                map_table, species, he)
                if not template_species: continue

                # get following  region
                if not next_seq_region:
                    next_seq_region = get_neighboring_region (cursor, ensembl_db_name, map_table, 
                                                              species, he, next[he])
                if not next_seq_region: continue

                # otherwise it is the last thing we found
                # the previous region is eyeballed from the next on
                # the previous and the  next region frame the search region
                prev_seq_region = left_region (next_seq_region, MAX_SEARCH_LENGTH)
                sought         += 1
                template_info   = [template_species, template_exon_id, template_dna, template_pepseq]
                matching_region = search_and_store(cursor, ensembl_db_name, acg, he, 
                                                   species,  gene_id, gene_coords, prev_seq_region, next_seq_region,
                                                   template_info, mitochondrial)
                if matching_region: 
                    found += 1
                    next_seq_region = matching_region
 
            # repeat the whole procedure on the right
            prev_seq_region = None
            for he in no_right:
                m =  map_table[species][he]
                # get template (known exon from the nearest species)
                [template_species, template_exon_id, template_exon_known,
                 template_dna, template_pepseq] = get_template (cursor, ensembl_db_name, 
                                                                map_table, species, he)
                if not template_species: continue

                # get following  region
                if not prev_seq_region:
                    prev_seq_region = get_neighboring_region (cursor, ensembl_db_name, map_table, 
                                                              species, he, previous[he])
                if not prev_seq_region: continue
                # otherwise it is the last thing we found
                    

                # the following region is eyeballed from the previous 
                next_seq_region = right_region (prev_seq_region, MAX_SEARCH_LENGTH)
                sought         += 1
                template_info   = [template_species,template_exon_id, template_dna, template_pepseq]
                matching_region = search_and_store (cursor, ensembl_db_name, acg, he, 
                                                    species, gene_id, gene_coords, prev_seq_region, next_seq_region,
                                                    template_info, mitochondrial)
                #print he.pepseq
                #print template_pepseq
                if matching_region: 
                    found += 1
                    prev_seq_region = matching_region
                    

            
            
        if verbose: print "done with ",  human_gene_id, human_stable, human_description 
        if verbose: print "sought", sought, " found", found
                


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


    species       = 'homo_sapiens'
    switch_to_db (cursor, ensembl_db_name[species])
    gene_list     = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
    cursor.close()
    db.close()

    parallelize (no_threads, find_missing_exons, gene_list, [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
	main()

'''
       for species in species_sorted_from_human:
            for he in human_exons:
                m =  map_table[species][he]
                if not m:
                    print "%6s " % 'none',
                else:
                    print "%6.2f " % m.similarity,
            print
'''
