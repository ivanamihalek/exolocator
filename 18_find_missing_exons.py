#!/usr/bin/python

import MySQLdb, subprocess, re, commands
from tempfile     import NamedTemporaryFile
from bisect       import bisect
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna
from difflib      import SequenceMatcher
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  Map, get_maps
from   el_utils.tree    import  species_sort
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.threads import  parallelize

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
    # if not rows: return
    return map(int,rows[0])

#########################################
def get_seq_region(cursor, ensembl_db_name, exon_id):
    switch_to_db (cursor, ensembl_db_name)
    qry  = "select name, file_name, seq_region_strand, seq_region_start, seq_region_end "
    qry += "from exon join seq_region on exon.seq_region_id = seq_region.seq_region_id "
    qry += "where exon_id=%d" % exon_id
    rows = search_db(cursor, qry)
    if not rows:   return False
    name, file_names, strand, start, end = rows[0]
    return Seq_Region(name, get_best_filename(file_names), int(strand), int(start), int(end))

#########################################
def get_alt_seq_region(cursor,ensembl_db_name, exon_id):

    switch_to_db (cursor, ensembl_db_name)
    qry  = "select name, file_name, seq_region_strand, exc_seq_region_start, exc_seq_region_end "
    qry += "from exon join assembly_exception as ae on ae.seq_region_id=exon.seq_region_id "
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

    update_fields['template_exon_id'] = int(template_exon_id)
    update_fields['template_species'] = template_species

    store_or_update (cursor, 'sw_exon', fixed_fields, update_fields)

    #################################
    sw_exon_id = find_sw_exon_id (cursor, human_exon_id)
    if not sw_exon_id: return False

    #################################
    fixed_fields  = {}
    update_fields = {}

    fixed_fields['exon_id']              = sw_exon_id
    fixed_fields['is_sw']                = 1

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
        pepseq = str(Seq(dnaseq,generic_dna).translate(table="Vertebrate Mitochondrial"))
    else:
        pepseq = str(Seq(dnaseq,generic_dna).translate())
  
    return pepseq

#########################################
def find_the_most_similar_frame (cursor, mitochondrial, template_pepseq, searchseq,  search_start, search_end):

    pepseq = translate(searchseq, search_start, search_end, mitochondrial)
    if pepseq==template_pepseq:
        return  [pepseq, search_start, search_end]

    s = SequenceMatcher(None, pepseq, template_pepseq)
    if s.ratio() > 0.5 and not '*' in pepseq:
        return  [pepseq, search_start, search_end]

    ratio = {}
    for b in range(-2, 3): # begin offset
        if search_start+b < 0: continue
        if search_start+b > len(searchseq): continue
        for e in range(-2, 3): # end offset
            if search_end + e < 0: continue
            if search_end + e > len(searchseq): continue
            
            pepseq = translate(searchseq, search_start+b, search_end+e,  mitochondrial)
            s = SequenceMatcher(None, pepseq, template_pepseq)
            key = (b,e)
            ratio[key] =  s.ratio()
    
    sorted_keys = sorted(ratio, key=lambda key: -ratio[key])
    
    key = sorted_keys[0]
    (b,e) = key
    if ratio[key] < 0.5:
        return ["", search_start, search_end]

    pepseq = translate(searchseq, search_start+b, search_end+e,  mitochondrial)

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
                e = search_start + b + 3*n - search_end  
                pepseq = pepseq[i:i+n]
        else:
            (i_b,j_b,n_b) = blocks[0]
            (i_e,j_e,n_e) = blocks[-1]
            matching_region_length = i_e+n_e-i_b
            if matching_region_length == len(template_pepseq):
                b += 3*i_b
                e = search_start + b + 3*(i_e+n_e-i_b)  - search_end
                pepseq = pepseq[i_b:i_e+n_e]
            

    return  [pepseq, search_start+b, search_end+e]
        

#########################################
def search_for_exons(human_gene_list, db_info):

    verbose = False

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
    found  = 0
    sought = 0
    #human_gene_list.reverse()
    for human_gene_id in human_gene_list:
    #for human_gene_id in [370495]: # Known hit
    #for human_gene_id in [378768]: #  p53
    #for human_gene_id in [412667]: #  wls
    #for human_gene_id in [374433]: # Known hit
    #for human_gene_id in [397321]: # nice example, finds 23 out of 61
	switch_to_db (cursor, ensembl_db_name['homo_sapiens'])

	# Get stable id and description of this gene -- DEBUG
	human_stable      = gene2stable    (cursor, human_gene_id)
	human_description = get_description(cursor, human_gene_id)

	# progress counter -- DEBUG
	gene_ct += 1
	if (not gene_ct%10): 
            print "processed ",   gene_ct, " out of ", len(human_gene_list), "genes"
            print "exons found: ", found, "out of ", sought, "sought"

	# find all human exons for this gene we are tracking in the database in order
	human_exons = [e for e in gene2exon_list(cursor, human_gene_id) if e.covering_exon < 0 and e.is_canonical]
	human_exons.sort(key=lambda exon: exon.start_in_gene)

	# get orthologs for each exon
	maps = []
        for he in human_exons:
           maps_for_exon =  get_maps(cursor, ensembl_db_name, he.exon_id, he.is_known) # exon data
           # should delete later:
           maps_for_exon = filter (lambda m: not m.source == 'sw_sharp',  maps_for_exon)
           maps.append(maps_for_exon)

	regions = dict() # region data
	for i,ms in enumerate(maps):
            for m in ms:
                if  m.exon_known_2 == 2: # sw exon
                    continue  # I should be able to use those
	        gene_id = exon_id2gene_id(cursor, ensembl_db_name[m.species_2], m.exon_id_2, m.exon_known_2)
		region  = get_seq_region (cursor, ensembl_db_name[m.species_2], m.exon_id_2)
                if not region: 
                    region  = get_alt_seq_region (cursor, ensembl_db_name[m.species_2], m.exon_id_2)
                if not region:
                    continue
		if m.species_2 not in regions: 
                    regions[m.species_2] = dict()
		if gene_id     not in regions[m.species_2]: 
                    regions[m.species_2][gene_id] = dict()
		if i not in regions[m.species_2][gene_id]: 
                    regions[m.species_2][gene_id][i] = []
		regions[m.species_2][gene_id][i].append(region)

	for species in all_species:
	#for species in ['felis_catus']:
            if species == "homo_sapiens" or species not in regions: continue # Skip if homo_sapiens or no exons for this species at all

	    # Get sorted list of other species to search for exons for this species
	    switch_to_db (cursor, ensembl_db_name[species])
	    sorted_species_list = species_sort(cursor,all_species,species)[1:]

	    #for i in xrange(1,len(human_exons)-1): # we cannot find first or last exon, so skip them -- DEBUG
	    for i in xrange(len(human_exons)):
		# skip if human exon not known or exon exists
		if any((m.exon_known_1 == 0 or m.species_2 == species for m in maps[i])): continue
		
		for gene_id, rs in regions[species].iteritems():
		    # Get preceding and following exon regions for this gene_id
		    positions = sorted(rs.keys())
		    p = bisect(positions, i)
		    if p == 0: prev_seq_regions = None
		    else:      prev_seq_regions = rs[positions[p-1]]
		    if p == len(positions): next_seq_regions = None
		    else:                   next_seq_regions = rs[positions[p]]

		    # Get search name, filename and strand
		    if prev_seq_regions is not None: r = prev_seq_regions[0]
		    else:                            r = next_seq_regions[0]
		    searchname, searchfile, searchstrand = r.name, r.filename, r.strand

		    # Get search start and end points
		    endsearchlength = 100000 # how many bps to search for start and end exons
		    if searchstrand == 1:
			start_regions = prev_seq_regions
			end_regions = next_seq_regions
		    else:
			start_regions = next_seq_regions
			end_regions = prev_seq_regions
		    if start_regions is not None:
			searchstart = max((sr.end for sr in start_regions)) + 1
		    else: searchstart = None
		    if end_regions is not None:
			searchend   = min((sr.start for sr in end_regions)) - 1
		    else: searchend = searchstart + endsearchlength
		    if searchstart is None: searchstart = searchend - endsearchlength
		    if searchstart < 0: searchstart = 0

		    # Skip if search region too short
		    if searchend and searchstart and searchend - searchstart < 10: continue

		    # Get template sequences
		    for species2 in sorted_species_list:
			nearestseqs  = [get_exon_seqs(cursor, m.exon_id_2, m.exon_known_2, ensembl_db_name[species2]) + [m.exon_id_2]
					for m in maps[i] if species2 == m.species_2]
			if nearestseqs:
			    nearestspecies = species2
			    break
		    if not nearestseqs: break # No template sequences, so skip (should never happen, due to homo_sapiens)

		    searchtmp = NamedTemporaryFile()
		    querytmp  = NamedTemporaryFile()
		    # Write out sequences to search
		    fastacmd = acg.generate_fastacmd_gene_command(species, searchname, searchfile, searchstrand, searchstart, searchend)
		    #fasta = subprocess.check_output(fastacmd, shell=True)
		    fasta  = commands.getoutput(fastacmd)
		    searchtmp.write(fasta)
		    searchtmp.flush()
		    # Write out query sequences
		    for seq in nearestseqs:
			querytmp.write(">{0}\n{1}\n".format(seq[0],seq[6])) # exon_seq_id, dna_seq
		    querytmp.flush()
		    # Perform SW# search
                    sought += 1
		    swsharpcmd = acg.generate_SW_nt(querytmp.name, searchtmp.name)
		    resultstr  = commands.getoutput (swsharpcmd)
		    searchtmp.close()
		    querytmp.close()
		    for r in (f.splitlines() for f in resultstr.split("#"*80+"\n")):
			if len(r) < 14: continue # Skip blank or malformed results

			# Parse result
			seqlen = min(int(re.split('\D+',r[1])[1]),int(re.split('\D+',r[3])[1]))
			identity, matchlen = map(int,re.split('\D+',r[7])[1:3])
			#similarity = int(re.split('\D+',r[8])[1])
			#gaps = int(re.split('\D+',r[9])[1])
			#score = float(r[10].split()[1])

			# Reject if identity too low or too short
			if identity < 0.4*seqlen or identity < 10: continue

			# FOUND AN EXON!
			search_start, search_end, template_start, template_end = map(int,re.split('\D+',r[6])[1:5])

			# Lengthen find to align to coding frame in template
                        # what _is_ this?
			start_adjust    = (template_start-1) % 3
			search_start   -= start_adjust + 1
			template_start -= start_adjust + 1

			end_adjust      = -template_end % 3
			search_end     -= end_adjust
			template_end   -= end_adjust

			template_name      = int(r[2].split('>')[1])
			template_seq       = [s for s in nearestseqs if s[0] == template_name][0] 
			template_gene_id   = template_seq[7]
			template_searchseq = template_seq[6]
			template_dnaseq    = template_searchseq[template_start:template_end]

                        # why are we doing this?
			#if is_mitochondrial(cursor, template_gene_id, ensembl_db_name[nearestspecies]):
			#    template_pepseq    = str(Seq(template_dnaseq,generic_dna).translate(table="Vertebrate Mitochondrial"))
			#else:
			#    template_pepseq    = str(Seq(template_dnaseq,generic_dna).translate())
                        template_pepseq   = template_seq[1]


			searchseq = "".join(fasta.splitlines()[1:])
                        # how different is translation for the translated template?
                        # align in all three frames and pick one which is the most similar to the template 
                        mitochondrial = is_mitochondrial(cursor, gene_id, ensembl_db_name[species])
                        [pepseq, search_start, search_end] = find_the_most_similar_frame (cursor, mitochondrial,
                                                                                          template_pepseq, searchseq, 
                                                                                          search_start, search_end)
                        if len(pepseq)<3: continue
			flanklen = 60
			left_flank = searchseq[search_start-flanklen if search_start>flanklen else 0 : search_start]
			right_flank = searchseq[search_end : search_end+flanklen if search_end+flanklen<len(searchseq) else len(searchseq)]
			has_NNN  = 1 if "NNN" in searchseq else 0
			has_stop = 1 if "*" in pepseq else 0

			prev_end_in_gene = max((get_exon(cursor,m.exon_id_2,db_name=ensembl_db_name[species]).end_in_gene
						for m in maps[positions[p-1]] if m.species_2 == species))
			start_in_gene = prev_end_in_gene + search_start + 1 # Do I need this last +1 ?
			end_in_gene   = prev_end_in_gene + search_end

                        dnaseq = searchseq [search_start:search_end]
                        found += 1
                        if not pepseq == translate(searchseq, search_start, search_end, mitochondrial):
                            print " ! "
                            print pepseq
                            print translate(searchseq, search_start, search_end, mitochondrial)
                            continue


			strand = -1 if searchstrand != 1 else 1
			sw_exon_id = store_sw_exon (cursor,ensembl_db_name[species], human_exons[i].exon_id, gene_id, 
				       start_in_gene, end_in_gene, strand, dnaseq, left_flank, right_flank, pepseq,
                                       has_NNN, has_stop, template_name, nearestspecies)

			if verbose:# Print out some debugging info
                            print "============================================"
                            print human_gene_id, human_stable, human_description 
                            print "found sequence for {0}:".format(species)
                            print dnaseq
                            print "left flank:  " + left_flank 
                            print "right flank: " + right_flank
                            print "translation: " + pepseq + (" (Stop codon)" if has_stop == 1 else "")
                            print "template:    " + template_pepseq
                            print "based on {0}, exon {1}, position {2} - {3}:".format(nearestspecies,template_name,
                                                                                       template_start,template_end)
                            print "storing to ", ensembl_db_name[species]
                            print "stored as exon ", sw_exon_id
                            print "exons found: ", found, "out of ", sought, "sought"
                            print 

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


    species       = 'homo_sapiens'
    switch_to_db (cursor, ensembl_db_name[species])
    gene_list     = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
    cursor.close()
    db.close()

    parallelize (no_threads, search_for_exons, gene_list[5000:-5000], [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
	main()


