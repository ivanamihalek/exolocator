#!/usr/bin/python

import MySQLdb, subprocess, re, commands
from tempfile     import NamedTemporaryFile
from bisect       import bisect
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

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
    # if not rows: return
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
    #if not rows: return
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
		   end_in_gene,  strand,  dnaseq, left_flank, right_flank, pepseq):

    switch_to_db(cursor, db_name)

    #################################
    fixed_fields  = {}
    update_fields = {}

    fixed_fields['maps_to_human_exon_id'] = human_exon_id

    update_fields['gene_id']       = gene_id
    update_fields['start_in_gene'] = start_in_gene
    update_fields['end_in_gene']   = end_in_gene
    update_fields['strand']        = strand

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

    return True


#########################################
def search_for_exons(human_gene_list, db_info):

    # 
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

    # find db ids and common names for each species db
    all_species, ensembl_db_name = get_species (cursor)

    # for each human gene
    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])

    gene_ct = 0
    for human_gene_id in human_gene_list:
    #for human_gene_id in [370495]: # Known hit
    #for human_gene_id in [378768,412667]: #  p53,wls
	switch_to_db (cursor, ensembl_db_name['homo_sapiens'])        

	gene_ct += 1
	if (not gene_ct%100): 
		print gene_ct, "out of ", len(human_gene_list)
	print human_gene_id, gene2stable(cursor, human_gene_id), get_description (cursor, human_gene_id)

	# find all human exons we are tracking in the database
	human_exons = [e for e in gene2exon_list(cursor, human_gene_id) if e.covering_exon < 0 and e.is_canonical]
	human_exons.sort(key=lambda exon: exon.start_in_gene)

	maps = [get_maps(cursor, ensembl_db_name, he.exon_id, he.is_known) for he in human_exons]
	regions = dict()
	for i,ms in enumerate(maps):
            for m in ms:
	        gene_id = exon_id2gene_id(cursor, ensembl_db_name[m.species_2], m.exon_id_2, m.exon_known_2)
		region  = get_seq_region (cursor, ensembl_db_name[m.species_2], m.exon_id_2)
		if m.species_2 not in regions: regions[m.species_2] = dict()
		if gene_id not in regions[m.species_2]: regions[m.species_2][gene_id] = dict()
		if i not in regions[m.species_2][gene_id]: regions[m.species_2][gene_id][i] = []
		regions[m.species_2][gene_id][i].append(region)

	for species in all_species:
            if species == "homo_sapiens" or species not in regions: 
                continue # Skip if homo_sapiens or no exons for this species at all

	    # Get sorted list of other species to search for exons
	    switch_to_db (cursor, ensembl_db_name[species])
	    sorted_species_list = species_sort(cursor,all_species,species)[1:]

	    for i in xrange(1,len(human_exons)-1): # we cannot find first or last exon, so skip them
		if any((species == m.species_2 for m in maps[i])): continue # exon exists, so skip
		for gene_id, rs in regions[species].iteritems():
			positions = sorted(rs.keys())
			p = bisect(positions, i)
			if p == 0 or p == len(positions): 
				continue # No preceding or following exon for this gene_id, so skip
			# Get preceding and following exon regions for this gene_id
			prev_seq_regions = rs[positions[p-1]]
			next_seq_regions = rs[positions[p]]
			# Get search name, filename and strand
			searchname, searchfile, searchstrand = \
			    prev_seq_regions[0].name, prev_seq_regions[0].filename, prev_seq_regions[0].strand
			# Get search start and end points
			if searchstrand == 1:
				searchstart = max((sr.end for sr in prev_seq_regions)) + 1
				searchend = min((sr.start for sr in next_seq_regions)) - 1
			else:
				searchstart = max((sr.end for sr in next_seq_regions)) + 1
				searchend = min((sr.start for sr in prev_seq_regions)) - 1
			if searchend - searchstart < 10: continue # Search region too short, so skip
			# Get template sequences
			for species2 in sorted_species_list:
				nearestseqs  = [get_exon_seqs(cursor, m.exon_id_2, 
							      m.exon_known_2, ensembl_db_name[species2])
						for m in maps[i] if species2 == m.species_2]
				if nearestseqs:
					nearestspecies = species2
					break
			if not nearestseqs: 
				break # No template sequences, so skip (should never happen, due to homo_sapiens)
			searchtmp = NamedTemporaryFile()
			querytmp = NamedTemporaryFile()
			# Write out sequences to search
			fastacmd = acg.generate_fastacmd_gene_command(species, searchname, 
								      searchfile, searchstrand, searchstart, searchend)
			#fasta = subprocess.check_output(fastacmd, shell=True)
			fasta  = commands.getoutput(fastacmd)
			searchtmp.write(fasta)
			searchtmp.flush()
			# Write out query sequences
			for seq in nearestseqs:
				querytmp.write(">{0}\n{1}\n".format(seq[0],seq[6]))
			querytmp.flush()
			# Perform SW# search
			swsharpcmd = acg.generate_SW_nt(querytmp.name, searchtmp.name)
			#resultstr = subprocess.check_output(swsharpcmd, shell=True)
			resultstr = commands.getoutput (swsharpcmd)
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
				if identity < 0.8 * seqlen or identity < 10: continue

				# FOUND AN EXON!
				search_start, search_end, template_start, template_end = map(int,re.split('\D+',r[6])[1:5])

				# Lengthen find to align to coding frame in template
				start_adjust = (template_start-1) % 3
				search_start -= start_adjust + 1
				template_start -= start_adjust + 1
				end_adjust = -template_end % 3
				search_end -= end_adjust
				template_end -= end_adjust

				template_name = int(r[2].split('>')[1])
				template_searchseq = [s[6] for s in nearestseqs if s[0] == template_name][0]
				template_dnaseq = template_searchseq[template_start:template_end]
				template_pepseq = str(Seq(template_dnaseq,generic_dna).translate())

				searchseq = "".join(fasta.splitlines()[1:])
				dnaseq = searchseq[search_start:search_end]
				pepseq = str(Seq(dnaseq,generic_dna).translate())
				flanklen = 60
				left_flank = searchseq[search_start-flanklen if search_start>flanklen else 0 : search_start]
				right_flank = searchseq[search_end : search_end+flanklen 
							if search_end+flanklen<len(searchseq) else len(searchseq)]
				has_NNN = 1 if "NNN" in searchseq else 0
				has_stop = 1 if "*" in pepseq else 0

				prev_end_in_gene = max((get_exon(cursor,m.exon_id_2,db_name=ensembl_db_name[species]).end_in_gene
							for m in maps[positions[p-1]] if m.species_2 == species))
				start_in_gene = prev_end_in_gene + search_start + 1 # Do I need this last +1 ?
				end_in_gene   = prev_end_in_gene + search_end

				# Print out some debugging info
				print
				print "found sequence for {0}:".format(species)
				print dnaseq
				print "left flank: " + left_flank 
				print "right flank: " + right_flank
				print "translation: " + pepseq + (" (Stop codon)" if has_stop == 1 else "")
				print 
				print "based on {0}, exon {1}, position {2} - {3}:".format(nearestspecies,
										       template_name,template_start,template_end)
				print template_dnaseq
				print "translation: " + template_pepseq
				print
				print resultstr
				print prev_end_in_gene
				print
				print "storing to ", ensembl_db_name[species]
				#raw_input("Press ENTER to continue")
				strand = -1 if searchstrand != 1 else 1
				store_sw_exon (cursor,ensembl_db_name[species], human_exons[i].exon_id, gene_id, 
					       start_in_gene, end_in_gene, strand, dnaseq, left_flank, right_flank, pepseq)

#########################################
def main():
    
    no_threads = 15

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

    parallelize (no_threads, search_for_exons, gene_list[5000:], [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
	main()


