#!/usr/bin/python

import MySQLdb, subprocess, re
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
def main():
	acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
	db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
	cursor = db.cursor()

	# find db ids and common names for each species db
	all_species, ensembl_db_name = get_species (cursor)

	# for each human gene
	switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
	human_gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
	for human_gene_id in human_gene_ids:
	#for human_gene_id in [378768,412667]: #  p53,wls
		print human_gene_id
		# find all human exons we are tracking in the database
		switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
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
			if species == "homo_sapiens" or species not in regions: continue # Skip if homo_sapiens or no exons for this species at all
			
			# Get sorted list of other species to search for exons
			switch_to_db (cursor, ensembl_db_name[species])
			sorted_species_list = species_sort(cursor,all_species,species)[1:]
			
			for i in xrange(1,len(human_exons)-1): # we cannot find first or last exon, so skip them
				if any((species == m.species_2 for m in maps[i])): continue # exon exists, so skip
				for gene_id, rs in regions[species].iteritems():
					positions = sorted(rs.keys())
					p = bisect(positions, i)
					if p == 0 or p == len(positions): continue # No preceding or following exon for this gene_id, so skip
					# Get preceding and following exon regions for this gene_id
					prev_seq_regions = rs[positions[p-1]]
					next_seq_regions = rs[positions[p]]
					# Get search name, filename and strand
					searchname, searchfile, searchstrand = prev_seq_regions[0].name, prev_seq_regions[0].filename, prev_seq_regions[0].strand
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
						nearestseqs  = [get_exon_seqs(cursor, m.exon_id_2, m.exon_known_2, ensembl_db_name[species2])
										for m in maps[i] if species2 == m.species_2]
						if nearestseqs:
							nearestspecies = species2
							break
					if not nearestseqs: break # No template sequences, so skip (should never happen, due to homo_sapiens)
					with NamedTemporaryFile() as searchtmp, NamedTemporaryFile() as querytmp:
						# Write out sequences to search
						fastacmd = acg.generate_fastacmd_gene_command(species, searchname, searchfile, searchstrand, searchstart, searchend)
						fasta = subprocess.check_output(fastacmd, shell=True)
						searchtmp.write(fasta)
						searchtmp.flush()
						# Write out query sequences
						for seq in nearestseqs:
							querytmp.write(">{}\n{}\n".format(seq[0],seq[6]))
						querytmp.flush()
						# Perform SW# search
						swsharpcmd = acg.generate_SW_nt(querytmp.name, searchtmp.name)
						resultstr = subprocess.check_output(swsharpcmd, shell=True)
					for r in (f.splitlines() for f in resultstr.split("#"*80+"\n")):
						if len(r) < 14: continue # Skip blank or malformed results

						# Parse result
						seqlen = min(int(re.split('\D+',r[1])[1]),int(re.split('\D+',r[3])[1]))
						s1start, s1end, s2start, s2end = map(int,re.split('\D+',r[6])[1:5])
						identity, matchlen = map(int,re.split('\D+',r[7])[1:3])

						#similarity = int(re.split('\D+',r[8])[1])
						#gaps = int(re.split('\D+',r[9])[1])
						#score = float(r[10].split()[1])
						#seq2 = "".join((s.split()[1] for s in r[13::3]))
						#assert(matchlen == int(re.split('\D+',r[8])[2]) and matchlen == int(re.split('\D+',r[9])[2]))

						# Reject if identity too low or too short
						if identity < 0.8 * seqlen or identity < 10: continue

						# FOUND AN EXON!
						searchseq = "".join(fasta.splitlines()[1:])
						dnaseq = "".join((s.split()[1] for s in r[12::3]))
						assert (dnaseq == searchseq[s1start-1:s1end])
						left_flank = searchseq[s1start-16 if s1start>16 else 0: s1start]
						right_flank
						print resultstr # DEBUG
						prev_end_in_gene = max((get_exon(cursor,m.exon_id_2,db_name=ensembl_db_name[species]) 
												for m in maps[positions[p-1]] if m.species_2 == species))
						start_in_gene = prev_end_in_gene + s1start
						end_in_gene   = prev_end_in_gene + s1end
						
						switch_to_db(cursor,ensembl_db_name[species])
						qry  = "insert into sw_exon (gene_id, start_in_gene, end_in_gene, maps_to_human_exon_id, "
						qry += "exon_seq_id, strand, phase, has_NNN, has_stop, has_3p_ss, has_5p_ss) values "
						qry += "('{}',{},{},{}".format(gene_id, start_in_gene, end_in_gene, human_exons[i])
						qry += ",{},{},{}".format('NULL' ,-1 if searchstrand != 1 else 1, 'NULL')
						qry += ",{},{},{},{})".format()

						qry  = "insert into exon_seq (exon_id, is_known, is_sw, dna_seq, left_flank, "
						qry += "right_flank, protein_seq, pepseq_transl_start, pepseq_transl_end) values "
						#qry += "({},0,1,'{}','{}',".format(,)
						#qry += "('{}','{}',{},{})".format()
						
						
						nearestid = r[2].split('>')[1]
						nearestseqinfo = [l for l in nearestseqs if l[0] == nearestid]
						nearestid, nearestpepseq, nearestpepseq_start, nearestpepseq_end, nearestleft_flank, nearestright_flank, nearestdna_seq = seq2info[0]
						switch_to_db(cursor, ensembl_db_name[nearestspecies])
						nearestexon = get_exon(cursor, nearest_id)
						gene_start, gene_end = get_gene_start(cursor, ensembl_db_name[species], gene_id)
						
						raw_input("Press ENTER to continue")
						
					

#########################################
if __name__ == '__main__':
	main()


