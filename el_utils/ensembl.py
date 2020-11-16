#!/usr/bin/python
from time import time

import MySQLdb
from el_utils.mysql import *
from el_utils.exon import Exon
import subprocess
import os

def get_trivial(cursor, species_names):
	trivial = {}
	for spec in species_names:
		qry = "select trivial_name from exolocator_meta.species_names where species='%s'" % spec
		ret = error_intolerant_search(cursor, qry)
		trivial[spec] = ret[0][0] if (ret and ret[0] and ret[0][0] and len(ret[0][0]) > 0) else spec
	return trivial


#########################################
def get_species_shorthand(cursor, species):
	switch_to_db(cursor, 'ensembl_meta')

	qry = "select shorthand from species_names where species='%s'" % species
	rows = search_db(cursor, qry)
	if not rows: return ""

	return rows[0][0]


#########################################
def get_current_db(cursor):
	return hard_landing_search(cursor, "SELECT DATABASE() FROM DUAL")[0][0]


#########################################
def canonical_transl_info(cursor, gene_id):
	qry = "select canonical_transcript_id from gene  where gene_id = %d " % gene_id
	rows = search_db(cursor, qry)
	if (not rows):
		search_db(cursor, qry, verbose=True)
		return []

	canonical_transcript_id = int(rows[0][0])
	qry = "select start_exon_id, seq_start, end_exon_id,  seq_end "
	qry += " from translation where transcript_id = %d " % canonical_transcript_id
	rows = search_db(cursor, qry)
	if (not rows):
		search_db(cursor, qry, verbose=True)
		return []

	return rows[0]


#########################################
def check_ccds(cursor, transcript_stable_id="", transcript_id=""):
	ccds = ""

	qry = "select dna_align_feature.hit_name "
	qry += "from dna_align_feature, transcript, transcript_supporting_feature "
	qry += "   where dna_align_feature.dna_align_feature_id =  transcript_supporting_feature.feature_id "
	qry += "   and transcript_supporting_feature.feature_type ='dna_align_feature' "
	qry += "   and transcript_supporting_feature.transcript_id =transcript.transcript_id "
	if (transcript_stable_id):
		qry += "   and transcript.stable_id = '%s' " % transcript_stable_id
	elif (transcript_id):
		qry += "   and transcript.transcript_id = '%s' " % transcript_id
	else:
		return ccds

	rows = search_db(cursor, qry)

	if not rows:
		return ccds

	for row in rows:
		if 'CCDS' in row[0]:
			ccds = row[0]

	return ccds


#########################################
def get_exons(cursor, species_db, gene_id, table):
	if table == 'exon':
		return get_ensembl_exons(cursor, species_db, gene_id)
	elif table == 'exolocator_exon':
		return get_novel_exons(cursor, gene_id, table)

	return []


#########################################
def get_pepseq_transl_range(cursor, exon_id, exon_known, species_db):
	qry = f"select pepseq_transl_start, pepseq_transl_end from {species_db}.exon_seq "
	qry += "where exon_id = %d " % int(exon_id)
	qry += "and is_known = %d " % int(exon_known)
	rows = error_intolerant_search(cursor, qry)
	return rows[0] if rows else []


#########################################
def get_canonical_exon_ids(cursor, canonical_transcript_id, species_db):
	qry = f"select exon_id from {species_db}.exon_transcript "
	qry += " where transcript_id = %d " % canonical_transcript_id
	rows = error_intolerant_search(cursor, qry)
	return [row[0] for row in rows] if rows else []


#########################################
def get_canonical_coordinates(cursor, canonical_transcript_id, species_db):
	qry  = "select seq_start, start_exon_id,  seq_end, end_exon_id "
	qry += f"from {species_db}.translation where transcript_id = {canonical_transcript_id} "
	rows = error_intolerant_search(cursor, qry)
	return rows[0] if rows else []


#########################################
def get_canonical_transcript_id(cursor, gene_id, species_db):
	qry  = "select canonical_transcript_id "
	qry += f"from  {species_db}.gene where gene_id={gene_id}"
	rows = error_intolerant_search(cursor, qry)
	return rows[0][0] if rows else ""


#########################################
def get_ensembl_exons(cursor, species_db, gene_id):
	exons = []

	qry  = "select distinct exon_transcript.exon_id "
	qry += f"from  {species_db}.exon_transcript, {species_db}.transcript "
	qry += "where exon_transcript.transcript_id = transcript.transcript_id "
	qry += f"and transcript.gene_id = {gene_id} "

	rows = error_intolerant_search(cursor, qry)

	if not rows: return exons
	exon_ids = [row[0] for row in rows]

	# get the region on the gene
	ret = get_gene_coordinates(cursor, species_db, gene_id)
	if ret:
		[gene_seq_id, gene_region_start, gene_region_end,
		 gene_region_strand] = ret
	else:
		print(f"region not retrived for {species_db}  {gene_id}")
		return []

	column_names = get_column_names(cursor, species_db, "exon")
	for exon_id in exon_ids:
		exon = Exon()
		exon.load_from_db(cursor, species_db, "exon", exon_id, column_names=column_names)
		if not exon.exon_id: continue
		exons.append(exon)

	return exons


#########################################
def get_predicted_exons(cursor, species_db, gene_id):
	exons = []

	print("check out the implentation of get_predicted_exons() in enembl.py")
	exit()
	# get the region on the gene
	ret = get_gene_coordinates(cursor, species_db, gene_id)
	if ret:
		[gene_seq_id, gene_region_start, gene_region_end, gene_region_strand] = ret
	else:
		print("region not retrived for ", species_db, gene_id)
		return []

	qry  = f"select * from {species_db}.prediction_exon  where seq_region_id = {gene_seq_id} "
	qry += f"and  seq_region_start >= {gene_region_start} and seq_region_start <= {gene_region_end} "
	qry += f"and  seq_region_end   >= {gene_region_start}  and seq_region_end   <= {gene_region_end} "
	rows = search_db(cursor, qry)

	if not rows: return []

	column_names = get_column_names(cursor, species_db, "prediction_exon")
	for row in rows:
		exon = Exon()
		exon.gene_id = gene_id
		#exon.load_from_db(cursor, species_db, "prediction_exon", exon_id, column_names=column_names)
		exons.append(exon)

	return exons


#########################################
def get_novel_exons(cursor, gene_id, table):
	exons = []

	qry = "select * from %s " % table
	qry += " where gene_id = %d " % int(gene_id)
	rows = search_db(cursor, qry)
	if not rows: return exons

	for row in rows:
		exon = Exon()
		exon.load_from_novel_exon(row, table)
		exons.append(exon)
	return exons


#########################################
def get_primary_seq_info(cursor, gene_id, species):
	# seq identifier from gene table
	qry = "select seq_region_id, seq_region_start, seq_region_end, seq_region_strand"
	qry += " from gene where gene_id = %d" % gene_id
	rows = search_db(cursor, qry)
	if (not rows):
		search_db(cursor, qry, verbose=True)
		exit(1)
	[seq_region_id, seq_region_start, seq_region_end, seq_region_strand] = rows[0]

	qry = "select name, file_name from seq_region "
	qry += " where seq_region_id= %d" % seq_region_id
	rows = search_db(cursor, qry)
	if (not rows):
		search_db(cursor, qry, verbose=True)
		return []
	[seq_name, file_names] = rows[0]
	# Ivana:
	# What is the proper way to find out whether the seq_region is mitochondrial?
	# EMily:
	# It's in the seq_region table, under name. Name will either be a chromosome
	# number, MT for mitochrondria or the name of a contig.
	mitochondrial = is_mitochondrial(cursor, gene_id)

	return [seq_name, file_names, seq_region_start,
			seq_region_end, seq_region_strand, mitochondrial]


#########################################
def get_alt_seq_info(cursor, gene_id, species):
	# seq identifier from gene table
	qry = "select seq_region_id, seq_region_start, seq_region_end,  seq_region_strand from gene where gene_id = %d" % gene_id
	rows = search_db(cursor, qry)
	if not rows: return []

	[seq_region_id, orig_seq_region_start, orig_seq_region_end, seq_region_strand] = rows[0]

	# check whether we have "assembly exception"
	# we do not want 'PAR' regions, though:
	'''
	The pseudo-autosomal regions are homologous DNA sequences on the (human) X and Y chromosomes. 
	They allow the pairing and crossing-over of these sex chromosomes the same way the autosomal 
	chromosomes do during meiosis. 
	As these genomic regions are identical between X and Y, they are oftentimes only stored once.
	The exception types we are interested in are PATCH_FIX and  PATCH_NOVEL
	'''
	qry = "select  seq_region_start,  seq_region_end, exc_seq_region_id, exc_seq_region_start,  exc_seq_region_end "
	qry += "from assembly_exception where seq_region_id = %d " % seq_region_id
	qry += "and  assembly_exception.exc_type  like 'PATCH_%'"
	rows = search_db(cursor, qry)
	if not rows: return []
	[seq_region_start, seq_region_end, exc_seq_region_id, exc_seq_region_start, exc_seq_region_end] = rows[0]

	qry = "select name, file_name from seq_region where seq_region_id= %d" % exc_seq_region_id
	rows = search_db(cursor, qry)
	if not rows: return []
	[seq_name, file_names] = rows[0]

	mitochondrial = is_mitochondrial(cursor, gene_id)

	return [seq_name, file_names, seq_region_start, seq_region_end, seq_region_strand, mitochondrial]


#########################################
def extract_gene_seq(acg, species, seq_name, file_names, seq_region_strand,
					 seq_region_start, seq_region_end):
	# now the question is, which file do I use if there are several options?
	first_choice = ""
	second_choice = ""
	for file_name in file_names.split(" "):
		if '.chromosome.' in file_name:
			first_choice = file_name
		elif '.toplevel.' in file_name:
			second_choice = file_name

	fasta_db_file = ""
	if first_choice:
		fasta_db_file = first_choice
	elif second_choice:
		fasta_db_file = second_choice

	if not fasta_db_file:
		print("failed to decide on fasta_db_file:")
		print(file_names)
		exit(1)

	# extract gene sequence  from fasta db
	fastacmd = acg.generate_fastacmd_gene_command(species, seq_name, fasta_db_file,
												  seq_region_strand, seq_region_start,
												  seq_region_end)
	ret = subprocess.getoutput(fastacmd)
	if not ret:
		#print("no refturn for fastacmd for", species, gene_id)
		print("fastacmd: ", fastacmd)
		exit(1)
	if ('ERROR' in ret):
		print("Error running fastacmd: ", fastacmd)
		print(ret)
		if 'Ignoring sequence location' in ret:
			print('will ignore')
			print()
		else:
			exit(1)
	gene_seq = ""
	reading = 0
	for line in ret.split("\n"):
		if ('>' in line):
			reading = 1
			continue
		if (not reading):
			continue
		line.rstrip()
		gene_seq += line

	return [gene_seq, fasta_db_file]


#########################################
def get_transcript_ids(cursor, species_db, gene_id):
	transcript_ids = []

	qry = f"select transcript_id, stable_id from {species_db}.transcript where gene_id = {gene_id}"
	rows = search_db(cursor, qry)

	if not rows:
		return transcript_ids

	for row in rows:
		transcript_ids.append(row)

	return transcript_ids


#########################################
def get_translation_coords(cursor, transcript_id, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return []

	qry = "select  seq_start, start_exon_id, seq_end, end_exon_id  "
	qry += " from translation where transcript_id = %d " % transcript_id

	rows = search_db(cursor, qry)

	if not rows:
		return []

	if 'error' in str(rows[0][0]).lower():
		return []

	return rows[0]


#########################################
def get_analysis_dict(cursor):
	source = {}
	qry = "select analysis_id, logic_name  from analysis"
	rows = search_db(cursor, qry)
	if (not rows):
		print("blah?")
		return False
	for row in rows:
		source[row[0]] = row[1]
	return source


#########################################
def exon_id2gene_id(cursor, ensembl_db_name, exon_id, is_known):
	switch_to_db(cursor, ensembl_db_name)
	if is_known == 3:  # sw_sharp exon
		qry = "select gene_id from usearch_exon where "
		qry += "exon_id = %d " % exon_id
	elif is_known == 2:  # sw_sharp exon
		qry = "select gene_id from sw_exon where "
		qry += "exon_id = %d " % exon_id
	else:
		qry = "select gene_id from gene2exon where "
		qry += "exon_id = %s and is_known = %s " % (exon_id, is_known)

	rows = search_db(cursor, qry)
	if (not rows or 'ERROR' in rows[0]):
		rows = search_db(cursor, qry, verbose=True)
		return ""

	return rows[0][0]


#########################################
def get_paras(cursor, gene_id):
	paras = []
	qry = "select cognate_gene_id from paralogue "
	qry += " where gene_id = %d " % gene_id

	rows = search_db(cursor, qry)
	if (not rows):
		return []

	for row in rows:
		para_gene_id = row[0]
		paras.append(para_gene_id)

	return paras


#########################################
def get_orthos(cursor, gene_id, table):
	orthos = []
	qry = "select cognate_gene_id, cognate_genome_db_id from " + table
	qry += " where gene_id = %d" % gene_id

	rows = search_db(cursor, qry)
	if (not rows):
		return []

	for row in rows:
		ortho_gene_id = row[0]
		ortho_genome_db = row[1]
		# note: the cursor will be pointing to compara db after this
		species = genome_db_id2species(cursor, ortho_genome_db)
		orthos.append([ortho_gene_id, species])

	return orthos


#########################################
def get_description(cursor, gene_id, db_name=None):
	if db_name and not switch_to_db(cursor, db_name):
		return False
	try:
		gene_id = int(gene_id)
		qry = "select description from gene where gene_id = %d " % gene_id
	except ValueError as verr:
		if type(gene_id) == str and gene_id[:3] == 'ENS':
			qry = "select description from gene where stable_id='%s' " % gene_id
		else:
			return False
	ret = error_intolerant_search(cursor, qry)
	return ret[0][0] if ret else False


#########################################
def get_logic_name(cursor, analysis_id, db_name=None):
	if analysis_id < 0:
		return ''

	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False
	qry = "SELECT logic_name FROM analysis WHERE analysis_id = %d" % analysis_id
	rows = search_db(cursor, qry)
	if (not rows):
		logic_name = ''
	else:
		logic_name = rows[0][0]
	return logic_name


#########################################
def is_coding(cursor, gene_id, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False

	qry = "select is_coding from gene where gene_id = %d " % int(gene_id)
	rows = search_db(cursor, qry)
	if not rows:
		return False

	return rows[0][0] > 0


#########################################
def get_gene_biotype(cursor, gene_id, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False
	if type(gene_id)==str and 'ENS' in gene_id:
		qry = "select biotype from gene where stable_id = '%s' " % gene_id
	else:
		qry = "select biotype from gene where gene_id = %d " % int(gene_id)
	rows = hard_landing_search(cursor, qry)

	return rows[0][0]


#########################################
def get_biotype(cursor, exon_id, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False

	qry = "select biotype from gene where gene_id = %d " % int(exon_id)
	rows = search_db(cursor, qry)
	if not rows:
		return ""

	return rows[0][0]


#########################################
def get_status(cursor, exon_id, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False

	qry = "select status from gene where gene_id = %d " % int(exon_id)
	rows = search_db(cursor, qry)
	if (not rows):
		return False

	return rows[0][0]


#########################################
def is_coding_exon(cursor, exon_id, is_known, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False
	qry = "select is_coding from gene2exon where exon_id = %d and is_known = %d" % (exon_id, is_known)
	rows = search_db(cursor, qry)
	if (not rows):
		return False

	return rows[0][0] > 0


#########################################
def get_gene_coordinates(cursor, db_name, gene_id):
	qry  = "select seq_region_id, seq_region_start, seq_region_end, seq_region_strand  "
	qry += f"from {db_name}.gene  where gene_id = {gene_id}"
	return hard_landing_search(cursor, qry)[0]  # we will not tolerate gene without coordinates


def get_gene_start(cursor, db_name, gene_id):
	return get_gene_coordinates(cursor, db_name, gene_id)[1]


#########################################
def is_mitochondrial(cursor, gene_id, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False

	qry = "select seq_region.name from seq_region, gene  "
	qry += " where seq_region.seq_region_id =  gene.seq_region_id "
	qry += " and gene.gene_id = %d" % gene_id
	rows = search_db(cursor, qry)
	if (not rows):
		search_db(cursor, qry, verbose=True)
		return []
	seq_name = rows[0][0]

	is_mitochondrial = (seq_name == 'MT')

	return is_mitochondrial


#########################################
def get_exon_pepseq(cursor, db_name, exon_id):

	qry = f"select protein_seq from {db_name}.exon_seq where exon_id = {exon_id}"
	rows = error_intolerant_search(cursor,qry)
	if not rows: return ""

	protein_seq = rows[0][0]
	if protein_seq is None: protein_seq = ""

	return protein_seq


#########################################
def exon_seq_id2exon_id(cursor, exon_seq_id, db_name=None):

	if db_name and not switch_to_db(cursor, db_name): return ""
	qry = "select exon_id, is_known from exon_seq where exon_seq_id = %d  " % int(exon_seq_id)
	rows = search_db(cursor, qry)
	if not rows: return ""

	return rows[0]


#########################################
def get_exon_dnaseq(cursor, exon, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return None

	if exon.analysis_id > 0:
		exon_id = exon.exon_id
		is_known = exon.is_known
		qry = "select dna_seq  "
		qry += " from exon_seq where exon_id = %d and is_known = %d" % (exon_id, is_known)
	else:
		exon_seq_id = exon.exon_seq_id
		qry = "select dna_seq "
		qry += " from  exon_seq where exon_seq_id = %d" % exon_seq_id

	rows = search_db(cursor, qry)
	if (not rows):
		# rows = search_db(cursor, qry, verbose = True)
		return None

	dna_seq = rows[0][0]

	return dna_seq


#########################################
def get_sw_seq_id(cursor, exon_id, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return -1
	qry = "select exon_seq_id "
	qry += " from sw_exon where exon_id = %d" % exon_id
	rows = search_db(cursor, qry)
	if not rows or not rows[0][0]:
		return -1

	return int(rows[0][0])


#########################################
def get_exon_seq_by_db_id(cursor, exon_seq_id, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False

	qry = "select exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
	qry += " left_flank, right_flank, dna_seq  "
	qry += " from exon_seq where exon_seq_id = %d" % exon_seq_id
	rows = search_db(cursor, qry)
	if (not rows):
		# rows = search_db(cursor, qry, verbose = True)
		return []

	[exon_seq_id, protein_seq, pepseq_transl_start,
	 pepseq_transl_end, left_flank, right_flank, dna_seq] = rows[0]
	if (protein_seq is None):
		protein_seq = ""
	if (left_flank is None):
		left_flank = ""
	if (right_flank is None):
		right_flank = ""
	if (dna_seq is None):
		dna_seq = ""

	return [exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, left_flank, right_flank, dna_seq]


#########################################
def get_exon_phase(cursor, exon_id, is_known, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False

	if is_known == 2:  # sw exon
		qry = "select phase from sw_exon"
	elif is_known == 1:
		qry = "select phase from exon"
	elif is_known == 0:
		qry = "select start_phase from prediction_exon"
	else:
		return None

	qry += " where exon_id = %d " % exon_id
	rows = search_db(cursor, qry)

	if not rows:
		rows = search_db(cursor, qry, verbose=True)
		return None

	return rows[0][0]


#########################################
def get_exon_seqs(cursor, exon_id, is_known, db_name=None):
	if (db_name):
		if not switch_to_db(cursor, db_name):
			return False

	if is_known == 2:  # sw exon
		qry = "select exon_seq.exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
		qry += " left_flank, right_flank, dna_seq  from  exon_seq "
		qry += " join sw_exon on exon_seq.exon_seq_id = sw_exon.exon_seq_id  "
		qry += " where sw_exon.exon_id = %d " % exon_id

	elif is_known == 3:  # usearch exon
		qry = "select exon_seq.exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
		qry += " left_flank, right_flank, dna_seq  from  exon_seq "
		qry += " join usearch_exon on exon_seq.exon_seq_id = usearch_exon.exon_seq_id  "
		qry += " where usearch_exon.exon_id = %d " % exon_id

	else:
		qry = "select exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
		qry += " left_flank, right_flank, dna_seq  "
		qry += " from  exon_seq where exon_id = %d and is_known = %d" % (exon_id, is_known)

	rows = search_db(cursor, qry)

	if not rows or not len(rows[0]) == 7:
		# rows = search_db(cursor, "select database()")
		# print "using db ", rows[0]
		# rows = search_db(cursor, qry, verbose = True)
		return []

	[exon_seq_id, protein_seq, pepseq_transl_start,
	 pepseq_transl_end, left_flank, right_flank, dna_seq] = rows[0]
	if (protein_seq is None):
		protein_seq = ""
	if (left_flank is None):
		left_flank = ""
	if (right_flank is None):
		right_flank = ""
	if (dna_seq is None):
		dna_seq = ""

	return [exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, left_flank, right_flank, dna_seq]


#########################################
def get_exon(cursor, exon_id, is_known=None, db_name=None):
	exon = Exon()

	if (db_name):
		if not switch_to_db(cursor, db_name):
			return exon

	if is_known == 2:
		# sw# exon
		qry = "select * from sw_exon where exon_id = %d" % exon_id
		rows = search_db(cursor, qry, verbose=False)
		if (not rows):
			return exon
		exon.load_from_novel_exon(rows[0], "sw_exon")
	elif is_known == 3:
		# sw# exon
		qry = "select * from usearch_exon where exon_id = %d" % exon_id
		rows = search_db(cursor, qry, verbose=False)
		if (not rows):
			return exon
		exon.load_from_novel_exon(rows[0], "usearch_exon")
	else:
		qry = "select * from gene2exon where exon_id = %d" % exon_id
		if is_known: qry += " and is_known = %s " % is_known
		rows = search_db(cursor, qry, verbose=False)
		if (not rows):
			return exon
		exon.load_from_gene2exon(rows[0])

	return exon



#########################################
def get_selenocysteines(cursor, gene_id):
	selenoC_pos = []

	canonical_transl_id = gene2canon_transl(cursor, gene_id)

	qry = "select value from translation_attrib "
	qry += " where attrib_type_id = 12 and  translation_id = %d " % canonical_transl_id

	rows = search_db(cursor, qry)
	if (not rows):
		return []

	for row in rows:
		blah = row[0].split(" ")
		start = int(blah[0])
		end = int(blah[1])
		for pos in range(start, end + 1):
			selenoC_pos.append(pos - 1)

	return selenoC_pos


#########################################
def gene2stable_canon_transl_id(cursor, gene_id, db_name=None):

	current_db = get_current_db(cursor)
	if db_name and not switch_to_db(cursor, db_name):
		return False

	qry  = "select translation.stable_id  from translation, gene "
	qry += "where gene.canonical_transcript_id = translation.transcript_id "
	qry += "and gene.gene_id = %d " % gene_id
	rows = search_db(cursor, qry, verbose=False)

	if not rows:
		rows = search_db(cursor, qry, verbose=True)
		return ""

	switch_to_db(cursor, current_db)
	return rows[0][0]


#########################################
def gene2canon_transl(cursor, gene_id, db_name=None, stable=False):
	if (db_name and not switch_to_db(cursor, db_name)):
		return False
	field = "translation.stable_id" if stable else "translation.translation_id"
	qry = "select %s from translation, gene " % field
	qry += " where gene.canonical_transcript_id = translation.transcript_id "
	qry += " and gene.gene_id = %d " % gene_id
	rows = search_db(cursor, qry, verbose=False)

	if (not rows):
		rows = search_db(cursor, qry, verbose=True)
		return ""

	return rows[0][0]


def gene_id2hgnc(cursor, gene_id):
	qry = "select h.approved_symbol from identifier_maps.hgnc h, gene g "
	qry += f"where h.ensembl_gene_id=g.stable_id and g.gene_id={gene_id}"
	ret = error_intolerant_search(cursor, qry)
	return ret[0][0] if ret else "unk"


def gene_stable2hgnc(cursor, gene_stable):
	qry = f"select approved_symbol from  identifier_maps.hgnc where ensembl_gene_id='{gene_stable}'"
	ret = error_intolerant_search(cursor, qry)
	return ret[0][0] if ret else "anon"

########
def stable2member(cursor, stable_id):
	# member_id refers to compara db
	# of which we need to have one
	# qry = "select  member_id from member where stable_id = '%s'" % stable_id
	# since version 76 the table is called gene_member, and the main id is gene_member_id
	qry = "select  gene_member_id from gene_member where stable_id = '%s'" % stable_id
	rows = search_db(cursor, qry)
	if (not rows or 'ERROR' in rows[0]):
		rows = search_db(cursor, qry, verbose=True)
		exit(1)
		return ""

	return int(rows[0][0])


########
def member2stable(cursor, member_id):
	# member_id refers to compara db
	# of which we need to have one
	qry = "select  stable_id from gene_member where gene_member_id = %d" % member_id
	rows = search_db(cursor, qry)
	if (not rows):
		rows = search_db(cursor, qry, verbose=True)
		return ""

	return rows[0][0]


################
def time_qry(cursor, qry, verbose=True):
	time1 = time0 = 0
	if verbose: time0 = time()
	ret = error_intolerant_search(cursor,qry)
	if verbose: time1 = time()
	if verbose: print("\n%s\ndone in %.3f s" % (qry, float(time1-time0)))
	return ret


################
# this was way too slow - check out the hack in 18_orthologues to go orund thins
# ./el_utils/kernprof.py -l <calling script>.py
# python3 -m line_profiler <calling script>.py.lprof
# @profile
# def get_orthologues(cursor, compara_db, gene_member_id, verbose=False):
# 	switch_to_db(cursor, compara_db)
# 	qry = "select  homology.homology_id, homology.description from  homology,  homology_member"
# 	qry += " where homology_member.gene_member_id =%d " % gene_member_id
# 	qry += " and homology.homology_id = homology_member.homology_id "
# 	ret = time_qry(cursor,qry, verbose)
# 	if not ret or len(ret)==0: return {}
# 	# hom type can be ortholog_one2one, ortholog_one2many, apparent_ortholog_one2one
# 	# ortholog_one2many, ortholog_many2many
# 	homology_type = dict(ret)
# 	orthos = {}
# 	for hom_id, hom_type in homology_type.items():
# 		if not "ortho" in hom_type:continue
# 		if not hom_type in orthos: orthos[hom_type] = []
#
# 		qry = "select gene_member_id from homology_member "
# 		qry += " where homology_id = %d" % hom_id
# 		qry += " and not  gene_member_id = %d" % gene_member_id
# 		ortho_id = time_qry(cursor,qry, verbose)[0][0] # there should be ony one other member of the pair
#
# 		qry = "select  gene_member.stable_id, genome_db.name, genome_db.genome_db_id "
# 		qry += " from gene_member, genome_db "
# 		qry += " where gene_member.gene_member_id = %d " % ortho_id
# 		qry += " and genome_db.genome_db_id = gene_member.genome_db_id"
# 		orthos[hom_type].append(time_qry(cursor, qry, verbose)[0])
# 	return orthos


########################
def get_orthologues_from_species(cursor, ensembl_db_name, ortho_type, gene_member_id, species):
	# the ortho_type is one of the following: 'ortholog_one2one',
	# 'ortholog_one2many', 'ortholog_many2many', 'possible_ortholog', 'apparent_ortholog_one2one'
	orthos = []

	# find genome db_id
	genome_db_id = species2genome_db_id(cursor, species)

	# make the cursor point to compara database - should be the responsibility of each function
	switch_to_db(cursor, get_compara_name(cursor))

	qry = "select homology.homology_id from homology_member, homology "
	qry += " where homology_member.gene_member_id =%d " % gene_member_id
	qry += " and homology.homology_id = homology_member.homology_id "
	qry += " and  homology.description = '%s' " % ortho_type
	rows = search_db(cursor, qry)

	if (not rows):
		return []  # no orthologs here

	# for each homology id find the other member id
	# print qry
	# print member_id, ortho_type, species, genome_db_id
	# print rows
	for row in rows:
		homology_id = row[0]
		# print "\t homology id:", homology_id
		switch_to_db(cursor, get_compara_name(cursor))
		qry = "select gene_member_id from homology_member "
		qry += " where homology_id = %d" % int(homology_id)
		qry += " and not  gene_member_id = %d" % gene_member_id

		rows2 = search_db(cursor, qry, verbose=False)
		if (not rows2):
			# print "\t ",
			# rows2 = search_db (cursor, qry, verbose = True)
			continue
		for row2 in rows2:
			ortho_id = row2[0]
			# print "\t\t ortho id:", ortho_id
			qry = "select  stable_id  from gene_member  "
			qry += " where gene_member_id = %d " % ortho_id
			qry += " and genome_db_id = %d " % genome_db_id
			rows3 = search_db(cursor, qry, verbose=False)
			if (not rows3):
				# print "\t\t ",
				# rows3 = search_db (cursor, qry, verbose = True)
				continue
			ortho_stable = rows3[0][0]
			# print "\t\t ortho stable:", ortho_stable
			orthos.append(ortho_stable)
	if orthos:
		switch_to_db(cursor, ensembl_db_name[species])
		orthos = [stable2gene(cursor, gene_id) for gene_id in orthos]
	# print 'orthos:', orthos
	return orthos


def canonical_protein2hgnc_symbol(cursor, protein_stable_id):

	qry  = "select  approved_symbol from identifier_maps.hgnc h, gene g, transcript mrna, translation prot "
	qry += "where h.ensembl_gene_id=g.stable_id and g.gene_id=mrna.gene_id and "
	qry += "mrna.canonical_translation_id=prot.translation_id "
	qry += f"and prot.stable_id='{protein_stable_id}'"
	ret = error_intolerant_search(cursor, qry)
	if not ret:
		print("[Warning] HGNC symbol not found for", protein_stable_id)
		return protein_stable_id

	return ret[0][0]


########
def stable2gene(cursor, stable_id=None, db_name=None):
	if (not stable_id):
		return 0

	if (db_name and not switch_to_db(cursor, db_name)):
		return False

	qry = "select gene_id from gene where stable_id='%s'" % stable_id
	rows = search_db(cursor, qry, verbose=False)

	if (not rows):
		rows = search_db(cursor, qry, verbose=True)
		return 0

	return int(rows[0][0])


########
def gene2stable(cursor, gene_id=None, db_name=None):
	if not gene_id: return ""

	if (db_name):
		qry = "use %s " % db_name
		rows = search_db(cursor, qry)
		if (rows):
			rows = search_db(cursor, qry, verbose=True)
			print(rows)
			exit(1)

	qry = "select stable_id from gene where gene_id=%d" % gene_id
	rows = search_db(cursor, qry, verbose=False)

	if (not rows):
		rows = search_db(cursor, qry, verbose=True)
		return ""

	return rows[0][0]


########
def gene_member2stable(cursor, gene_id=None, db_name=None):
	if not gene_id: return ""
	qry = "select stable_id "
	if not db_name:  # cursor is already set to compara db
		qry += "from gene_member "
	else:
		qry += "from %s.gene_member " % db_name
	qry += "where gene_member_id=%d" % gene_id
	return hard_landing_search(cursor, qry)[0][0]


########
def exon2stable(cursor, exon_id=None, db_name=None):
	if (not exon_id):
		return ""

	if (db_name):
		qry = "use %s " % db_name
		rows = search_db(cursor, qry)
		if (rows):
			rows = search_db(cursor, qry, verbose=True)
			print(rows)
			exit(1)

	qry = "select stable_id from exon where exon_id=%d" % exon_id
	rows = search_db(cursor, qry, verbose=False)

	if (not rows):
		rows = search_db(cursor, qry, verbose=True)
		return ""

	return rows[0][0]


########
def stable2exon(cursor, stable_id, db_name=None):
	if (db_name):
		qry = "use %s " % db_name
		rows = search_db(cursor, qry)
		if (rows):
			rows = search_db(cursor, qry, verbose=True)
			print(rows)
			exit(1)

	qry = "select exon_id from exon where stable_id='%s'" % stable_id
	rows = search_db(cursor, qry, verbose=False)

	if (not rows):
		rows = search_db(cursor, qry, verbose=True)
		return ""
	return rows[0][0]


########
def is_reference(cursor, gene_id, non_ref_id, db_name=None, stable=False):
	if (db_name):
		qry = "use %s " % db_name
		rows = search_db(cursor, qry)
		if (rows):
			rows = search_db(cursor, qry, verbose=True)
			print(rows)
			exit(1)

	if stable:
		qry = "select seq_region_id from gene where stable_id='%s'" % gene_id
	else:
		qry = 'select seq_region_id from gene where gene_id=%d' % int(gene_id)
	rows = search_db(cursor, qry)
	if not rows: return True

	seq_region_id = rows[0]

	qry = 'select attrib_type_id from  seq_region_attrib where seq_region_id=%s' % seq_region_id
	rows = search_db(cursor, qry)
	if not rows: return True

	for row in rows:
		attrib_type_id = int(row[0])
		if attrib_type_id == non_ref_id: return False

	return True


########
def get_gene_ids(cursor, db_name=None, biotype=None,  stable=False, ref_only=False):
	if db_name and not switch_to_db(cursor, db_name): exit()

	gene_ids = []

	qry = "select stable_id from gene " if stable else "select gene_id from gene "
	if biotype:
		qry += " where biotype='%s'" % biotype
	rows = search_db(cursor, qry, verbose=False)

	if not rows:
		search_db(cursor, qry, verbose=True)
		return []
	else:
		if 'Error' in rows[0]:
			print(qry)
			print(rows)
			return []

		# I don't want to hard code the id for the annotation "non_ref"
		# and I do no know where else to do it, so qe do it here
		non_ref_id = 0
		if ref_only:
			qry = "select attrib_type_id from attrib_type where code='non_ref'"
			rows2 = search_db(cursor, qry)
			if not rows2 or not type(rows2[0][0]) is int:
				ref_only = False
			else:
				non_ref_id = int(rows2[0][0])

		for row in rows:
			gene_id = row[0]
			if not stable: gene_id = int(gene_id)
			if ref_only and not is_reference(cursor, gene_id, non_ref_id, stable): continue
			gene_ids.append(gene_id)

	return gene_ids


########
def get_ensembl_species(cursor):
	exceptions = {'cebus_capucinus_imitator':'cebus_capucinus',
					'gorilla_gorilla_gorilla':'gorilla_gorilla'}
	ensembl_db_name = {}
	common_name = {}
	all_species = []

	# find the release number
	qry = "select value from exolocator_meta.parameters where name = 'ensembl_release_number'"
	release_number = hard_landing_search(cursor, qry)[0][0]

	qry = f"show databases like '%core_{release_number}%'"
	for row in hard_landing_search(cursor, qry):
		db_name = row[0]

		qry = f"select meta_value from {db_name}.meta where meta_key='species.scientific_name'"
		species = hard_landing_search(cursor, qry)[0][0].lower().replace(" ", "_").replace("'", "_")
		if species in exceptions: species = exceptions[species]
		ensembl_db_name[species] = db_name

		qry = f"select meta_value from {db_name}.meta where meta_key='species.common_name'"
		species_common = hard_landing_search(cursor, qry)[0][0].lower().replace(" ", "_").replace("'", "_")
		common_name[species] = species_common

		all_species.append(species)

	return all_species, ensembl_db_name, common_name


def get_species(cursor):
	qry = "select species_name, db_name from exolocator_meta.db_names "
	ensembl_db_name = dict(hard_landing_search(cursor, qry))
	return sorted(list(ensembl_db_name.keys())), ensembl_db_name

########
def get_compara_name(cursor):
	# find the release number
	qry = "select value from exolocator_meta.parameters where name = 'ensembl_release_number'"
	release_number = hard_landing_search(cursor, qry)[0][0]

	qry = f"show databases like '%compara_{release_number}%'"
	rows = search_db(cursor, qry)
	if not rows:
		search_db(cursor, qry, verbose=True)
		return ""

	return rows[0][0]


########
def species2taxid(cursor, species):
	qry = "select tax_id from exolocator_meta.species_names where species = '%s'" % species
	tax_id = hard_landing_search(cursor, qry)[0][0]
	return tax_id


########
def ncbi_species2taxid(cursor, species):
	species_txt = species.replace('_', ' ')
	qry  = "select taxon_id from ncbi_taxonomy.ncbi_taxa_name "
	qry += f"where name='{species_txt}' "
	tax_id = hard_landing_search(cursor, qry)[0][0]
	return tax_id


########
def species2genome_db_id(cursor, species):

	compara_name = get_compara_name(cursor)
	qry = f"select genome_db_id from {compara_name}.genome_db where name = '{species}'"

	rows = error_intolerant_search(cursor, qry)
	return int(rows[0][0]) if rows else 0


########
def genome_db_id2species(cursor, genome_db_id):

	compara_name = get_compara_name(cursor)
	qry = f"select name from {compara_name}.genome_db where genome_db_id ={int(genome_db_id)}"
	rows = error_intolerant_search(cursor, qry)
	return rows[0][0] if rows else ""

