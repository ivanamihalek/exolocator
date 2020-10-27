
import os
from el_utils.el_specific import *

# from  el_utils.py_alignment import  smith_waterman, exon_aware_smith_waterman
# pycharm does not recognize this, but python itself does
from c_utils.alignment import smith_waterman_context

#########################################################
class Map:    # this particular map is between exons

	#################################################
	# print
	def __str__(self):

		printstr = ""

		for attr, value in self.__dict__.items():
			if value is None:
				printstr += attr + "    None"
			else:
				printstr += attr + "    " + str(value)
			printstr += "\n"

		return printstr

	###################################
	# when something is defined as an exon ....
	def __init__(self):

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
def get_maps(cursor, species_db_name, exon_id, paralogue=False):

	maps = []
	table = "paralogue_exon_map" if paralogue else "exon_map"
	qry  = f"select * from {species_db_name}.{table} where exon_id = {exon_id} "
	rows = error_intolerant_search (cursor, qry)
	if not rows: return []

	for row in rows:
		map = Map()
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
	relevant_human_exons = find_relevant_exons_and_pepseqs (cursor, human_exons, ensembl_db_name['homo_sapiens'], ref_species=True)
	for he in relevant_human_exons:
		map = Map()
		map.species_1    = 'homo_sapiens'
		map.species_2    = 'homo_sapiens'
		map.exon_id_1    = he['exon_id']
		map.exon_id_2    = he['exon_id']

		map.exon_known_1 = he['is_known']
		map.exon_known_2 = he['is_known']
		map.similarity   = 1.0

		pepseq = get_exon_pepseq (cursor, he, verbose=True)
		map.cigar_line   = cigar_line(pepseq, pepseq)
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
def maps_evaluate (min_similarity, rep_species, rep_species_exons, ortho_exons, aligned_seq, exon_positions, verbose=False):

	maps = []

	if len(list(aligned_seq.keys())) > 2:
		print("right now the mapping implemented for two species only")
		return []

	other_species = [k for k in list(aligned_seq.keys()) if k!=rep_species][0] # this is two-seq_alignment

	if verbose:
		print('==============================================================')
		print(rep_species, "        ", other_species)


	for rep_exon_ct in range(len(rep_species_exons)):

		padded_count_human = "{0:03d}".format(rep_exon_ct+1)
		if padded_count_human not in exon_positions[rep_species]: continue
		[rep_start, rep_end] = exon_positions[rep_species][padded_count_human]

		for ortho_exon_ct in range(len(ortho_exons)):

			padded_count_ortho = "{0:03d}".format(ortho_exon_ct+1)
			if padded_count_ortho not in exon_positions[other_species]:
				continue
			[other_start, other_end] = exon_positions[other_species][padded_count_ortho]

			if overlap(rep_start, rep_end, other_start, other_end):

				map = Map()
				map.species_1    = rep_species
				map.species_2    = other_species

				map.exon_id_1    = rep_species_exons[rep_exon_ct]['exon_id']
				map.exon_id_2    = ortho_exons[ortho_exon_ct]['exon_id']

				map.exon_known_1 = rep_species_exons[rep_exon_ct]['is_known']
				map.exon_known_2 = ortho_exons[ortho_exon_ct]['is_known']

				if ortho_exons[ortho_exon_ct]['is_known'] == 2:
					map.source   = 'sw_sharp'
				elif ortho_exons[ortho_exon_ct]['is_known'] == 3:
					map.source   = 'usearch'
				else:
					map.source   = 'ensembl'

				exon_seq_human   = aligned_seq[rep_species][rep_start:rep_end].replace('#','-')
				exon_seq_other   = aligned_seq[other_species][other_start:other_end].replace('#','-')
				[seq_human, seq_other] = pad_the_alnmt (exon_seq_human,rep_start,
														exon_seq_other, other_start)
				seq = {'rep':seq_human, 'other':seq_other}
				seq = strip_gaps(seq)
				if not seq:
					c=inspect.currentframe()
					print(" in %s:%d" % ( c.f_code.co_filename, c.f_lineno))
					return []

				map.similarity = pairwise_tanimoto(seq['rep'], seq['other'])

				if verbose:
					print('===============================')
					print(seq['rep'])
					print(seq['other'])
					print(map)

				if map.similarity < min_similarity: continue

				map.cigar_line = cigar_line(seq['rep'], seq['other'])

				maps.append(map)

	return maps


##############################
def get_exon_pepseqs_batch(cursor, exon_list, db_name=None):
	exon_string = ",".join([str(exon['exon_id']) for exon in exon_list])
	rows = error_intolerant_search(cursor, f"select exon_id, protein_seq from {db_name}.exon_seq where exon_id in ({exon_string})")
	if not rows: return {}
	return dict(rows)


#########################################
def find_relevant_exons_and_pepseqs (cursor, all_exons, db_name, ref_species):

	relevant_exons = []
	protein_seq    = []

	# 1) choose exons that I need
	for exon in all_exons:
		# let's just stick to canonical exons for now - there is too musch stuff going on
		if exon['is_canonical'] and not exon['is_coding']:
			continue
		# if exon['covering_exon'] > 0:
		# 	continue
		# if ref_species and not exon['is_coding']:
		# 	continue
		# if not exon['exon_seq_id']:
		# 	continue
		relevant_exons.append(exon)

	# 2) sort them by their start position in the gene
	to_remove = []
	pepseqs = get_exon_pepseqs_batch(cursor, relevant_exons, db_name=db_name)
	if not pepseqs: return relevant_exons

	relevant_exons.sort(key=lambda exon: exon['start_in_gene'])
	for i in range(len(relevant_exons)):
		exon   = relevant_exons[i]
		pepseq = pepseqs.get(exon['exon_id'], None)
		if not pepseq:
			to_remove.append(i)
			continue
		pepseq = pepseq.replace('X', '')
		if not pepseq:
			to_remove.append(i)
		else:
			exon['pepseq'] = pepseq

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
def decorate_and_concatenate (exons):
	if not exons: return
	# assumes that the exons are sorted in the order of appearance on main strand
	# also assume that they are all coded on the same strand
	if exons[0]['strand']<0: exons.reverse()
	# get rid of the stop codon
	if exons[-1].get('pepseq', None):
		if exons[-1]['pepseq'][-1]=='*': exons[-1]['pepseq'] =  exons[-1]['pepseq'][:-1]
	# formatting mini-language https://docs.python.org/3.8/library/string.html#formatspec
	# this is adding the exon number to the alignment (am I hacking here or what)
	# for exaple, for pepseq = "JKHHG", count=5
	# 'B'+padded_count+pepseq+'Z' is 'B005JKHHGZ'
	decorated_seq = ""
	count = 1
	for exon in exons:
		pepseq = exon.get('pepseq', '')
		padded_count = "{0:03d}".format(count)
		decorated_seq += 'B'+padded_count+pepseq+'Z'
		count += 1

	return decorated_seq

#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./el_utils/kernprof.py -l this.py
# followed by
# python3 -m line_profiler this.py.lprof
# @profile
#########################################
def make_maps(cursor, ensembl_db_name, rep_species, ortho_species, relevant_exons_rep_species, ortho_exons, verbose=False):

	maps = []
	if not relevant_exons_rep_species: return maps

	if verbose: print("##############################", relevant_exons_rep_species, ortho_species)
	switch_to_db(cursor,  ensembl_db_name[ortho_species])
	relevant_ortho_exons = find_relevant_exons_and_pepseqs(cursor, ortho_exons,  ensembl_db_name[ortho_species], ref_species=False)
	if not relevant_ortho_exons: return maps

	#
	rep_species_seq = decorate_and_concatenate(relevant_exons_rep_species)
	ortho_seq = decorate_and_concatenate(relevant_ortho_exons)
	if verbose: print(rep_species_seq)
	if verbose: print(ortho_seq)
	#
	if not rep_species_seq or not ortho_seq: return maps

	aligned_seq = {rep_species: smith_waterman_context(rep_species_seq, ortho_seq, -5, -3)[0],
	               ortho_species: smith_waterman_context(rep_species_seq, ortho_seq, -5, -3)[1]}
	if verbose: print()
	if verbose: print(f">{rep_species}\n{aligned_seq[rep_species]}")
	if verbose: print(f">{ortho_species}\n{aligned_seq[ortho_species]}")
	if not aligned_seq[rep_species] or not aligned_seq[ortho_species]: return []


	# find the positions of the exons in the alignment
	exon_positions = {}
	for species, seq in aligned_seq.items():
		# move B to beginning of each exon sequence
		seq = moveB(seq)
		# beginning and end of each exon in the alignment
		exon_positions[species] = find_exon_positions(seq)

		if verbose: print()
		if verbose: print(species, exon_positions[species])
	# fill in the actual map values
	# TODO what should be the min exon similarity? store somehwere in config
	# maps_evaluate (min_similarity, rep_species, rep_species_exons, ortho_exons, aligned_seq, exon_positions, verbose=False)
	maps = maps_evaluate(0.8, rep_species, relevant_exons_rep_species, relevant_ortho_exons, aligned_seq, exon_positions, verbose=False)
	if verbose:
		for map in maps:
			print(map)

	return maps


# if __name__=="__main__":
# 	return
