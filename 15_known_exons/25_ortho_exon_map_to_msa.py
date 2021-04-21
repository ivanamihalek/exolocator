#!/usr/bin/python3 -u
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, subprocess, re, sys
from el_utils.utils import erropen, output_fasta
from el_utils.map import get_maps, map2exon, Map, get_sorted_canonical_exons, find_relevant_exons_and_pepseqs
from el_utils.special_gene_sets import *
from el_utils.processes import parallelize, get_process_id
from config import Config
from random import sample


from time import time
from Bio import SeqIO
from bitstring import Bits

verbose = False


#########################################
def print_description(cursor, gene_id, rep_genome_db_id, representative_species_db, human_db):
	stable_id = gene2stable(cursor, gene_id)
	qry  = f"select gene_id from {human_db}.orthologues "
	qry += f"where cognate_genome_db_id= {rep_genome_db_id} and cognate_gene_id={gene_id}"
	human_cognate = hard_landing_search(cursor, qry)[0][0]
	description = get_description(cursor, gene_id, representative_species_db)

	hgnc = gene_stable2hgnc(cursor, gene2stable(cursor, human_cognate, human_db))

	exons_in_human = get_sorted_canonical_exons(cursor, human_db, human_cognate)
	relevant_exons_rep_species = find_relevant_exons_and_pepseqs(cursor, exons_in_human, human_db, ref_species=True)
	print(stable_id, hgnc, f"exons in human: {len(relevant_exons_rep_species)}  ", description, )


#########################################
def concatenate_exons(cursor, ensembl_db_name, sequences, exons_per_species):
	concatenated = {}

	# are there multiple candidates from the same s pecies?
	for species, exon_labels in exons_per_species.items():
		if len(exon_labels) < 2: continue
		# if there are multiple candidates from the same species - are they from the same gene?
		exons_per_gene = {}
		for [exon_id, exon_known_code] in exon_labels:
			gene_id = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known_code)
			if not gene_id in list(exons_per_gene.keys()):
				exons_per_gene[gene_id] = []
			exons_per_gene[gene_id].append([exon_id, exon_known_code])
		# if yes - do they overlap in the gene? ... I so need to do this whole crap differently
		# the whole idea was no to be doing this here, but the alignment progs (mafft) can screw up here
		# big time
		switch_to_db(cursor, ensembl_db_name[species])
		for gene_id, exons_from_gene in exons_per_gene.items():
			if (len(exons_from_gene) < 2): continue
			# how robust should I be here? how many fragments should I worry about?
			# how about some combinatorial pearls, like 3 exons non overlapping, but 4 overlapping 1 or more ...?
			# and then exons could be overlapping when the translation regions are not ...
			# for now, I'll only offer this patch for the cases when the pieces are non-overlapping
			exons = []
			for [exon_id, exon_known_code] in exons_from_gene:
				exon = get_exon(cursor, exon_id, exon_known_code)
				if not exon: continue
				exons.append(exon)
			# sort by translation start
			exons.sort(key=lambda exon: exon.start_in_gene)
			# is transl_Start < transl_end of the previous exon
			overlap = False
			exon_prev = exons[0]
			for exon in exons[1:]:
				if exon_prev.pepseq_transl_end is None:
					end_prev = exon_prev.end_in_gene
				else:
					end_prev = exon_prev.start_in_gene + exon_prev.pepseq_transl_end
				if exon.pepseq_transl_start is None:
					start_this = exon.start_in_gene
				else:
					start_this = exon.start_in_gene + exon.pepseq_transl_start
				# print species, "prev end:", end_prev, "this start:", start_this
				if end_prev > start_this:
					# yes => overlap
					overlap = True
					break
				exon_prev = exon
			# if they overlap, do nothing - ther are already both in the fasta set
			if overlap: continue
			# if they do not not overlap, concatenate them, and mark them as concatenated
			new_name = species + "_concat_" + str(len(concatenated))
			concatenated[new_name] = []
			concat_seq = ""
			for exon in exons:  # note 'exons' are sorted, 'exons_from_gene' are not
				old_name = "{0}:{1}:{2}".format(species, exon.exon_id, exon.is_known)
				if old_name not in sequences:
					# print "no key ", old_name, "in the original set (?) "
					continue

				concatenated[new_name].append(old_name)
				if concat_seq: concat_seq += "Z"
				concat_seq += sequences[old_name]
				# also remove the original seqs from the alignment
				del sequences[old_name]
			sequences[new_name] = concat_seq

	return concatenated


#########################################
def get_exon_maps(cursor, gene_id, column_names):
	exon_map = []
	# get all exons
	qry = f"select e.* from  gene2exon g, exon_map e where g.exon_id=e.exon_id and g.gene_id={gene_id}"
	ret = error_intolerant_search(cursor, qry)
	if ret: exon_map = [dict(zip(column_names, line)) for line in ret]
	return exon_map


#########################################
def split_concatenated_exons(sequences, concatenated):
	for seq_name in concatenated:
		if not 'concat' in seq_name: continue
		pieces = sequences[seq_name].split('Z')
		del sequences[seq_name]
		if not len(concatenated[seq_name]) == len(pieces): continue
		for piece_name in concatenated[seq_name]:
			piece_seq = ""
			for piece in pieces:
				if piece_seq: piece_seq += '-'  # to make up for 'Z' we have lost in the splitting op above
				if pieces.index(piece) == concatenated[seq_name].index(piece_name):
					piece_seq += piece
				else:
					piece_seq += '-' * len(piece)
			sequences[piece_name] = piece_seq


#########################################
def multiple_alignment_exon(cursor, ensembl_db_name, ref_genome_db_id, alignable_exon_maps):

	pepseq = {}
	if not alignable_exon_maps: return pepseq
	# at this point, all  maps should belong to the same reference exon
	map = alignable_exon_maps[0]
	ref_species = genome_db_id2species(cursor, ref_genome_db_id)
	ref_genome_db = ensembl_db_name[ref_species]
	pepseq[ref_species] = get_exon_pepseq(cursor, ref_genome_db, map['exon_id'])
	print(ref_species, pepseq[ref_species])
	for map in alignable_exon_maps:
		other_species = genome_db_id2species(cursor,  map['cognate_genome_db_id'])
		other_db = ensembl_db_name[other_species]
		if map['similarity'] < Config.min_accptbl_exon_sim:
			pepseq[other_db] = ""
		else:
			pepseq[other_db] = get_exon_pepseq(cursor, other_db, map['cognate_exon_id'])
		print(map)
		print(pepseq[other_db])
	exit()

	# 		pepseq = get_exon_pepseq(cursor, exon)
	# 		if (not pepseq):
	# 			continue
	# 		if map.source == 'sw_sharp':
	# 			exon_known_code = 2
	# 			hassw = True
	# 		elif map.source == 'usearch':
	# 			exon_known_code = 3
	# 			hassw = True
	# 		else:
	# 			exon_known_code = map.exon_known_2
	# 		seqname = "{0}:{1}:{2}".format(map.species_2, map.exon_id_2, exon_known_code)
	# 		headers.append(seqname)
	# 		sequences[seqname] = pepseq
	# 		# for split exon concatenation (see below)
	# 		if not map.species_2 in list(exons_per_species.keys()):
	# 			exons_per_species[map.species_2] = []
	# 		exons_per_species[map.species_2].append([map.exon_id_2, exon_known_code]);
	#
	# 	if (len(headers) <= 1):
	# 		if verbose: print("single species in the alignment")
	# 		no_orthologues += 1
	# 		continue
	#
	# 	# concatenate exons from the same gene - the alignment program might go wrong otherwise
	#   # for example, for some species the exon is split in two - might be a mistake, but still this is what we have
	# 	concatenated = concatenate_exons(cursor, ensembl_db_name, sequences, exons_per_species)
	#
	# 	fasta_fnm = "{0}/{1}.fa".format(cfg.dir_path['scratch'], human_exon.exon_id)
	# 	output_fasta(fasta_fnm, list(sequences.keys()), sequences)
	#
	# 	# align
	# 	afa_fnm = "{0}/{1}.afa".format(cfg.dir_path['scratch'], human_exon.exon_id)
	# 	mafftcmd = acg.generate_mafft_command(fasta_fnm, afa_fnm)
	# 	ret = subprocess.getoutput(mafftcmd)
	#
	# 	if (verbose): print('almt to', afa_fnm)
	#
	# 	# read in the alignment
	# 	inf = erropen(afa_fnm, "r")
	# 	aligned_seqs = {}
	# 	for record in SeqIO.parse(inf, "fasta"):
	# 		aligned_seqs[record.id] = str(record.seq)
	# 	inf.close()
	# 	# split back the concatenated exons
	# 	if concatenated: split_concatenated_exons(aligned_seqs, concatenated)
	#
	# 	human_seq_seen = False
	# 	for seq_name, sequence in aligned_seqs.items():
	# 		# if this is one of the concatenated seqs, split them back to two
	#
	# 		### store the alignment as bitstring
	# 		# Generate the bitmap
	# 		bs = Bits(bin='0b' + re.sub("[^0]", "1", sequence.replace('-', '0')))
	# 		# The returned value of tobytes() will be padded at the end
	# 		# with between zero and seven 0 bits to make it byte aligned.
	# 		# I will end up with something that looks like extra alignment gaps, that I'll have to return
	# 		msa_bitmap = bs.tobytes()
	# 		# Retrieve information on the cognate
	# 		cognate_species, cognate_exon_id, cognate_exon_known = seq_name.split(':')
	# 		if cognate_exon_known == '2':
	# 			source = 'sw_sharp'
	# 		elif cognate_exon_known == '3':
	# 			source = 'usearch'
	# 		else:
	# 			source = 'ensembl'
	# 		if (cognate_species == 'homo_sapiens'):
	# 			human_seq_seen = True
	# 		cognate_genome_db_id = species2genome_db_id(cursor, cognate_species)  # moves the cursor
	# 		switch_to_db(cursor, ensembl_db_name['homo_sapiens'])  # so move it back to homo sapiens
	# 		# Write the bitmap to the database
	# 		# if (cognate_species == 'homo_sapiens'):
	# 		if verbose:  # and (source=='sw_sharp' or source=='usearch'):
	# 			print("storing")
	# 			print(human_exon.exon_id, human_exon.is_known)
	# 			print(cognate_species, cognate_genome_db_id, cognate_exon_id, cognate_exon_known, source)
	# 			print(sequence)
	# 			if not msa_bitmap:
	# 				print("no msa_bitmap")
	# 				continue
	# 	# store_or_update(cursor, "exon_map",    {"cognate_genome_db_id":cognate_genome_db_id,
	# 	#    "cognate_exon_id":cognate_exon_id   ,"cognate_exon_known"  :cognate_exon_known,
	# 	#    "source": source, "exon_id" :human_exon.exon_id, "exon_known":human_exon.is_known},
	# 	#   {"msa_bitstring":MySQLdb.escape_string(msa_bitmap)})
	#
	# 	ok += 1
	# 	subprocess.getoutput("rm " + afa_fnm + " " + fasta_fnm)
	#


#########################################
def multiple_alignment_gene(cursor, ensembl_db_name, ref_genome_db_id, exon_maps):

	# exon map has exon_id (reference exon) and cognate_exon_id - exon from other species in the taxonomical group
	# here we are only fishing out the sequences
	exon_ids = set([m['exon_id'] for m in exon_maps])
	for exon_id in exon_ids:
		alignable = [m for m in exon_maps if m['exon_id'] == exon_id]
		multiple_alignment_exon(cursor, ensembl_db_name, ref_genome_db_id, alignable)


#########################################
def multiple_alignment_genes(gene_list, other_args):

	print("process pid: %d, length of gene list: %d" % (get_process_id(), len(gene_list)))
	[reference_species, ensembl_db_name] = other_args

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	# I need the following two  to fetch the gene description
	ref_genome_db_id = species2genome_db_id(cursor, reference_species)
	human_db = ensembl_db_name['homo_sapiens']

	t0 = time()
	report_chunk = 100
	column_names = get_column_names(cursor, db_name=ensembl_db_name[reference_species], table_name="exon_map")
	switch_to_db(cursor, ensembl_db_name[reference_species])
	for gene_ct, gene_id in enumerate(gene_list):
		if not gene_ct % report_chunk:
			print(f"{gene_ct} genes out of {len(gene_list)}; the last {report_chunk} in %.2f min" % ((time()-t0)/60))
			t0 = time()
		exon_maps = get_exon_maps(cursor, gene_id, column_names)
		if not exon_maps: continue
		# what is the gene that we are looking at
		print_description(cursor, gene_id, ref_genome_db_id, ensembl_db_name[reference_species], human_db)

		# make multiple alignment
		multiple_alignment_gene(cursor, ensembl_db_name, ref_genome_db_id, exon_maps)

		break

	cursor.close()
	db.close()


#########################################
def main():
	# TODO: cplit codons
	no_processes = 1
	rep_species = 'monodelphis_domestica'

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	[all_species, ensembl_db_name] = get_species(cursor)
	# switch_to_db(cursor,  ensembl_db_name[rep_species])
	# qry = "select distinct(g.gene_id) from exon_map e left join gene2exon g on e.exon_id = g.exon_id"
	# genes_with_maps = [ret[0] for ret in hard_landing_search(cursor, qry)]
	# print(f"genes with maps: {len(genes_with_maps)}")
	# gene_list = sample(genes_with_maps, 10)
	gene_list = [8979] # NARS; only three exons out of 14, and not the same species as the ones for which I have the full seqeuence

	cursor.close()
	db.close()

	parallelize(no_processes, multiple_alignment_genes, gene_list, [rep_species, ensembl_db_name])

	return True


#########################################
if __name__ == '__main__':
	main()

'''
	#for gene_id in [412667]: #  wls
	#for gene_id in [378768]: #  p53
	#for gene_id in [378766]: #  dynein
'''
