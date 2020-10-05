#!/usr/bin/python3 -u

from random import sample
from   el_utils.el_specific import  *
from el_utils.processes import parallelize
from config import Config

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

verbose = True


####################
# phase: position where the codon is split, not the offset
# phase can be -1 in UTR exons
def find_offset(exon, exon_length):
	phase = max(exon['phase'],0)
	# it seems that the phase always refers to cdna, not the main strand
	offset = (3-phase)%3 # complement of 3, but 0 stays 0
	full_codons = (exon_length-offset)//3*3 # integer division in python3

	return offset, full_codons


#########################################
def reconstruct_exon_seqs (gene_region_dna, sorted_exons, mitochondrial):

	flank_length = Config.exon_flanking_region_length

	exon_seq = {}
	left_flank  = {}
	right_flank = {}
	exon_pepseq = {}
	phase = {}

	# the hope is we have checled some place else
	# fo the bizarre possibilty that not all exons are marked as being on the same strand
	reverse = sorted_exons[0]['strand']<0

	last_coding_exon = sorted_exons[-1]
	for exon in sorted_exons if reverse else reversed(sorted_exons):
		if exon['is_coding']:
			last_coding_exon = exon
			break

	for exon in sorted_exons:
		if not exon['is_coding']: continue
		exon_id = exon['exon_id']
		start = exon['start_in_gene'] if exon['canon_transl_start'] < 0 else exon['canon_transl_start']
		end   = exon['end_in_gene'] if exon['canon_transl_end'] < 0 else exon['canon_transl_end']
		exon_length = end - start + 1

		left   = gene_region_dna[max(start - flank_length, 0):start]
		coding = gene_region_dna[start:end+1]
		right  = gene_region_dna[end+1:end+1 + flank_length]

		if reverse:
			exon_seq[exon_id]    = str(Seq(coding,generic_dna).reverse_complement())
			left_flank[exon_id]  = str(Seq(right,generic_dna).reverse_complement())
			right_flank[exon_id] = str(Seq(left,generic_dna).reverse_complement())
		else:
			exon_seq[exon_id]    = coding
			left_flank[exon_id]  = left
			right_flank[exon_id] = right

		offset, full_codons = find_offset(exon, exon_length)
		if mitochondrial:
			exon_pepseq[exon_id] = str(Seq(exon_seq[exon_id][offset:offset+full_codons],
			                               generic_dna).translate(table="Vertebrate Mitochondrial"))
		else:
			exon_pepseq[exon_id] = str(Seq(exon_seq[exon_id][offset:offset+full_codons], generic_dna).translate())

		# this now is to be understood as the phase in the reading direction of the gene
		phase[exon_id] = (3-offset)%3 # complement of 3, but 0 stays 0

		# sanity
		# it can happen (eg human,  gene id 597950 ENSG00000163239) that an exon contains only a part of a single codon
		# waht we do not want to see, however, is a stop codon in the middle of the seqeunce
		if len(exon_pepseq[exon_id])>0 and \
				('*' in exon_pepseq[exon_id][:-1] or exon!=last_coding_exon and exon_pepseq[exon_id][-1]=='*'):
			return f"Error: stop codon while trying to recosntruct exon id {exon_id}."
			# print(f"problem when translating exon {sorted_exons.index(exon)} {len(sorted_exons)}")
			# print(f"mitochondrial: {mitochondrial}, strand: {exon['strand']},  phase reported: {exon['phase']}, phase calculated: { phase[exon_id]}")
			# print(f"  exon length: {exon_length}, exon_length%3: {exon_length%3}")
			# print(exon_pepseq[exon_id])
			# for offset in range(3):
			# 	print(Seq(exon_seq[exon_id][offset:offset+full_codons], generic_dna).translate())
			# print()
			# exit()

	return [exon_seq, left_flank, right_flank, exon_pepseq]


#########################################
def store_exon_seqs(cursor, exons, exon_seq, left_flank, right_flank, exon_pepseq):

	for exon in exons:
		if not exon['is_coding']: continue
		exon_id = exon['exon_id']
		#####
		fixed_fields  = {}
		update_fields = {}
		fixed_fields['exon_id'] = exon_id
		fixed_fields['by_exolocator']   = 0
		update_fields['phase']  = exon['phase']
		if exon_seq[exon_id]:
			update_fields['dna_seq']     = exon_seq[exon_id]
		if left_flank[exon_id]:
			update_fields['left_flank']  = left_flank[exon_id]
		if right_flank[exon_id]:
			update_fields['right_flank'] = right_flank[exon_id]
		if exon_id in exon_pepseq and exon_pepseq[exon_id]:
			update_fields['protein_seq'] = exon_pepseq[exon_id]
		store_or_update (cursor, 'exon_seq', fixed_fields, update_fields, primary_key='exon_seq_id')


def store_exon_seqs_gene(cursor, gene_id, species, ensembl_db_name):
	db_name = ensembl_db_name[species]
	# print(f" {gene_id} {gene2stable(cursor, gene_id, db_name=db_name)}")
	# extract raw gene  region - bonus return from checking whether the
	# sequence is correct: translation of canonical exons
	gene_region_dna = get_gene_dna(cursor, species, db_name, gene_id)
	if not gene_region_dna:
		return f"Error: gene equence not found for {gene_id} {gene2stable(cursor, gene_id, db_name=db_name)}"

	# get _all_ exons
	ret = get_sorted_canonical_exons(cursor, db_name, gene_id)
	if not ret:
		return f"Error: no exons found."
	if type(ret)==str and 'Error' in ret: # some sort of failure msg - ignore warnings
		return ret
	# if everything went ok, we have the sorted list of exons
	sorted_exons = ret

	# cdna - use this to chek if it is translateble into protein
	cdna = exons2cdna(gene_region_dna, sorted_exons)
	# sanity
	mitochondrial = is_mitochondrial(cursor, gene_id, db_name=db_name)
	if mitochondrial:
		protein = Seq(cdna, generic_dna).translate(table="Vertebrate Mitochondrial")
	else:
		protein = Seq(cdna, generic_dna).translate()
	if '*' in protein[:-1]:
		return f"Error: stop codon in the protein sequence."

	# (the return are four dictionaries, with exon_ids as keys)
	ret = reconstruct_exon_seqs(gene_region_dna, sorted_exons, mitochondrial)
	if type(ret)==str: return ret # there was a problem
	[exon_seq, left_flank, right_flank, exon_pepseq] = ret
	store_exon_seqs(cursor, sorted_exons, exon_seq, left_flank, right_flank, exon_pepseq)
	return "ok"

#########################################
def store_exon_seqs_species(species_list, other_args):

	[ensembl_db_name] = other_args
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")

	for species in species_list:
		print()
		print("############################")
		print(species)
		switch_to_db(cursor, ensembl_db_name[species])
		gene_ids = get_gene_ids (cursor, biotype='protein_coding')
		print(f"nummber of genes: {len(gene_ids)}")

		###########################
		tot = 0
		fail_ct = 0
		# gene_ids = [596919] # this is ABCA4
		# gene_ids = sample(gene_ids, 10) # testing
		for gene_id in gene_ids:
			tot += 1
			if tot%1000 == 0: print(species, "tot genes:", tot, " fail:", fail_ct)
			ret = store_exon_seqs_gene(cursor, gene_id, species, ensembl_db_name)
			if ret != "ok":
				# print(gene_id, ret)
				fail_ct += 1
				store_problem(cursor, ensembl_db_name[species], gene_id, f"When storing exon seqs: {ret}")

		print(species, "done; tot:", tot, " fail:", fail_ct)

	cursor.close()
	db    .close()


#########################################
def main():

	"""
	Main entry point, but in reality does nothing except taking care of the parallelization.
	The parallelization here is per-species.
	"""

	no_threads = 8

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species (cursor)
	#all_species = ['homo_sapiens']
	#all_species = sample(all_species, 10)

	cursor.close()
	db    .close()

	parallelize (no_threads, store_exon_seqs_species, all_species, [ensembl_db_name])



#########################################
if __name__ == '__main__':
	main()
'''
	# Ivana:
	# What is the proper way to find out whether the seq_region is mitochondrial?
	# EMily:
	# It's in the seq_region table, under name. Name will either be a chromosome
	# number, MT for mitochrondria or the name of a contig.
Hi Ivana

To find which contigs are mitochondrial, use the assembly table.

select s1.name from seq_region s1, assembly asm, seq_region s2 where
s1.seq_region_id = asm.cmp_seq_region_id and s2.seq_region_id =
asm.asm_seq_region_id and s2.name = "MT" ;
in human, returns
+----------------+
| name |
+----------------+
| NC_012920 |
| J01415.1.16569 |
+----------------+

Emily

Hi Ivana

The reason for this is that patches are not always the same length as the
genomic region they patch over. In most cases, a patch corrects sequencing
errors but the number of bp stays the same, but in some cases, a patch adds a
chunk of sequence. We anchor the 5' (wrt the chromosome orientation) end of the
patch to the identical reference coordinates, and allow the 3' end to differ
slightly from the reference coordinates. We don't change the complete genomic
coordinates every time a patch is added (this happens when we bring out a new
assembly) as this would be more hassle than it's worth and most people don't
notice it anyway. However the one thing it does affect is the coordinates of
genes that overlap the 3' end of the patch.

ENSG00000261899, and presumably the other genes you've had a similar issue
with, overlaps the 3' end of a patch, so its coordinates on the patch are
shifted compared to its coordinates on the reference genome.

Here's the data from the assembly_exception table for the particular patch that
affects ENSG00000261899:
+-----------------------+---------------+------------------+----------------+-------------+------------------+-----+
| assembly_exception_id | seq_region_id | seq_region_start | seq_region_end |exc_type | exc_seq_region_id | exc_seq_region_start | exc_seq_region_end | ori
|
+-----------------------+---------------+------------------+----------------+-------------+-------------------+----------------------+--------------------+-----+
| 67                    | 1000057054    | 36453102 | 36596491 | PATCH_NOVEL | 27508 | 36453102 |36590458 | 1 |
+-----------------------+---------------+------------------+----------------+-------------+-------------------+----------------------+--------------------+-----+

To get over this you need to link over the assembly_exception table. This will
give you the exact relationship between the patch you are looking at and the
chromosome it is linked to. If you then use the
seq_region_start/exc_seq_region_start relationship, this will cover all cases,
whether there is a shift or not.

Hope this helps,

Emily
'''
