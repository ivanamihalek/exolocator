#!/usr/bin/python3

# somebody somehwere wrote a buggy program and its output became a norm:
# sometimes, but no systematically, for a gene coded in reverse strand
# when mutations are within the exons, their position is numbered from left to right

# c.3607G>T means that the position 3607 on the coding DNA (counted from the left) is mutated to T
# however, c.3523-2A>G means that at position  3523 counted from *right* there is mutation in the intron,
# 2 positions to the right(?) whichever imbecile came up with that

import os
import subprocess

from el_utils.el_specific import *
from el_utils.ensembl import *
from config import Config
from Bio.Seq  import Seq
from Bio.Alphabet import generic_dna, generic_protein
import pprint


def exon_flanks():
	return

##########
def get_seq_data(cursor, species, db_name, gene_id, sorted_exons):

	splice_length = Config.exon_flanking_region_length

	gene_region_dna = get_gene_dna(cursor, species, db_name, gene_id)

	switch_to_db(cursor, db_name)
	qry =  f"select seq_region_id, seq_region_start,  seq_region_strand from gene where gene_id={gene_id}"
	[gene_region_id, gene_region_start,  gene_region_strand] = hard_landing_search(cursor, qry)[0]
	reverse = gene_region_strand<0

	qry =  f"select name from seq_region where seq_region_id={gene_region_id}"
	seq_region_name = hard_landing_search(cursor, qry)[0][0]

	cdna = exons2cdna(gene_region_dna, sorted_exons)
	protein = str(Seq(cdna, generic_dna).translate())
	cdna_length = len(cdna)

	print(f"cdna_length = {cdna_length}, cdna_length%3 = {cdna_length%3}, cdna_length/3 = {cdna_length/3}")
	print(protein)
	print()

	acceptor_splice = {}
	donor_splice = {}
	cdna2gdna = {}
	cumulative_length = 0

	for exon in reversed(sorted_exons) if reverse else sorted_exons:

		# these coordinates are forward
		start = exon['start_in_gene'] if exon['canon_transl_start'] < 0 else exon['canon_transl_start']
		end   = exon['end_in_gene'] if exon['canon_transl_end'] < 0 else exon['canon_transl_end']
		exon_length = end - start + 1

		# exon start and end on the cdna sequence - this one runs backward on reverse strand
		cdna2gdna[cumulative_length+1] = gene_region_start + (end if reverse else start)
		cdna2gdna[cumulative_length+exon_length] = gene_region_start + (start if reverse else end)
		# store the splice site fragments
		if reverse:
			flank = Seq(gene_region_dna[end+1:end + splice_length], generic_dna)
			acceptor_splice[cumulative_length+1] = str(flank.reverse_complement())
			flank = Seq( gene_region_dna[max(start-splice_length, 0):start], generic_dna)
			donor_splice[cumulative_length+exon_length] = str(flank.reverse_complement())
		else:
			acceptor_splice[cumulative_length+1] = gene_region_dna[max(start-splice_length, 0):start]
			donor_splice[cumulative_length+exon_length] = gene_region_dna[end+1:end + splice_length]
		# how far did we move along the cDNA?
		cumulative_length += exon_length


	# for cdnapos in sorted(cdna2gdna.keys()):
	# 	print(cdnapos, cdna2gdna[cdnapos])
	# 	if cdnapos in acceptor_splice:
	# 		print(f"\tacceptor splice: {acceptor_splice[cdnapos]}")
		# if cdnapos in donor_splice:
		# 	print(f"\tdonor splice: {donor_splice[cdnapos]}")
	return seq_region_name, cdna, protein, acceptor_splice, donor_splice, cdna2gdna


output_format = '''

SYMBOL_chromosome = "{}"

SYMBOL_cdna = \'\'\'
{}
\'\'\'

SYMBOL_protein =  \'\'\'
{}
\'\'\'

SYMBOL_donor_splice = {}

SYMBOL_acceptor_splice = {}

SYMBOL_cdna2gdna = {}

def get_cdna():
	return SYMBOL_cdna.replace("\\n", "")
	
def get_protein():
	return SYMBOL_protein.replace("\\n", "")
	
def get_codons():
	return [SYMBOL_cdna[i:i+3] for i in range(0,len(SYMBOL_cdna),3)]
	
def get_backward_numbered_acceptor_splice(bdry_pos):
	cumulative_length = max(SYMBOL_acceptor_splice.keys())
	return SYMBOL_acceptor_splice.get(cumulative_length-bdry_pos, None)
	
def get_backward_numbered_donor_splice(bdry_pos):
	cumulative_length = max(SYMBOL_donor_splice.keys())
	return SYMBOL_donor_splice.get(cumulative_length-bdry_pos, None)

'''

########################################
def main():

	global output_format

	symbol = "ABCA4"
	species = "homo_sapiens"

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)

	# ensmebl id from gene name
	qry = f"select ensembl_gene_id from identifier_maps.hgnc where approved_symbol='{symbol}'"
	stable_id = hard_landing_search(cursor, qry)[0][0]
	gene_id = stable2gene(cursor, stable_id,  ensembl_db_name[species])
	tr_id = get_canonical_transcript_id(cursor, gene_id, ensembl_db_name[species])
	qry=f"select stable_id from {ensembl_db_name[species]}.transcript where transcript_id={tr_id}"
	canonical_tr_stable_id = hard_landing_search(cursor,qry)[0][0]
	print(symbol, stable_id, gene_id, canonical_tr_stable_id)

	# coordinates of canonical exons
	# we will just extract positions of exon boundaries on cdna
	sorted_exons = get_sorted_canonical_exons(cursor, ensembl_db_name[species], gene_id)
	if type(sorted_exons)==str:
		print(sorted_exons) # some sort of failure; this is the errmsg
		exit()

	#
	# canonical cdna and its translation
	ret = get_seq_data(cursor, species, ensembl_db_name[species], gene_id, sorted_exons)
	[seq_region_name, cdna, protein, acceptor_splice, donor_splice, cdna2gdna] = ret
	pp = pprint.PrettyPrinter(indent=4)
	with open(f"{symbol.lower()}_gene.py", "w") as outf:
		newlined_cdna = "\n".join([cdna[i:i+100] for i in range(0,len(cdna),100)])
		newlined_protein = "\n".join([protein[i:i+100] for i in range(0,len(protein),100)])
		output_format =  output_format.replace("SYMBOL", symbol.lower())
		outf.write(output_format.format(seq_region_name, newlined_cdna, newlined_protein, pp.pformat(donor_splice),
										pp.pformat(acceptor_splice),  pp.pformat(cdna2gdna)))

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

