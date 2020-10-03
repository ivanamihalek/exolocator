#!/usr/bin/python3

# somebody somehwere wrote a buggy program and its output became a norm:
# sometimes, but no systematically, for a gene coded in reverse strand
# when mutations are within the exons, their position is numbered from left to right

# c.3607G>T means that the position 3607 on the coding DNA (counted from the left) is mutated to T
# however, c.3523-2A>G means that at position  3523 counted from *right* there is mutation in the intron,
# 2 positions to the right(?) whichever imbecile came up with that

import os
import subprocess

from el_utils.mysql import *
from el_utils.ensembl import *
from config import Config
from Bio.Seq  import Seq
from Bio.Alphabet import generic_dna, generic_protein
import pprint

# class Exon:
# 	def __init__ (self, exon_id, seq_region_start, seq_region_end):
# 		self.id = exon_id
# 		self.start = seq_region_start
# 		self.end = seq_region_end
# 		self.is_first = False
# 		self.is_last = False
# 		self.donor_splice = None
# 		self.acceptor_splice = None
# 	def to_string(self):
# 		return f"id:{self.id}   start:{self.start}   end:{self.end}"
#

def get_canonical_exon_bdries(cursor, db_name, gene_id):

	column_names = get_column_names(cursor, db_name, "gene2exon")
	exons = []
	for ret in hard_landing_search(cursor, f"select * from {db_name}.gene2exon where gene_id={gene_id}"):
		exon = dict(zip(column_names, ret))
		if not exon['is_canonical']: continue
		exons.append(exon)

	sorted_exons = sorted(exons,key=lambda x: x['start_in_gene'])
	# sanity check
	total_length = 0
	for exon in sorted_exons:
		#print(exon)
		reading = exon['is_coding']
		if reading:
			start = exon['start_in_gene'] if exon['canon_transl_start']<0 else  exon['canon_transl_start']
			end = exon['end_in_gene'] if exon['canon_transl_end']<0 else  exon['canon_transl_end']
			total_length += end - start  + 1

	if total_length%3 != 0:
		print(f"total length {total_length} not divisible by 3")
		exit()
	print(f"total_length = {total_length}, total_length%3 = {total_length%3}, total_length/3 = {total_length/3}")
	print()

	return sorted_exons


def get_file_name(cursor, seq_region_id):
	# for human I do not expect that a gene would be split across different regions,
	# but for other species not so sure
	qry =  f"select file_ids from seq_region2file where seq_region_id={seq_region_id}"
	file_ids = hard_landing_search(cursor,qry)[0][0]
	qry = f"select file_name from file_names where file_id in ({file_ids})"
	file_names = [ret[0] for ret in hard_landing_search(cursor, qry)]

	if len(file_names)>1:
		print(f"multiple file names found for seq_region_id {seq_region_id}")
		exit()
	return file_names[0]


#########################################
def extract_cdna_from_fasta(fasta_full_path, seq_region_name, seq_region_start,seq_region_end):
	tmpfile = f"tmp.{os.getpid()}.fasta"
	if os.path.exists(tmpfile): os.remove(tmpfile)
	cmd  = f"{Config.blastdbcmd} -db {fasta_full_path} -dbtype nucl -entry {seq_region_name} "
	cmd += f"-range {seq_region_start}-{seq_region_end} -out {tmpfile} -outfmt %s"
	subprocess.call(["bash", "-c", cmd])
	if not os.path.exists(tmpfile):
		print(f"{tmpfile} not produced")
		exit()
	with open(tmpfile) as inf:
		inseq = inf.read().replace("\n", "")
	os.remove(tmpfile)
	return inseq


##########
def get_cdna(cursor, species, db_name, gene_id, sorted_exons):

	switch_to_db(cursor, db_name)

	qry =  f"select seq_region_id, seq_region_start, seq_region_end, seq_region_strand from gene where gene_id={gene_id}"
	[seq_region_id, seq_region_start, seq_region_end, seq_region_strand] = hard_landing_search(cursor, qry)[0]

	qry =  f"select name from seq_region where seq_region_id={seq_region_id}"
	seq_region_name = hard_landing_search(cursor, qry)[0][0]

	file_name = get_file_name(cursor, seq_region_id)
	print(seq_region_id, seq_region_start, seq_region_end, seq_region_strand, seq_region_name, file_name)

	# convention over configuration :}
	fasta_full_path = f"{Config.fasta_repo}/{species}/dna/{file_name}"
	if not os.path.exists(fasta_full_path):
		print(f"{fasta_full_path} not found")
		exit()
	cdna = extract_cdna_from_fasta(fasta_full_path, seq_region_name, seq_region_start,seq_region_end)
	print(cdna[:50])
	flank_length = Config.exon_flanking_region_length
	reverse = seq_region_strand<0
	#
	# seq = ""
	# for exon in exons:
	# 	with open(tmpfile) as inf:
	# 		inseq = inf.read().replace("\n", "")
	# 	seq += inseq[splice_length:-splice_length]
	# 	if reverse:
	# 		exon.donor_splice    = str(Seq(inseq[:splice_length]).reverse_complement())
	# 		exon.acceptor_splice = str(Seq(inseq[-splice_length:]).reverse_complement())
	# 	else:
	# 		# if the exon is not the first, technically this is not the splice region
	# 		exon.acceptor_splice = inseq[:splice_length]
	# 		exon.donor_splice    = inseq[-splice_length:]
	# 	# I do not know wtf this is, but it correposnds to the way variants are reported in the literature
	#
	# acceptor_splice = {}
	# donor_splice = {}
	# cdna2gdna = {}
	# cumulative_length = 0
	# exon_boundary = 0
	#
	# for exon in sorted(exons, key=lambda e: e.start, reverse=reverse):
	# 	acceptor_splice[exon_boundary+1] = exon.acceptor_splice
	# 	cdna2gdna[exon_boundary+1] = exon.end if reverse else exon.start
	# 	cumulative_length += exon.end-exon.start+1
	#
	# 	exon_boundary  = cumulative_length
	#
	# 	donor_splice[exon_boundary] = exon.donor_splice
	# 	cdna2gdna[exon_boundary] = exon.start if reverse else exon.end
	#
	# # for eb in sorted(cdna2gdna.keys()):
	# # 	print(f"{eb}    {cdna2gdna[eb]}")
	# # exit()
	# print(cumulative_length, cumulative_length%3, cumulative_length/3)
	# biopython_dna = Seq(seq, generic_dna)
	# if reverse: biopython_dna = biopython_dna.reverse_complement()
	# cdna = str(biopython_dna)
	# protein = str(biopython_dna.translate())
	# # print(biopython_dna)
	# # print("============")
	# # print(protein)
	# # codon = [cdna[i:i+3] for i in range(0,len(cdna),3)]
	# # for i in range(len(protein)):
	# # 	print(i, protein[i], codon[i])
	# # exit()
	# return seq_region_name, cdna, protein, acceptor_splice, donor_splice, cdna2gdna
	#
	return [0]*6


output_format = '''

abca4_chromosome = "{}"

abca4_cdna = \'\'\'
{}
\'\'\'

abca4_protein =  \'\'\'
{}
\'\'\'

abca4_donor_splice = {}

abca4_acceptor_splice = {}

abca4_cdna2gdna = {}

def get_cdna():
	return abca4_cdna.replace("\\n", "")
	
def get_protein():
	return abca4_protein.replace("\\n", "")
	
def get_codons():
	return [abca4_cdna[i:i+3] for i in range(0,len(abca4_cdna),3)]
	
def get_backward_numbered_acceptor_splice(bdry_pos):
	cumulative_length = max(abca4_acceptor_splice.keys())
	return abca4_acceptor_splice.get(cumulative_length-bdry_pos, None)
	
def get_backward_numbered_donor_splice(bdry_pos):
	cumulative_length = max(abca4_donor_splice.keys())
	return abca4_donor_splice.get(cumulative_length-bdry_pos, None)

'''

########################################
def main():

	symbol = "ABCA4"
	species = "homo_sapiens"

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)

	# ensmebl id from gene name
	qry = f"select ensembl_gene_id from identifier_maps.hgnc where approved_symbol='{symbol}'"
	stable_id = hard_landing_search(cursor, qry)[0][0]
	gene_id = stable2gene(cursor, stable_id,  ensembl_db_name[species])
	print(symbol, stable_id, gene_id)

	# coordinates of canonical exons
	# we will just extract positions of exon boundaries on cdna
	sorted_exons = get_canonical_exon_bdries(cursor, ensembl_db_name[species], gene_id)
	#
	# canonical cdna and its translation
	ret = get_cdna(cursor, species, ensembl_db_name[species], gene_id, sorted_exons)
	[seq_region_name, cdna, protein, acceptor_splice, donor_splice, cdna2gdna] = ret
	# pp = pprint.PrettyPrinter(indent=4)
	# with open(f"{gene.lower()}_gene.py", "w") as outf:
	# 	newlined_cdna = "\n".join([cdna[i:i+100] for i in range(0,len(cdna),100)])
	# 	newlined_protein = "\n".join([protein[i:i+100] for i in range(0,len(protein),100)])
	# 	outf.write(output_format.format(seq_region_name, newlined_cdna, newlined_protein, pp.pformat(donor_splice),
	# 									pp.pformat(acceptor_splice),  pp.pformat(cdna2gdna)))
	#
	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

