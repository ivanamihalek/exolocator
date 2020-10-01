#!/usr/bin/python3
# make sure the search db is formatted for blast search
# /usr/bin/makeblastdb -in canonical_peptides.fa -dbtype prot -parse_seqids
import os

from config import Config
from el_utils.ensembl import  *
from el_utils.mysql import  *
from el_utils.processes import *
# after this, to proceed to gene level paralogues need to know exon structure
# so do the exons first

####################################################
def process_blast_search(stable_translation_id, blastout, paralogue_file):
	inf = open(blastout,"r")
	hits = [stable_translation_id] # should include query itself
	for line in inf:
		if line[0]=="#": continue
		hit = line.split("\t")[1]
		if hit not in hits: hits.append(hit)

	outf = open(paralogue_file,"a")
	outf.write(f"{stable_translation_id}\t{','.join(hits[1:])}\n")

	outf.close()
	inf.close()

####################
def find_paralogues(stable_translation_ids, other_args):
	[in_fasta, out_dir] = other_args
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	paralogue_file = f"{out_dir}/paralogues.{os.getpid()}.tsv"
	cmd = "touch {}".format(paralogue_file)
	subprocess.call(["bash","-c", cmd])

	for stable_translation_id in stable_translation_ids:
		#if visited[stable_translation_id]: continue
		#########
		# show some signs of life
		if stable_translation_ids.index(stable_translation_id)%100==0:
			print("###", stable_translation_ids.index(stable_translation_id), "out of", len(stable_translation_ids))

		#########
		# extract the query sequence
		blastin = "{}.fa".format(stable_translation_id)
		# %f is sequence with fasta header
		cmd = "{} -db {} -dbtype prot -out {} -outfmt %f -entry {} -logfile blastcmd.log".\
				format(Config.blastdbcmd, in_fasta, blastin,  stable_translation_id)
		subprocess.call(["bash","-c", cmd])
		if not os.path.exists(blastin) or os.path.getsize(blastin)==0: continue

		#########
		# the actual search
		blastout = "{}.blastp".format(stable_translation_id)
		if  os.path.exists(blastout): os.remove(blastout)
		# 6 is tabular foramt, 7 tabular with comment lines
		cmd = "{} -db {} -out {} -outfmt 7 -query {} -evalue 0.00001".\
				format(Config.blastp, in_fasta, blastout,  blastin)
		subprocess.call(["bash","-c", cmd])

		#########
		# store paralogues, mark them as visited
		process_blast_search(stable_translation_id, blastout, paralogue_file)
		for tmpfile in [blastin, blastout]:
			if os.path.exists(tmpfile): os.remove(tmpfile)

	cursor.close()
	db.close()


####################################################
def main():
	# we'll do this for human only, for now
	species = 'homo_sapiens'

	fasta_path = "/storage/databases/ensembl-{}/fasta/{}/pep".format(Config.release_number, species)
	in_fasta = "{}/canonical_peptides.fa".format(fasta_path)
	for dep in [fasta_path, in_fasta, Config.blastdbcmd]:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	db_name = ensembl_db_name[species]
	switch_to_db (cursor, ensembl_db_name[species])

	out_dir = "raw_tables"
	if not os.path.exists(out_dir): os.mkdir(out_dir)
	for fnm in os.listdir(out_dir):
		if "paralogues" not in fnm: continue
		os.remove(f"{out_dir}/{fnm}")

	#######################
	print("collecting peptide ids")
	stable_translation_ids = []
	visited = {}
	for gene_id in get_gene_ids(cursor, db_name, biotype="protein_coding"):
		canonical_translation_stable_id = gene2canon_transl(cursor, gene_id, stable=True)
		if canonical_translation_stable_id and "ENSP0" in canonical_translation_stable_id:
			stable_translation_ids.append(canonical_translation_stable_id)
			visited[canonical_translation_stable_id] = False
	cursor.close()
	db.close()

	#######################
	number_of_chunks = 8
	parallelize(number_of_chunks, find_paralogues, stable_translation_ids, [in_fasta, out_dir])

#####################################################
if __name__=="__main__":
	main()
