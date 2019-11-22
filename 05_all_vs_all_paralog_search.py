#!/usr/bin/python3
# make sure the search db is formatted for blast search
# /usr/bin/makeblastdb -in canonical_peptides.fa -dbtype prot -parse_seqids
import os

from config import Config
from el_utils.ensembl import  *
from el_utils.mysql import  *

####################################################
def process_blast_search(blastout, visited, paralogue_file):
	inf = open(blastout,"r")
	hits = [] # should include query itself
	for line in inf:
		if line[0]=="#": continue
		hit = line.split("\t")[1]
		hits.append(hit)
		visited[hit] = True
	outf = open(paralogue_file,"a")
	hits = list(set(hits))
	outf.write("\t".join(hits)+"\n")

	outf.close()
	inf.close()

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

	#######################
	print("collecting peptide ids")
	stable_translation_ids = []
	visited = {}
	for gene_id in get_gene_ids(cursor, db_name, biotype="protein_coding"):
		canonical_translation_stable_id = gene2canon_transl(cursor, gene_id, stable=True)
		if canonical_translation_stable_id and "ENSP0" in canonical_translation_stable_id:
			stable_translation_ids.append(canonical_translation_stable_id)
			visited[canonical_translation_stable_id] = False

	paralogue_file = "paralogues.txt"
	if  not os.path.exists(paralogue_file): os.remove(paralogue_file)
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
		process_blast_search(blastout, visited, paralogue_file)
		for tmpfile in [blastin, blastout]:
			if os.path.exists(tmpfile): os.remove(tmpfile)
		#exit()

	cursor.close()
	db.close()


#####################################################
if __name__=="__main__":
	main()
