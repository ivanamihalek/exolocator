#!/usr/bin/python3

import subprocess, os

from config import Config

from el_utils.ensembl import  *


####################################################


####################################################
def main():
	# we'll do this for human only, for now
	species = 'homo_sapiens'

	fasta_path = "/storage/databases/ensembl-{}/fasta/{}/pep".format(Config.release_number, species)
	in_fasta = "{}/Homo_sapiens.GRCh38.pep.all.fa".format(fasta_path)
	for dep in [fasta_path, in_fasta, Config.blastdbcmd]:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()

	out_fasta = "{}/canonical_peptides.fa".format(fasta_path)

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	db_name = ensembl_db_name[species]
	switch_to_db (cursor, ensembl_db_name[species])

	#######################
	print("collecting peptide ids")
	stable_translation_ids = []
	for gene_id in get_gene_ids(cursor, db_name, biotype="protein_coding"):
		canonical_translation_stable_id = gene2canon_transl(cursor, gene_id, stable=True)
		if canonical_translation_stable_id and "ENSP0" in canonical_translation_stable_id:
			stable_translation_ids.append(canonical_translation_stable_id)


	#######################
	print("extracting sequences")

	tmpfile = "tmp.fa"

	for stable_translation_id in stable_translation_ids:
		version = 1
		seq_not_found = True
		if os.path.exists(tmpfile): os.remove(tmpfile)
		while seq_not_found and version<12:
			# %s is sequence without the header
			cmd = "{} -db {} -dbtype prot -out {} -outfmt %s -entry {}.{} -logfile balstcmd.log".\
				format(Config.blastdbcmd, in_fasta, tmpfile,  stable_translation_id, version)
			subprocess.call(["bash","-c", cmd])
			if os.path.exists(tmpfile) and os.path.getsize(tmpfile)>0:
				seq_not_found = False
				cmd = "echo '>{}' >> {}".format(stable_translation_id, out_fasta)
				subprocess.call(["bash","-c", cmd])
				cmd = "cat {} >> {}".format(tmpfile, out_fasta)
				subprocess.call(["bash","-c", cmd])
			else:
				version += 1

	if os.path.exists(tmpfile): os.remove(tmpfile)

	cursor.close()
	db.close()


#####################################################
if __name__=="__main__":
	main()
