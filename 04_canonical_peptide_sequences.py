#!/usr/bin/python3 -u


from config import Config

from el_utils.ensembl import  *
from el_utils.processes import *


####################################################
def prioritize_versioning(fastasfile):
	headerct = 0
	index_versioning = 0
	with open(fastasfile) as inf:
		for line in inf:
			if not '>' in line: continue
			headerct +=1
			name = line.split()[0]
			if "." in name:
				index_versioning += 1
			if headerct==100: break
	if index_versioning/headerct>0.5:
		return [f".{v}" for v in range(1,13)] + [""]
	else:
		return [""] + [f".{v}" for v in range(1,13)]


###########
def canon_seqs_for_species(cursor, species, ensembl_db_name):

	print(f"#######################\n{species}")

	# we'll do this for human only, for now
	switch_to_db (cursor, ensembl_db_name[species])

	fasta_path = "/storage/databases/ensembl-{}/fasta/{}/pep".format(Config.release_number, species)
	ext = "pep.all.fa"
	in_fastas = [fnm for fnm in  os.listdir(fasta_path) if fnm[-(len(ext)):] == ext]
	if len(in_fastas) == 0:
		print(f"*.{ext} not found")
		exit()
	if len(in_fastas) > 1:
		print("multiple *.pep.all.fa  found:", in_fastas)
		exit()
	in_fasta = in_fastas[0]
	# this is very unsystematic, but some species such as Mus spretus do not have the version number
	# so in that case look for that name format first
	versioning = prioritize_versioning(f"{fasta_path}/{in_fasta}")

	out_fasta = "{}/canonical_peptides.fa".format(fasta_path)

	#######################
	print(f"{species} collecting peptide ids")
	time0 = time()
	stable_translation_ids = []
	# get_gene_ids(cursor, db_name=None, biotype=None,  stable=False, ref_only=False):
	for gene_id in get_gene_ids(cursor, ensembl_db_name[species], biotype="protein_coding"):
		canonical_translation_stable_id = gene2stable_canon_transl_id(cursor, gene_id)
		if canonical_translation_stable_id: stable_translation_ids.append(canonical_translation_stable_id)
	mins = (time()-time0)/60
	print(f"{species} collecting peptide ids done in %.1f min" % mins)


	#######################
	print(f"{species} extracting sequences")
	time0 = time()
	tmpfile = f"tmp.{species}.fa"


	count = 0
	for stable_translation_id in stable_translation_ids:
		count += 1
		if count%1000 == 0: print(f"{species} {count} of {len(stable_translation_ids)}")

		if os.path.exists(tmpfile): os.remove(tmpfile)

		logfile =  f"blastcmd.{species}.log"
		for version in versioning:
			# keep only the last log
			# %s format (as a blast option) is sequence without the header
			cmd  = f"{Config.blastdbcmd} -db {fasta_path}/{in_fasta} -dbtype prot -out {tmpfile} -outfmt %s "
			cmd += f"-entry {stable_translation_id}{version} -logfile {logfile}"
			subprocess.call(["bash","-c", cmd])
			if os.path.exists(tmpfile) and os.path.getsize(tmpfile)>0:
				cmd = "echo '>{}' >> {}".format(stable_translation_id, out_fasta)
				subprocess.call(["bash","-c", cmd])
				cmd = "cat {} >> {}".format(tmpfile, out_fasta)
				subprocess.call(["bash","-c", cmd])
				break
	# blast also produces some emoty junkfile with extension perf
	perfile =  f"blastcmd.{species}.perf"
	if os.path.exists(perfile): os.remove(perfile)
	mins = (time()-time0)/60
	print("%s extracting seqs  done in %.1f min" % (species, mins))

	if os.path.exists(tmpfile): os.remove(tmpfile)



###########
def canon_seqs(species_chunk, other_args):
	[ensembl_db_name] = other_args
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	for species in species_chunk:
		canon_seqs_for_species(cursor, species, ensembl_db_name)
	cursor.close()
	db.close()



####################################################
def main():
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	all_species.remove("homo_sapiens")
	cursor.close()
	db.close()

	number_of_chunks = 8

	print(f"number of species: {len(all_species)}, number of chunks: {number_of_chunks}")
	parallelize(number_of_chunks, canon_seqs, all_species, [ensembl_db_name])



#####################################################
if __name__=="__main__":
	main()
