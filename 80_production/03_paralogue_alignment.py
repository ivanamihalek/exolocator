#!/usr/bin/python3 -u
from el_utils.mysql import  *
from el_utils.ensembl import *
from config import Config
from el_utils.tree import species_sort


def write_to_fasta(home, species, stable_transl_id, tmpfile, logfile, out_fasta):

	fasta_path = "/storage/databases/ensembl-{}/fasta/{}/pep".format(Config.release_number, species)
	assert os.path.exists(fasta_path), f"{fasta_path} not found"
	os.chdir(fasta_path)

	canonical_fasta = f"canonical_peptides.fa"
	assert os.path.exists(canonical_fasta), f"{canonical_fasta} not found"
	cmd  = f"{Config.blastdbcmd} -db {canonical_fasta} -dbtype prot -out {home}/{tmpfile} -outfmt %s "
	cmd += f"-entry {stable_transl_id} -logfile {logfile}"
	subprocess.call(["bash","-c", cmd])
	if os.path.exists(logfile.replace(".log",".perf")): os.remove(logfile.replace(".log",".perf"))
	os.chdir(home)

	if os.path.exists(tmpfile) and os.path.getsize(tmpfile)>0:
		name = stable_transl_id
		cmd = "echo '>{}' >> {}".format(name, out_fasta)
		subprocess.call(["bash","-c", cmd])
		cmd = "cat {} >> {}".format(tmpfile, out_fasta)
		subprocess.call(["bash","-c", cmd])
		return True

	return False


def reorder_seqs(in_afa, ids_sorted, hgnc_symbol, out_afa):
	seq = {}
	name = 'anon'
	with open(in_afa) as inf:
		for line in inf:
			line = line.strip()
			if len(line)==0:continue
			if line[0]=='>':
				name = line[1:]
				seq[name] = ""
			else:
				seq[name] += line

	with open(out_afa, "w") as outf:
		for name in ids_sorted:
			outf.write(f">{hgnc_symbol[name]}\n")
			for i in range(0,len(seq[name]), 100):
				outf.write(f"{seq[name][i:i+100]}\n")


#########################################
def main():

	if len(sys.argv)<2:
		print("usage: %s <gene symbol> " % sys.argv[0])
		exit()
	gene_name = sys.argv[1]
	ref_species = 'homo_sapiens'  # the paralog table is filled only here, for the moment

	out_fasta = f"{gene_name}.paras.fasta"
	out_afa = f"{gene_name}.paras.afa"
	tmpfile = "tmp.fa"
	logfile = "tmp.log"
	for fnm in [out_fasta, out_afa, tmpfile, logfile]:
		if os.path.exists(fnm): os.remove(fnm)
	home = os.getcwd()

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	qry = "select ensembl_gene_id  from identifier_maps.hgnc where approved_symbol='%s'" % gene_name
	ensembl_stable_gene_id = hard_landing_search(cursor, qry)[0][0]

	[all_species, ensembl_db_name] = get_species(cursor)

	switch_to_db(cursor, ensembl_db_name[ref_species])
	qry = "select gene_id from gene where stable_id='%s'" % ensembl_stable_gene_id
	gene_id = hard_landing_search(cursor, qry)[0][0]

	ref_stable_transl_id = gene2stable_canon_transl_id(cursor, gene_id,  ensembl_db_name[ref_species])
	print(gene_name, ref_stable_transl_id)
	write_to_fasta(home, ref_species, ref_stable_transl_id, tmpfile, logfile, out_fasta)
	paras_in_the_almt = []
	qry = f"select stable_ids from paralogue_groups  where bait_stable_id='{ref_stable_transl_id}'"
	ret = error_intolerant_search(cursor, qry)
	if not ret:
		print(f"no paralogues found for {ref_stable_transl_id}")
	else:
		hgnc_symbol = {ref_stable_transl_id:gene_name}
		stable_ids_sorted = ret[0][0].split(",") # these should be in order of similarity to bait (see 07_all_vs_all_paralogue_search.py)
		for stable_transl_id in stable_ids_sorted:
			hgnc_symbol[stable_transl_id] = canonical_protein2hgnc_symbol(cursor, stable_transl_id)

			ok = write_to_fasta(home, ref_species, stable_transl_id, tmpfile, logfile, out_fasta)
			if ok: paras_in_the_almt.append(stable_transl_id)

		if os.path.exists(tmpfile): os.remove(tmpfile)

		cmd  = f"{Config.muscle} -in {out_fasta} -out tmp.afa"
		subprocess.call(["bash","-c", cmd])

		# this will also rename them by their HGNC symbol
		reorder_seqs('tmp.afa', [ref_stable_transl_id] + stable_ids_sorted, hgnc_symbol, out_afa)
		if os.path.exists('tmp.afa'): os.remove('tmp.afa')

	cursor.close()
	db.close()
	return True


#########################################
if __name__ == '__main__':
	main()
