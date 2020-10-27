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
		# name = "_".join([sp[:3] for sp in species.split("_")]) # +[gene_name])
		name = species
		cmd = "echo '>{}' >> {}".format(name, out_fasta)
		subprocess.call(["bash","-c", cmd])
		cmd = "cat {} >> {}".format(tmpfile, out_fasta)
		subprocess.call(["bash","-c", cmd])
		return True

	return False


def reorder_seqs(in_afa, species_sorted, out_afa, trivial=None, name_prefix=None):
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
		for name in species_sorted:
			outname = trivial[name] if trivial else name
			if name_prefix: outname = f"{name_prefix}_{outname}"
			outf.write(f">{outname}\n")
			for i in range(0,len(seq[name]), 100):
				outf.write(f"{seq[name][i:i+100]}\n")


#########################################
def main():

	if len(sys.argv)<2:
		print("usage: %s <gene symbol> [trivial] [prepend]" % sys.argv[0])
		print("trivial = use trivial species name; prepend = prepend gene name")
		exit()
	gene_name = sys.argv[1]
	trivial = "trivial" in sys.argv
	prepend = "prepend" in sys.argv # prepends geen synbol to gene name
	ref_species = 'homo_sapiens'  # the orthologue table is filled only here, for the moment

	out_fasta = f"{gene_name}.orthos.fasta"
	out_afa = f"{gene_name}.orthos.afa"
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
	# species_sort(cursor, all_species, 'homo_sapiens')
	# exit()

	switch_to_db(cursor, ensembl_db_name[ref_species])
	qry = "select gene_id from gene where stable_id='%s'" % ensembl_stable_gene_id
	gene_id = hard_landing_search(cursor, qry)[0][0]

	ref_stable_transl_id = gene2stable_canon_transl_id(cursor, gene_id,  ensembl_db_name[ref_species])
	write_to_fasta(home, ref_species, ref_stable_transl_id, tmpfile, logfile, out_fasta)

	print(gene_name, ensembl_stable_gene_id, gene_id, ref_stable_transl_id)
	species_in_the_almt = [ref_species]
	qry = "select  cognate_gene_id, cognate_genome_db_id from orthologues where gene_id=%d" % gene_id
	for line in error_intolerant_search(cursor, qry):
		[cognate_gene_id, cognate_genome_db_id] = line
		qry = f"select db_name from exolocator_meta.db_names where genome_db_id={cognate_genome_db_id}"
		db_name = hard_landing_search(cursor, qry)[0][0]
		stable_transl_id = gene2stable_canon_transl_id(cursor, cognate_gene_id, db_name)
		species = db_name.split("core")[0].rstrip("_")
		if species not in all_species: continue
		print(db_name,  species, cognate_gene_id, stable_transl_id)
		ok = write_to_fasta(home, species, stable_transl_id, tmpfile, logfile, out_fasta)
		if ok: species_in_the_almt.append(species)
	if os.path.exists(tmpfile): os.remove(tmpfile)

	cmd  = f"{Config.muscle} -in {out_fasta} -out tmp.afa"
	subprocess.call(["bash","-c", cmd])

	species_sorted = species_sort(cursor, species_in_the_almt, ref_species)
	trivial_names = get_trivial(cursor, species_sorted) if trivial else None
	name_prefix = gene_name if prepend else None
	reorder_seqs('tmp.afa', species_sorted, out_afa, trivial_names, name_prefix)
	if os.path.exists('tmp.afa'): os.remove('tmp.afa')

	cursor.close()
	db.close()
	return True


#########################################
if __name__ == '__main__':
	main()
