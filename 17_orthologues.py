#!/usr/bin/python3 -u
# we are collecting orthologues for human only
# the possibility that "A is ortholog of B, and B is ortholog of C, but A is not ortholog of C"
# is something we don't want to know about right now
from time import time

from config import Config
from el_utils.ensembl import *
from el_utils.processes import *


############
def get_orthologues(cursor, compara_db, tmp_table, gene_stable_id, verbose=False):
	orthos = {}
	qry  = "select gm.stable_id, gm.genome_db_id, hml.description "
	qry += "from %s.gene_member as gm, " % compara_db
	qry += "%s.homology_member as hm, " % compara_db
	qry += "%s.homology as hml, " % compara_db
	qry += "%s.%s as s " % (compara_db, tmp_table)
	qry += "where gm.gene_member_id = hm.gene_member_id and hm.homology_id = s.homology_id "
	qry += "and hm.homology_id = hml.homology_id "
	qry += "and s.stable_id = '%s' " % gene_stable_id
	qry += "and gm.stable_id != '%s' " % gene_stable_id
	ret = error_intolerant_search(cursor, qry)
	if not ret: return []
	for line in ret:
		[orthologue_stable_id, genome_db_id, homology_description] = line
		qry = "select db_name from exolocator_meta.db_names where genome_db_id=%d" % genome_db_id
		ret = error_intolerant_search(cursor,qry)
		if not ret: continue
		db_name = ret[0][0]
		qry = f"select gene_id from {db_name}.gene where stable_id='{orthologue_stable_id}'"
		ortho_gene_id = hard_landing_search(cursor, qry)[0][0]
		# print(orthologue_stable_id, db_name, homology_description)
		if not homology_description in orthos:  orthos[homology_description] = []
		orthos[homology_description].append([ortho_gene_id, genome_db_id])

	return orthos


#########################################
# writing directly to db way too slow
def write_orthologues(ortho_file_handle, gene_id, orthos):
	for ortho in orthos:
		[ortho_gene_id, ortho_genome_db_id] = ortho
		ortho_file_handle.write("{}\t{}\t{}\t{}\n".format(gene_id, ortho_gene_id, ortho_genome_db_id, 'ensembl'))
		ortho_file_handle.flush()
	return


#########################################
# ./el_utils/kernprof.py -l <calling script>.py
# python3 -m line_profiler <calling script>.py.lprof
#@profile
def collect_orthologues(stable_gene_id_list, dummy):

	[ensembl_compara_name, ensembl_db_name, tmp_table] = dummy
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor, ensembl_db_name['homo_sapiens'])
	search_db(cursor, "set autocommit=1")

	ortho_table = {'ortholog_one2one': 'orthologue', 'apparent_ortholog_one2one': 'orthologue',
					'possible_ortholog': 'unresolved_ortho', 'ortholog_one2many': 'unresolved_ortho',
					'ortholog_many2many': 'unresolved_ortho'}

	pid= os.getpid()
	filehandle = {'orthologue':open("raw_tables/orthologue.{}.tsv".format(pid), "w"),
					'unresolved_ortho':open("raw_tables/unresolved_ortho.{}.tsv".format(pid), "w")}

	ct = 0
	time0 = time()
	for stable_id in stable_gene_id_list:
		gene_id = stable2gene(cursor, stable_id,  ensembl_db_name['homo_sapiens'])
		# in compara table, get everything that homology has to say about
		# the possible orthologues
		# find all orthologous pairs suggested for this gene
		orthos = get_orthologues(cursor, ensembl_compara_name, tmp_table, stable_id, verbose=False)

		for ortho_type in ['ortholog_one2one','possible_ortholog', 'apparent_ortholog_one2one',
							'ortholog_one2many','ortholog_many2many']:
			if (not ortho_type in orthos) or (not orthos[ortho_type]): continue
			# the triple returned for each ortho type is [ortho_stable, species,  genome_db_id]
			table = ortho_table[ortho_type]
			write_orthologues(filehandle[table], gene_id, orthos[ortho_type])
		ct += 1
		if not ct%100:
			print(os.getpid(), ct, "out of ", len(stable_gene_id_list), "  last 100: %d mins" % ((time() - time0) / 60))
			time0 = time()

	filehandle['orthologue'].close()
	filehandle['unresolved_ortho'].close()

	cursor.close()
	db.close()


#########################################
def main():
	no_threads = 1
	species = 'homo_sapiens'

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	ensembl_compara_name = get_compara_name(cursor)
	switch_to_db(cursor, ensembl_compara_name)

	tmp_table = "scratch_homology_table"

	if not check_table_exists(cursor, ensembl_compara_name, tmp_table):
		# we don't want it temp, because we want the subprocesses to use it
		qry  = "create  table  %s  engine=MYISAM  as " % tmp_table
		qry += "select homology_member.homology_id, homology_member.gene_member_id, subtable.stable_id "
		qry += "from homology_member right join "
		qry += "(select gene_member_id, stable_id from gene_member where stable_id like 'ENSG0%') as subtable "
		qry += "on homology_member.gene_member_id=subtable.gene_member_id"
		time_qry(cursor,qry, verbose=False)

		print(hard_landing_search(cursor,"select count(*) from %s"%tmp_table))

	# this check fro existence before attmepting create
	create_index(cursor, ensembl_compara_name,  tmp_table, "stable_id_idx", ["stable_id"], verbose=False)

	[all_species, ensembl_db_name] = get_species(cursor)
	db_name = ensembl_db_name[species]
	switch_to_db (cursor, ensembl_db_name[species])
	stable_gene_id_list = get_gene_ids(cursor, db_name, biotype="protein_coding", stable=True)
	# for gene_id in stable_gene_id_list[:100]:
	# 	qry = "select homology_id from %s.%s where stable_id='%s'" % (ensembl_compara_name, tmp_table, gene_id)
	# 	ret = error_intolerant_search(cursor, qry)
	# 	print (gene_id, len(ret))

	cursor.close()
	db.close()
	arguments =[ensembl_compara_name, ensembl_db_name, tmp_table]
	parallelize(no_threads, collect_orthologues, stable_gene_id_list, arguments)

	return True


#########################################
if __name__ == '__main__':
	main()
