#!/usr/bin/python3 -u
# we are collecting orthologues for human only
# the possibility that "A is ortholog of B, and B is ortholog of C, but A is not ortholog of C"
# is something we don't want to know about right now
from time import time

from config import Config
from el_utils.ensembl import *
from el_utils.processes import *


#########################################
# writing directly to db way too slow
def write_orthologues(ortho_file_handle, cursor, all_species, ensembl_db_name,  gene_id, orthos):
	for ortho in orthos:
		[ortho_stable, species, cognate_genome_db_id] = ortho
		if species not in all_species: continue
		ortho_gene_id = stable2gene (cursor, ortho_stable, ensembl_db_name[species])
		ortho_file_handle.write("{}\t{}\t{}\t{}\n".format(gene_id, ortho_gene_id, cognate_genome_db_id, 'ensembl'))
		ortho_file_handle.flush()
	return


#########################################
# ./el_utils/kernprof.py -l <calling script>.py
# python3 -m line_profiler <calling script>.py.lprof
#@profile
def collect_orthologues(gene_list, dummy):

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species (cursor)

	db_human  = connect_to_mysql(Config.mysql_conf_file)
	cursor_human = db_human.cursor()
	switch_to_db(cursor_human, ensembl_db_name['homo_sapiens'])
	search_db(cursor_human, "set autocommit=1")

	ensembl_compara_name = get_compara_name(cursor)

	db_compara     = connect_to_mysql(Config.mysql_conf_file)
	cursor_compara = db_compara.cursor()
	switch_to_db(cursor_compara, ensembl_compara_name)

	ortho_table = {'ortholog_one2one': 'orthologue', 'apparent_ortholog_one2one': 'orthologue',
					'possible_ortholog': 'unresolved_ortho', 'ortholog_one2many': 'unresolved_ortho',
					'ortholog_many2many': 'unresolved_ortho'}

	pid= os.getpid()
	filehandle = {'orthologue':open("raw_tables/orthologue.{}.tsv".format(pid), "w"),
					'unresolved_ortho':open("raw_tables/unresolved_ortho.{}.tsv".format(pid), "w")}

	ct = 0
	time0 = time()
	for gene_id in gene_list:
		# find stable
		stable_id = gene2stable(cursor_human, gene_id=gene_id)
		# memebr id refers to entries in compara db
		member_id = stable2member(cursor_compara, stable_id)

		orthos = get_orthologues(cursor_compara, ensembl_compara_name, member_id, verbose=False)

		# in compara table, get everything that homology has to say about
		# the possible orthologues
		# find all orthologous pairs suggested for this gene
		for ortho_type in ['ortholog_one2one','possible_ortholog', 'apparent_ortholog_one2one',
							'ortholog_one2many','ortholog_many2many']:
			if (not ortho_type in orthos) or (not orthos[ortho_type]): continue
			# the triple returned for each ortho type is [ortho_stable, species,  genome_db_id]
			table = ortho_table[ortho_type]

			#write_orthologues(filehandle[table], cursor, all_species, ensembl_db_name,  gene_id, orthos[ortho_type])
		ct += 1
		if not ct%100:
			print(os.getpid(), ct, "out of ", len(gene_list), "  last 100: %d mins"%((time()-time0)/60))
			time0 = time()

	filehandle['orthologue'].close()
	filehandle['unresolved_ortho'].close()


	cursor.close()
	db.close()
	cursor_human.close()
	db_human.close()
	cursor_compara.close()
	db_compara.close()


#########################################
def main():

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
		time_qry(cursor,qry, verbose=True)

		print(hard_landing_search(cursor,"select count(*) from %s"%tmp_table))


	[all_species, ensembl_db_name] = get_species(cursor)
	db_name = ensembl_db_name[species]
	switch_to_db (cursor, ensembl_db_name[species])
	gene_list = get_gene_ids(cursor, db_name, biotype="protein_coding", stable=True)
	for gene_id in gene_list[:100]:
		qry = "select homology_id from %s where stable_id=%s" % (tmp_table, gene_id)
		ret = error_intolerant_search(cursor, qry)
		print (gene_id, len(ret))


	cursor.close()
	db.close()

	#parallelize(no_threads, collect_orthologues, gene_list, [ensembl_db_name])

	return True


#########################################
if __name__ == '__main__':
	main()
