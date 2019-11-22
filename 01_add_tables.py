#!/usr/bin/python3

from config import Config

from el_utils.ensembl import get_species
from el_utils.ncbi    import *


# for Incorrect datetime value: '0000-00-00 00:00:00' for column 'created_date' at row 1
# when trying to do alter table, see https://stackoverflow.com/questions/35565128/mysql-incorrect-datetime-value-0000-00-00-000000/35565866


#########################################
def make_novel_exon_table (cursor, table):
	# this is a silly way to create tables - needs rewrite

	# if maps_to_human_exon_id is 0
	# the region was searched, but nothing was found
	# has NNN refers to the fact that the searched region  contains NNN stretch,
	#     indicating that the exon might not have been sequenced
	# 5p_ss and 3p_ss refer to canonical splice sites -r' and t'  refer to the intron

	qry  = "CREATE TABLE " + table + "  (exon_id INT(10) PRIMARY KEY AUTO_INCREMENT)"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	for column_name in ['gene_id', 'start_in_gene',
						'end_in_gene', 'maps_to_human_exon_id', 'exon_seq_id', 'template_exon_seq_id']:
		qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['template_species']:
		qry = "ALTER TABLE %s  ADD %s VARCHAR(120)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['strand', 'phase', 'end_phase',  'has_NNN', 'has_stop']:
		qry = "ALTER TABLE %s  ADD %s tinyint" %  (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False
	for column_name in ['has_3p_ss', 'has_5p_ss']:
		qry = "ALTER TABLE %s  ADD %s blob" %  (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False


#########################################
def make_gene2exon_table (cursor):

	table = 'gene2exon'

	qry  = "CREATE TABLE " + table + " (gene2exon_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	for column_name in ['gene_id', 'exon_id', 'start_in_gene', 'end_in_gene',
						'canon_transl_start', 'canon_transl_end', 'exon_seq_id']:
		qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['strand', 'phase',  'is_known',
						'is_coding', 'is_canonical', 'is_constitutive']:
		qry = "ALTER TABLE %s  ADD %s tinyint" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['covering_exon']:
		qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['covering_is_known']:
		qry = "ALTER TABLE %s  ADD %s tinyint" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['analysis_id']:
		qry = "ALTER TABLE %s  ADD %s INT" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

#########################################
def make_exon_seq_table (cursor):


	table = 'exon_seq'

	qry  = "CREATE TABLE " + table + "  (exon_seq_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	for column_name in ['exon_id']:
		qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False


	for column_name in ['is_known', 'is_sw']:
		qry = "ALTER TABLE %s  ADD %s tinyint" %  (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['dna_seq', 'left_flank', 'right_flank', 'protein_seq']:
		qry = "ALTER TABLE %s  ADD %s blob" %  (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	qry  = "alter TABLE " + table + "  add column pepseq_transl_start int(10)"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	qry  = "alter TABLE " + table + "  add column pepseq_transl_end int(10)"
	rows = search_db (cursor, qry)
	if (rows):
		return False


#########################################
def make_coding_region_table(cursor):


	table = 'coding_region'

	qry  = "CREATE TABLE " + table + " (gene_id INT PRIMARY KEY)"
	rows = search_db (cursor, qry)
	if (rows):
		return False


	# make the columns
	columns = ['start', 'end']
	for column in columns:
		qry = "ALTER TABLE %s ADD %s  INT(10)" % (table, column)
		rows = search_db (cursor, qry)
		if (rows):
			return False


#########################################
def make_orthologue_table (cursor, table):


	# if congate_gene_id is 0, and source is 'rbh'
	# means that the reciprocal-best-hit was attempted but nothing was found

	qry  = "CREATE TABLE " + table + "  (orth_pair_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	for column_name in ['gene_id', 'cognate_gene_id']:
		qry = "ALTER TABLE %s  ADD %s INT(10)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['cognate_genome_db_id']:
		qry = "ALTER TABLE %s  ADD %s INT" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False


	for column_name in ['source']:
		qry = "ALTER TABLE %s  ADD %s VARCHAR(20)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False


#########################################
def modify_exon_map_table (cursor):

	table = 'exon_map'

	for column_name in ['warning']:
		qry = "ALTER TABLE %s add %s blob" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False


#########################################
def make_exon_map_table (cursor):

	table = 'exon_map'

	qry  = "CREATE TABLE " + table + "  (exon_map_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	for column_name in ['exon_id', 'cognate_exon_id']:
		qry = "ALTER TABLE %s add %s INT(10)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['exon_known', 'cognate_exon_known']:
		qry = "ALTER TABLE %s add %s tinyint" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['cognate_genome_db_id']:
		qry = "ALTER TABLE %s add %s INT" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['cigar_line']:
		qry = "ALTER TABLE %s add %s blob" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['similarity']:
		qry = "ALTER TABLE %s add %s float" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['source']:
		qry = "ALTER TABLE %s add %s VARCHAR(20)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['msa_bitstring']:
		qry = "ALTER TABLE %s add %s varbinary(1000)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False


	for column_name in ['warning']:
		qry = "ALTER TABLE %s add %s blob" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False


#########################################
def make_para_exon_map_table (cursor):

	table = 'para_exon_map'

	qry  = "CREATE TABLE " + table + "  (exon_map_id INT(10)  PRIMARY KEY AUTO_INCREMENT)"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	for column_name in ['exon_id', 'cognate_exon_id']:
		qry = "ALTER TABLE %s add %s INT(10)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['exon_known', 'cognate_exon_known']:
		qry = "ALTER TABLE %s add %s tinyint" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['cigar_line']:
		qry = "ALTER TABLE %s add %s blob" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['similarity']:
		qry = "ALTER TABLE %s add %s float" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['source']:
		qry = "ALTER TABLE %s add %s VARCHAR(20)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

	for column_name in ['msa_bitstring']:
		qry = "ALTER TABLE %s add %s varbinary(1000)" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False


#########################################
def make_para_groups_table(cursor, db_name):
	switch_to_db (cursor, db_name)
	table = 'paralogue_groups'
	if check_table_exists(cursor, db_name, table):
		qry = "drop table " + table
		search_db(cursor, qry, verbose=True)

	qry = ""
	qry += "CREATE TABLE  %s (" % table
	qry += "     group_id INT NOT NULL AUTO_INCREMENT, "
	qry += "  	 stable_ids text NOT NULL, "
	qry += "	 PRIMARY KEY (group_id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)


#########################################
def make_table (cursor, db_name, table):

	qry = "use %s" % db_name
	rows = search_db (cursor, qry, verbose=False)
	if (rows):
		return False

	if   table == 'gene2exon':
		make_gene2exon_table (cursor)
	elif table == 'exon_seq':
		make_exon_seq_table (cursor)
	elif table == 'sw_exon':
		make_novel_exon_table (cursor, table)
	elif table == 'usearch_exon':
		make_novel_exon_table (cursor, table)
	elif table == 'coding_region':
		make_coding_region_table (cursor)
	elif table in ['orthologue', 'unresolved_ortho']:
		make_orthologue_table (cursor, table)
	elif table == 'paralogue_groups':
		make_para_groups_table(cursor, db_name)
	elif table == 'exon_map':
		make_exon_map_table (cursor)
	elif table == 'para_exon_map':
		make_para_exon_map_table (cursor)

	else:
		print("I don't know how to make table '%s'" % table)



#########################################
def add_filename_column (cursor, db_name):

	qry = "alter table seq_region add file_name blob"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	return True


#########################################
def modify_filename_column (cursor, db_name):

	qry = "alter table seq_region modify column  file_name blob"
	rows = search_db (cursor, qry)
	if (rows):
		return False

	return True


#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	# the 000 date problem
	qry = " SET @@sql_mode :='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,ERROR_FOR_DIVISION_BY_ZERO,NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION'"
	error_intolerant_search(cursor,qry)

	# add exon tables to all species
	for species in all_species:

		db_name = ensembl_db_name[species]
		switch_to_db (cursor, ensembl_db_name[species])

		for table in ['gene2exon', 'exon_seq', 'sw_exon', 'usearch_exon', 'coding_region']:
			if check_table_exists(cursor, db_name, table):
				print(table, " found in ", db_name)
			else:
				print(table, " not found in ", db_name)
				make_table (cursor, db_name, table)

		print("optimizing gene2exon")
		qry = "optimize table gene2exon"
		print(search_db(cursor, qry))
		create_index (cursor, db_name, 'eg_index',    'gene2exon', ['exon_id', 'gene_id'])
		create_index (cursor, db_name, 'gene_id_idx', 'gene2exon', ['gene_id'])
		create_index (cursor, db_name, 'ek_index',    'exon_seq',  ['exon_id', 'is_known'])
		create_index (cursor, db_name, 'seq_index',   'exon_seq',  ['exon_seq_id'])
		print("optimizing exon_seq")
		qry = "optimize table exon_seq"
		print(search_db(cursor, qry))

	# add file_name column to seq_region table (seq_region table  already exists in ensembl schema)
	for species in all_species:
		print(species)
		db_name = ensembl_db_name[species]

		if column_exists(cursor, db_name, "seq_region", "file_name"):
			print("file_name found in seq_region, ", db_name)

		else:  # modify_filename_column (cursor, db_name)
			print("file_name  not found in seq_region, ", db_name, "(making the column)")
			add_filename_column (cursor, db_name)

	# add orthologue table to human - we are human-centered here
	# ditto for map (which exons from other species map onto human exons)
	print("adding orthologue to human")
	species = 'homo_sapiens'
	db_name = ensembl_db_name[species]
	for table in ['orthologue', 'unresolved_ortho', 'paralogue', 'exon_map']:
		if check_table_exists(cursor, db_name, table):
			print(table, " found in ", db_name)
		else:
			print(table, " not found in ", db_name)
			make_table (cursor, db_name, table)
		if table == 'exon_map':
			create_index(cursor, db_name,'gene_index', table, ['exon_id'])
			create_index(cursor, db_name,'exon_index', table, ['exon_id', 'exon_known'])
			create_index(cursor, db_name,'cognate_exon_index', table, ['cognate_exon_id', 'cognate_exon_known', 'cognate_genome_db_id'])
		else:
			create_index (cursor, db_name,'gene_index', table, ['gene_id'])

	# for other species, add paralogue map pnly
	for species in all_species:
		db_name = ensembl_db_name[species]
		# add_column checks whether column exists
		add_column(cursor, db_name, 'gene', 'paralogue_group_ids', "text")
		for table in ['paralogue_groups', 'para_exon_map']:
			if check_table_exists(cursor, db_name, table):
				print(table, " found in ", db_name)
			else:
				print(table, " not found in ", db_name)
				make_table(cursor, db_name, table)
			if table == 'para_exon_map':
				create_index(cursor, db_name,'gene_index', table, ['exon_id'])


	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
