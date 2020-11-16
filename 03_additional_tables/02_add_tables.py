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
def make_gene2exon_table(cursor):

	table = 'gene2exon'
	# no auto_increment for the gene2exon table - we will be reading it from  the tsv
	qry = ''' CREATE TABLE `gene2exon` (
			  `gene2exon_id` int NOT NULL,
			  `gene_id` int NOT NULL,
			  `exon_id` int NOT NULL,
			  `start_in_gene` int DEFAULT NULL,
			  `end_in_gene` int DEFAULT NULL,
			  `canon_transl_start` int DEFAULT NULL,
			  `canon_transl_end` int DEFAULT NULL,
			  `exon_seq_id` int DEFAULT NULL,
			  `strand` tinyint DEFAULT NULL,
			  `phase` tinyint DEFAULT NULL,
			  `provenance` tinyint DEFAULT NULL,
			  `is_coding` tinyint DEFAULT NULL,
			  `is_canonical` tinyint DEFAULT NULL,
			  `is_constitutive` tinyint DEFAULT NULL,
			  `covering_exon` int DEFAULT NULL,
			  `covering_provenance` tinyint DEFAULT NULL,
			  `analysis_id` int DEFAULT NULL,
			  PRIMARY KEY (`gene2exon_id`),
			  KEY `exon_id_idx` (`exon_id`)
		)  ENGINE=MyISAM '''
	error_intolerant_search(cursor, qry)


#########################################
def make_exon_seq_table(cursor):
	# error_intolerant_search(cursor, "drop table if exists exon_seq")
	qry = ""
	qry += "  CREATE TABLE  exon_seq ("
	qry += "     exon_seq_id INT PRIMARY KEY AUTO_INCREMENT, "
	qry += "     exon_id INT, "
	qry += "     phase tinyint, "
	qry += "     by_exolocator tinyint, "
	qry += "	 dna_seq text, "
	qry += "	 left_flank text, "
	qry += "	 right_flank text, "
	qry += "	 protein_seq text "
	qry += ") ENGINE=MyISAM"
	error_intolerant_search(cursor, qry)


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
def make_problems_table(cursor):

	qry = ""
	qry += "  CREATE TABLE  problems ("
	qry += "     gene_id INT PRIMARY KEY, "
	qry += "	 description text"
	qry += ") ENGINE=MyISAM"
	error_intolerant_search(cursor, qry)

#########################################
def modify_exon_map_table (cursor):

	table = 'exon_map'

	for column_name in ['warning']:
		qry = "ALTER TABLE %s add %s blob" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False

#########################################
def make_table (cursor, db_name, table):

	qry = "use %s" % db_name
	rows = search_db (cursor, qry, verbose=False)
	if rows:
		return False

	if table == 'gene2exon':
		make_gene2exon_table (cursor)
	elif table == 'exon_seq':
		make_exon_seq_table (cursor)
	elif table == 'sw_exon':
		make_novel_exon_table (cursor, table)
	elif table == 'usearch_exon':
		make_novel_exon_table (cursor, table)
	elif table == 'coding_region':
		make_coding_region_table (cursor)
	elif table == 'problems':
		make_problems_table (cursor)
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

	# add exon tables to all species
	for species in all_species:
		print(species)
		db_name = ensembl_db_name[species]
		switch_to_db (cursor, db_name)
		#make_exon_seq_table(cursor)

		#for table in ['gene2exon', 'exon_seq', 'sw_exon', 'usearch_exon', 'coding_region', 'problems']:
		for table in ['gene2exon']:
			check_and_drop_table(cursor, db_name, table)
			make_table (cursor, db_name, table)
			# if check_table_exists(cursor, db_name, table):
			# 	print(table, " found in ", db_name)
			# else:
			# 	print(table, " not found in ", db_name)
			# 	make_table (cursor, db_name, table)

		print("optimizing gene2exon")
		qry = "optimize table gene2exon"
		print(search_db(cursor, qry))
		# (cursor, db_name, table, index_name, columns, verbose=False)
		create_index (cursor, db_name,  'gene2exon',  'eg_index',   ['exon_id', 'gene_id'], verbose=True)
		create_index (cursor, db_name, 'gene2exon', 'gene_id_idx', ['gene_id'], verbose=True)
		# create_index (cursor, db_name, 'ek_index',    'exon_seq',  ['exon_id', 'is_known'])
		# create_index (cursor, db_name, 'seq_index',   'exon_seq',  ['exon_seq_id'])
		# print("optimizing exon_seq")
		# qry = "optimize table exon_seq"
		# print(search_db(cursor, qry))


	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
