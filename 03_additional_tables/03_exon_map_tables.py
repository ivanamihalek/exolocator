#!/usr/bin/python3

from config import Config

from el_utils.ensembl import get_species
from el_utils.ncbi    import *


# for Incorrect datetime value: '0000-00-00 00:00:00' for column 'created_date' at row 1
# when trying to do alter table, see https://stackoverflow.com/questions/35565128/mysql-incorrect-datetime-value-0000-00-00-000000/35565866


#########################################
########################################
from el_utils.tree import species_tree


def modify_exon_map_table (cursor):

	table = 'exon_map'

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
def make_exon_map_table (cursor, db_name):

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
		if rows:
			return False


	for column_name in ['warning']:
		qry = "ALTER TABLE %s add %s blob" % (table, column_name)
		rows = search_db (cursor, qry)
		if (rows):
			return False



#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)

	# add orthologue table to human - we are human-centered here
	# ditto for map (which exons from other species map onto human exons)
	table = "exon_map"
	qry = "select representative_species from exolocator_meta.taxonomy_groups"
	for representative_species in [line[0] for line in hard_landing_search(cursor, qry)]:
		print(f"adding exon_map to {representative_species}")
		db_name = ensembl_db_name[representative_species]

		if check_table_exists(cursor, db_name, table):
			print(table, " found in ", db_name)
		else:
			print(table, " not found in ", db_name)
			make_exon_map_table(cursor, db_name)
		#            cursor, db_name, table, index_name, columns,
		create_index(cursor, db_name, table, 'gene_index', ['exon_id'])
		create_index(cursor, db_name, table, 'exon_index', ['exon_id', 'exon_known'])
		create_index(cursor, db_name, table, 'cognate_exon_index', ['cognate_exon_id', 'cognate_exon_known', 'cognate_genome_db_id'])


	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
