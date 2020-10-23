#!/usr/bin/python3 -u
# we are collecting orthologues for human only
# the possibility that "A is ortholog of B, and B is ortholog of C, but A is not ortholog of C"
# is something we don't want to know about right now


from time import time

from config import Config
from el_utils.ensembl import *
from el_utils.processes import *

def main():

	number_of_chunks = 1
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	ensembl_compara_name = get_compara_name(cursor)
	[all_species, ensembl_db_name] = get_species(cursor)

	# make sure we have the indices we need - does nothing if the index exists
	# create_index(cursor, db_name, table, index_name, columns, verbose=False)
	create_index(cursor, ensembl_compara_name, 'homology_member_human', 'gene_member_id', ['gene_member_id'], verbose=True)
	create_index(cursor, ensembl_compara_name, 'homology_member_human', 'homology_id', ['homology_id'], verbose=True)
	# strangely enough, this helps, even though gene_member_id is primary key
	# is it becasue I am using it in "for gene_meber in" construct,
	# and thies index dose not require uniquenss of the key (?) (see "show index from gene_member")
	create_index(cursor, ensembl_compara_name, 'gene_member', 'gene_member_id', ['gene_member_id'], verbose=True)

	# accomodate change between mysql versions 5.8 and later in allowed default datetime  values:
	# (if date is set to 0000-00-00 trying to create an index if date  0000-00-00 results in
	# Incorrect datetime value: '0000-00-00 00:00:00' for column 'created_date' at row 1)
	error_intolerant_search(cursor, "SET SQL_MODE='ALLOW_INVALID_DATES'")

	for db in ensembl_db_name.values():
		# stable id_id already exists, inherited form Ensembl itself,
		# defined on stable_id and version
		create_index(cursor, db, 'gene', 'stable_id', ['stable_id'], verbose=True)

	# this takes about 15 mins for ensembl 101
	table = "homology_human"
	if not check_table_exists(cursor, ensembl_compara_name, table):
		switch_to_db(cursor, ensembl_compara_name)
		qry = f"create table {table} engine=MYISAM  as select  h.homology_id, h.description "
		qry += "from homology h right join homology_member_human hh on hh.homology_id = h.homology_id"
		error_intolerant_search(cursor, qry)
	create_index(cursor, ensembl_compara_name, table, 'homology_id', ['homology_id'], verbose=True)

	cursor.close()
	db.close()

	return True


#########################################
if __name__ == '__main__':
	main()

'''
none of this worked:
create  table  scratch_gene_mem_subtable  engine=MYISAM  as select gene_member_id, 
stable_id from gene_member where taxon_id=9606  and biotype_group='coding' ;

Query OK, 23471 rows affected (0.66 sec)
Records: 23471  Duplicates: 0  Warnings: 0

'''

