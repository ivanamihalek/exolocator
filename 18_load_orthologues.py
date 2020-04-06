#! /usr/bin/python3

from el_utils.mysql import *
from el_utils.ensembl import *
from config import Config
import os

#########################################
# problem:
# "The MySQL server is running with the --secure-file-priv option so it cannot execute this statement"
# if you do not wish to  to reconfigure and restart the server,
# find the privileged directory with
# mysql> SHOW VARIABLES LIKE "secure_file_priv";
# | Variable_name    | Value                 |
# +------------------+-----------------------+
# | secure_file_priv | /var/lib/mysql-files/ |
# move files that need to be loaded tothat dir
# remove later


def main():
	print("careful, this script deletes table contents")
	exit()
	no_threads = 1
	species = 'homo_sapiens'
	indir = "/var/lib/mysql-files/"
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")

	[spec_names, ensembl_db_name] = get_species(cursor)
	switch_to_db(cursor, ensembl_db_name[species])

	for table_name in ["unresolved_ortho", "orthologue"]:
		print("deleting from ", table_name)
		error_intolerant_search(cursor, f"delete from {table_name}")
		print("resetting auto increment ", table_name)
		error_intolerant_search(cursor, f"alter table {table_name} auto_increment=1")
		for fnm in os.listdir(indir):
			if not table_name in fnm: continue
			print("loading", fnm)
			qry  = f"LOAD DATA INFILE '{indir}/{fnm}' INTO TABLE  "
			qry += f"{table_name} (gene_id, cognate_gene_id, cognate_genome_db_id, source)"
			print(qry)
			error_intolerant_search(cursor, qry)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
