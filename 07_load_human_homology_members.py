#! /usr/bin/python3

from el_utils.mysql import *
from el_utils.ensembl import *
from config import Config
import os

#########################################/home/ivana/.tcga_conf
# problem:
# if try to load with
# load data infile from within mysql
# "The MySQL server is running with the --secure-file-priv option so it cannot execute this statement"
# if you do not wish to  to reconfigure and restart the server,
# find the privileged directory with
# mysql> SHOW VARIABLES LIKE "secure_file_priv";
# | Variable_name    | Value                 |
# +------------------+-----------------------+
# | secure_file_priv | /var/lib/mysql-files/ |
# move files that need to be loaded tothat dir
# remove later

# cleaner solution:
# for this to work conf file must have the local_infile=1 line
# and in mysql one has to run
# mysql>  SET GLOBAL local_infile = 'ON';
# check
# mysql> SHOW GLOBAL VARIABLES LIKE 'local_infile';
# not sure how thid variable gets flipped back to OFF - on server restart?
# table name and tsv name must match
# Then load with mysqlimport (-L stands for local) - table and the filename without the extension have to match
# "mysqlimport --login-path=tcga --fields_escaped_by=\\\\ $db -L *.txt"

def make_human_hom_member_table(cursor, db_name, table):
	check_and_drop_table(cursor, db_name, table)
	switch_to_db(cursor, db_name)
	qry = ""
	qry += f"CREATE TABLE  {table} ("
	qry += "     homology_id bigint unsigned not null, "
	qry += "  	 gene_member_id int unsigned not null, "
	qry += "	 PRIMARY KEY (homology_id, gene_member_id) "
	qry += ") ENGINE=MyISAM"
	error_intolerant_search(cursor, qry)
	return



def main():
	print("careful, this script deletes table contents")
	exit()
	# in ensembl 101 this take about 1 min to read in
	# check out, though, the previous scrip for the time to produce this table
	home = os.getcwd()
	table_name = "homology_member_human"
	infile = f"{home}/raw_tables/{table_name}.tsv"
	if not os.path.exists(infile):
		print(f"{infile} not found")
		exit()

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")

	ensembl_compara_name = get_compara_name(cursor)
	make_human_hom_member_table(cursor, ensembl_compara_name, table_name)
	cursor.close()
	db.close()

	cmd = f"mysqlimport --login-path=tcga --fields_escaped_by=\\\\ {ensembl_compara_name} -L {infile}"
	subprocess.call(["bash","-c", cmd])



#########################################
if __name__ == '__main__':
	main()
