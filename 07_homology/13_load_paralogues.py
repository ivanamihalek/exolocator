#!/usr/bin/python3


from config import Config
from el_utils.ensembl import *

#########################################
# problem:
# "The MySQL server is running with the --secure-file-priv option so it cannot execute this statement"
# if you do not wish to  to reconfigure and restart the server,

# for this to work conf file must have the local_infile=1 line
# and in mysql one has to run
# mysql>  SET GLOBAL local_infile = 'ON';
# check
# mysql> SHOW GLOBAL VARIABLES LIKE 'local_infile';
# not sure how thid variable gets flipped back to OFF - on server restart, apparently
# table name and tsv name must match <<<< !!! NOTE THIS
# Then load with mysqlimport (-L stands for local) - table and the filename without the extension have to match
# "mysqlimport --login-path=tcga --fields_escaped_by=\\\\ $db -L *.txt"


#########################################
def make_para_groups_table(cursor, db_name):
	switch_to_db (cursor, db_name)
	table = 'paralogue_groups'
	if check_table_exists(cursor, db_name, table):
		qry = "drop table " + table
		search_db(cursor, qry, verbose=True)

	qry = ""
	qry += "CREATE TABLE  %s (" % table
	qry += "     bait_stable_id varchar (128), "
	qry += "  	 stable_ids text NOT NULL, "
	qry += "	 PRIMARY KEY (bait_stable_id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)



####################################################
def main():
	# print("careful, this script deletes a table")
	# exit()

	# we'll do this for human only, for now
	species = 'homo_sapiens'
	in_dir  = "raw_tables"
	table_name = "paralogue_groups"
	fnm = f"{table_name}.tsv"
	for dep in [in_dir, f"{in_dir}/{fnm}"]:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	[all_species, ensembl_db_name] = get_species(cursor)

	switch_to_db (cursor, ensembl_db_name[species])

	make_para_groups_table(cursor, ensembl_db_name['homo_sapiens'])
	cursor.close()
	db.close()

	print("loading", fnm)
	cmd = f"mysqlimport --login-path=tcga --fields_escaped_by=\\\\ {ensembl_db_name[species]} -L {in_dir}/{fnm}"
	subprocess.call(["bash","-c", cmd])


#####################################################
if __name__=="__main__":
	main()
