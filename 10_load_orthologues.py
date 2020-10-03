#!/usr/bin/python3


from config import Config
from el_utils.ensembl import *

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

# cleaner solution:
# for this to work conf file must have the local_infile=1 line
# and in mysql one has to run
# mysql>  SET GLOBAL local_infile = 'ON';
# check
# mysql> SHOW GLOBAL VARIABLES LIKE 'local_infile';
# not sure how thid variable gets flipped back to OFF - on server restart, apparently
# table name and tsv name must match <<<< !!! NOTE THIS
# Then load with mysqlimport (-L stands for local) - table and the filename without the extension have to match
# "mysqlimport --login-path=tcga --fields_escaped_by=\\\\ $db -L *.txt"


def make_orthologues_table(cursor, db_name):
	check_and_drop_table(cursor, db_name, "orthologues")
	switch_to_db(cursor, db_name)
	qry = ""
	qry += f"CREATE TABLE  orthologues ("
	qry += "     gene_id int unsigned not null, "
	qry += "  	 cognate_gene_id int unsigned not null, "
	qry += "  	 cognate_genome_db_id int unsigned not null"
	#qry += "	 PRIMARY KEY (gene_id) " # I have duplicate gene_ids here!
	qry += ") ENGINE=MyISAM"
	error_intolerant_search(cursor, qry)
	return


####################################################
def main():
	print("careful, this script deletes orthologues table")
	exit()

	# we'll do this for human only, for now
	species = 'homo_sapiens'
	in_dir = "raw_tables"
	for dep in [in_dir, f"{in_dir}/orthologues.tsv"]:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	[all_species, ensembl_db_name] = get_species(cursor)

	make_orthologues_table(cursor, ensembl_db_name[species])
	cursor.close()
	db.close()

	cmd = f"mysqlimport --login-path=tcga --fields_escaped_by=\\\\ {ensembl_db_name[species]} -L {in_dir}/orthologues.tsv"
	subprocess.call(["bash","-c", cmd])


#####################################################
if __name__=="__main__":
	main()
