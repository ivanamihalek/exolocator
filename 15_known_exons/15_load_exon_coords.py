#!/usr/bin/python3


from config import Config
from el_utils.ensembl import *

#########################################
# problem:
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
	print("careful, this script deletes contents of gene2exon table")
	exit()
	in_dir = "raw_tables"
	table  = "gene2exon"
	if not os.path.exists(in_dir):
		print(in_dir,"not found")
		exit()

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	for species in all_species:
		print(species)
		infile = f"{in_dir}/{species}/{table}.tsv"
		for dep in [f"{in_dir}/{species}", infile]:
			if not os.path.exists(dep):
				print(dep,"not found")
				exit()
		error_intolerant_search(cursor, f"delete from {ensembl_db_name[species]}.{table}")
	cursor.close()
	db.close()

	for species in all_species:
		infile = f"{in_dir}/{species}/{table}.tsv"
		cmd = f"mysqlimport --login-path=tcga --fields_escaped_by=\\\\ {ensembl_db_name[species]} -L {infile}"
		print(cmd)
		subprocess.call(["bash","-c", cmd])


#####################################################
if __name__=="__main__":
	main()
