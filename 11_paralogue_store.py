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

####################################################
def main():
	print("careful, this script deletes table contents")
	exit()

	# we'll do this for human only, for now
	species = 'homo_sapiens'
	in_dir = "/var/lib/mysql-files/"
	for dep in [in_dir]:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	[all_species, ensembl_db_name] = get_species(cursor)

	switch_to_db (cursor, ensembl_db_name[species])

	table_name = "paralogue_groups"

	print("deleting from ", table_name)
	error_intolerant_search(cursor, f"delete from {table_name}")

	for fnm in os.listdir(in_dir):
		if "paralogues" not in fnm: continue
		print("loading", fnm)
		qry  = f"LOAD DATA INFILE '{in_dir}/{fnm}' INTO TABLE  "
		qry += f"{table_name} (bait_stable_id, stable_ids)"
		print(qry)
		error_intolerant_search(cursor, qry)

	cursor.close()
	db.close()

#####################################################
if __name__=="__main__":
	main()
