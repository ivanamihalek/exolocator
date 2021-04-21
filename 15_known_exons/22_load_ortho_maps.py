#!/usr/bin/python3 -u

from el_utils.map import *
#########################################

#########################################
def main():

	indir = "raw_tables"
	if not os.path.exists(indir):
		print(f"{indir} not found")
		exit()

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	search_db(cursor, "set global local_infile = 'ON'")

	[all_species, ensembl_db_name] = get_species (cursor)
	members = {}
	qry = "select representative_species, members from exolocator_meta.taxonomy_groups "
	for rep_species, members_str in hard_landing_search(cursor, qry):
		members[rep_species] = members_str.split(",")

	table = "exon_map"
	for rep_species in members.keys():
		if len(members[rep_species])==1: continue
		infile = f"{indir}/{rep_species}/{table}.tsv"
		dbname = ensembl_db_name[rep_species]
		print(f"{rep_species} number of members: {len(members[rep_species])}  writing to {dbname}.{table}")
		if not os.path.exists(infile):
			print(f"input file {infile} not found")
			exit()
		if not os.path.getsize(infile) > 0:
			print(f"input file {infile} appears to be empty")
			exit()

		cmd = f"mysqlimport --login-path=tcga --fields_escaped_by=\\\\ {dbname} -L {infile}"
		print(cmd)
		subprocess.call(["bash","-c", cmd])

	cursor.close()
	db.close()

	return True


#########################################
if __name__ == '__main__':
	main()

