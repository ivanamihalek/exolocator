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



####################################################
def main():
	# print("careful, this script deletes contents of exon_seq table")
	# exit()
	in_dir = "raw_tables"
	table  = "exon_seq"
	if not os.path.exists(in_dir):
		print(in_dir,"not found")
		exit()

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	cursor.close()
	db.close()
	for species in all_species:
		infile = f"{in_dir}/{species}/{table}.tsv"
		for dep in [f"{in_dir}/{species}", infile]:
			if not os.path.exists(dep):
				print(dep,"not found")
				exit()

	for species in all_species[1:]:
		print(species)
		infile = f"{in_dir}/{species}/{table}.tsv"
		outfile =  f"{in_dir}/{species}/{table}.2.tsv"
		outf = open(outfile, "w")
		with open(infile) as inf:
			for line in inf:
				fields = line.split("\t")
				newfields = fields[:5]+fields[6:]
				outf.write("\t".join(newfields))
		outf.close()
		os.rename(outfile, infile)

#####################################################
if __name__=="__main__":
	main()
