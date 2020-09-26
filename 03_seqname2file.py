#!/usr/bin/python3

import os

from config import Config

from el_utils.ensembl       import  *



def collect_fasta_file_content(dna_path):
	fasta_files = []
	for r,d,files in os.walk(dna_path):
		for file in files:
			if not file[-3:] == ".fa":
				continue
			fasta_files.append(file)

	name2file = {}
	for file in fasta_files:
		print(dna_path, file)
		# fgrep is enough if not looking for regexp
		cmd = f"LC_ALL=C fgrep '>' {dna_path}/{file}"
		ret = subprocess.getoutput(cmd)
		headers = ret.split("\n")
		print("\tnumber of headers: ", len(headers))
		for hdr in headers:
			fields = hdr.split(" ")
			name = fields[0].replace (">", "")
			# print name
			if (name not in name2file):
				name2file[name] = []
			name2file[name].append(file)

	return name2file


###################
def make_filenames_table(cursor, db_name, table):
	check_and_drop_table(cursor, db_name, table)
	switch_to_db(cursor, db_name)
	qry = ""
	qry += f"CREATE TABLE  {table} ("
	qry += "     file_id int unsigned AUTO_INCREMENT, "
	qry += "  	 file_name text NOT NULL, "
	qry += "	 PRIMARY KEY (file_id) "
	qry += ") ENGINE=MyISAM"
	error_intolerant_search(cursor, qry)
	return


def store_filenames(cursor, name2file, table):
	file_ids = {}
	allnames = set()
	for names in name2file.values():
		allnames.update(names)
	for fnm in allnames:
		file_id = store_or_update(cursor, table, fixed_fields={'file_name':fnm},
									update_fields={}, primary_key='file_id')
		if len(file_id)>1:
			print(f"unexpected multiple table entries for file name {fnm}")
			exit()
		file_ids[fnm] = file_id[0]

	return file_ids


###################
def output_the_mapping(cursor, tmpfile, name2file, file_id):
	outf = open(tmpfile, "w")
	qry = "select seq_region_id, name from seq_region"
	for seq_region_id, region_name in hard_landing_search(cursor, qry):
		if region_name not in name2file:
			# print(f"region name {region_name} not found")
			file_ids = "not found"
		else:
			file_ids = ",".join([str(file_id[f]) for f in name2file[region_name]])
		outf.write(f"{seq_region_id}\t{file_ids}\n")

	outf.close()


###################
def make_mapping_table(cursor, db_name, table):
	check_and_drop_table(cursor, db_name, table)
	switch_to_db(cursor, db_name)
	qry = ""
	qry += f"CREATE TABLE  {table} ("
	qry += "     seq_region_id int unsigned, "
	qry += "     file_ids varchar(255), "
	qry += "	 PRIMARY KEY (seq_region_id) "
	qry += ") ENGINE=MyISAM"
	error_intolerant_search(cursor, qry)
	return


def slurp_in(db_name, table):
	# this fails:
	# "Loading local data is disabled; this must be enabled on both the client and server sides"
	# qry = f"load data local infile '{dirpath}/{tmpfile}' into table {db_name}.{table}"
	# error_intolerant_search(cursor, qry)
	# for this to work conf file must have the local_infile=1 line
	# and in mysql one has to run
	# mysql>  SET GLOBAL local_infile = 'ON';
	# check
	# mysql> SHOW GLOBAL VARIABLES LIKE 'local_infile';
	# not sure how thid variable gets flipped back to OFF - on server restart?
	# table name and tsv name must match
	cmd = f"mysqlimport --login-path=tcga --fields_escaped_by=\\\\ {db_name} -L {table}.tsv"
	subprocess.call(["bash","-c", cmd])
	print("out of slurping")

	return


####################################################
def main():

	db = connect_to_mysql(Config.mysql_conf_file)

	cursor = db.cursor()
	fasta_path = f"/storage/databases/ensembl-{Config.release_number}/fasta"

	[all_species, ensembl_db_name] = get_species (cursor)

	for species in all_species:

		time0 = time()

		print("\n========================\n"+species)
		dna_path = f"{fasta_path}/{species}/dna"
		if not os.path.exists(dna_path):
			print("problem:", dna_path, "not found")
			exit(1)

		# find the file each region is stored in
		name2file = collect_fasta_file_content(dna_path)

		# store each file and return the assigned id
		table = "file_names"
		make_filenames_table(cursor, ensembl_db_name[species], table)
		file_ids = store_filenames(cursor, name2file, table)

		table = "seq_region2file"
		# table name and tsv name must match
		tmpfile = f"{table}.tsv"
		output_the_mapping(cursor, tmpfile, name2file, file_ids)

		make_mapping_table(cursor, ensembl_db_name[species], table)
		slurp_in(ensembl_db_name[species], table)
		# cleanup
		os.remove(tmpfile)

		print("%50s stored in %.1f mins"%(species, (time()-time0)/60))


	cursor.close()
	db    .close()



####################################################
if __name__ == '__main__':
	main()
