#!/usr/bin/python3

import os

from config import Config

from el_utils.ensembl       import  *


####################################################
def store_seq_filenames (cursor, db_name, seq_region_name, file_names):
	qry = f"update {db_name}.seq_region set file_name='{file_names}' where name='{seq_region_name}'"
	error_intolerant_search(cursor, qry)
	return

####################################################
def main():

	db = connect_to_mysql(Config.mysql_conf_file)

	cursor = db.cursor()
	fasta_path = f"/storage/databases/ensembl-{Config.release_number}/fasta"

	[all_species, ensembl_db_name] = get_species (cursor)

	for species in all_species:

		print("\n========================\n"+species)
		dna_path = f"{fasta_path}/{species}/dna"
		if not os.path.exists(dna_path):
			print("problem:", dna_path, "not found")
			exit(1)

		fasta_files = []
		for r,d,files in os.walk(dna_path):
			for file in files:
				if not file[-3:] == ".fa":
					continue
				fasta_files.append(file)

		name2file = {}
		for file in fasta_files:
			print(dna_path, file)
			# fgreo is enough if not looking for regexp
			cmd = f"LC_ALL=C fgrep '>' {dna_path}/{file}"
			ret = subprocess.getoutput(cmd)
			headers = ret.split("\n")
			print("\tnumber of headers: ", len(headers))
			for hdr in headers:
				fields = hdr.split(" ")
				name = fields[0].replace (">", "")
				#print name
				if (name not in name2file):
					name2file[name] = []
				name2file[name].append(file)

		# TODO move this to a new table and read it in from ccsv
		switch_to_db(cursor, ensembl_db_name[species])
		time1 = time0 = time()
		ct = 0
		for seq_region_id, name in hard_landing_search(cursor, "select seq_region_id, name from seq_region"):
			if name not in name2file:
				print(f"name {name} not found in {ensembl_db_name[species]}.seq_region")
				exit()
			files = " ".join(name2file[name])
			qry = f"update seq_region set file_name='{files}' where seq_region_id={seq_region_id}"
			error_intolerant_search(cursor, qry)
			ct += 1
			if ct%1000==0:
				print(f"{ct} in   %.1f mins " %((time()-time1)/60))
				time1=time()
		print("%50s stored in %.1f mins"%(species, (time()-time0)/60))
	cursor.close()
	db    .close()



####################################################
if __name__ == '__main__':
	main()
