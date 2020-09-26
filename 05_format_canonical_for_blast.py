#!/usr/bin/python3 -u


from config import Config

from el_utils.ensembl import  *
from el_utils.processes import *

formatting_tool = "/usr/bin/makeblastdb"


def format_for_blast(species):
	fasta_path = "/storage/databases/ensembl-{}/fasta/{}/pep".format(Config.release_number, species)
	canonical_fasta = f"{fasta_path}/canonical_peptides.fa"
	for dep in [fasta_path, canonical_fasta]:
		if not os.path.exists(dep):
			print(f"{dep} not found")
			return
	cmd = f"{formatting_tool} -in {canonical_fasta} -dbtype prot -parse_seqids"
	subprocess.call(["bash","-c", cmd])


def format_for_blast_loop(species_chunk, other_args):
	for species in species_chunk:
		format_for_blast(species)


####################################################
def main():

	if not os.path.exists(formatting_tool):
		print(f"{formatting_tool} not found")
		exit()

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	cursor.close()
	db.close()

	number_of_chunks = 8

	print(f"number of species: {len(all_species)}, number of chunks: {number_of_chunks}")
	parallelize(number_of_chunks, format_for_blast_loop, all_species, [])



#####################################################
if __name__=="__main__":
	main()
