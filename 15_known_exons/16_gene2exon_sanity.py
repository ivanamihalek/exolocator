#!/usr/bin/python3

from random import sample

from config import Config
from el_utils.ensembl import *


def sampling_test(cursor, ensembl_db_name,all_species):
	for species in ['homo_sapiens'] + sample(all_species, 10):
		species_db = ensembl_db_name[species]
		gene_ids = get_gene_ids(cursor, biotype='protein_coding', db_name=species_db)
		print()
		print(f"{species}  number of genes {len(gene_ids)}")
		for gene_id in sample(gene_ids, 10):
			number_of_exons = len(get_ensembl_exons(cursor, species_db, gene_id))
			qry = f"select count(*) from gene2exon where gene_id={gene_id}"
			gene2exon_ct = hard_landing_search(cursor, qry)[0][0]
			warn = "" if number_of_exons==gene2exon_ct else "    <====="
			print(f"{gene_id}     {number_of_exons}    {gene2exon_ct}    {warn}")


def exhaustive_test(cursor, ensembl_db_name,all_species):
	# I know some transcript are missing exon in mus caroli, paharus, and spretus
	for species in all_species:
		species_db = ensembl_db_name[species]
		gene_ids = get_gene_ids(cursor, biotype='protein_coding', db_name=species_db)
		print(f"{species}  number of genes {len(gene_ids)}")
		for gene_id in sample(gene_ids, 1000): # still not exhaustive - too boring
			number_of_exons = len(get_ensembl_exons(cursor, species_db, gene_id))
			qry = f"select count(*) from gene2exon where gene_id={gene_id}"
			gene2exon_ct = hard_landing_search(cursor, qry)[0][0]
			if number_of_exons!=gene2exon_ct :
				print(f"\t{gene_id}     {number_of_exons}    {gene2exon_ct}  ")
				#exit()


def main():

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	#all_species = ["mus_musculus"]
	#all_species.remove('homo_sapiens')

	#sampling_test(cursor, ensembl_db_name,all_species)
	exhaustive_test(cursor, ensembl_db_name,all_species)

	cursor.close()
	db    .close()

	return


#########################################
if __name__ == '__main__':
	main()
