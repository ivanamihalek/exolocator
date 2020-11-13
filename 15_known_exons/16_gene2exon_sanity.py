#!/usr/bin/python3

from random import sample

from config import Config
from el_utils.ensembl import *


def gene_questions(cursor, species_db, gene_id):
	number_of_exons = len(get_ensembl_exons(cursor, species_db, gene_id))
	qry = f"select count(*) from {species_db}.gene2exon where gene_id={gene_id}"
	gene2exon_ct = hard_landing_search(cursor, qry)[0][0]
	warn = "" if number_of_exons==gene2exon_ct else "    <====="

	qry = f"select count(*) from {species_db}.gene2exon where gene_id={gene_id} and  is_canonical=1"
	gene2exon_canonical = hard_landing_search(cursor, qry)[0][0]

	print(f"gene id: {gene_id}      number_of_exons: {number_of_exons}", end="   ")
	print(f"gene2exon_ct: {gene2exon_ct}    canonical: {gene2exon_canonical}   {warn}")


def specific_test(cursor, species, ensembl_db_name, gene_id):
	species_db = ensembl_db_name[species]
	gene_ids = get_gene_ids(cursor, biotype='protein_coding', db_name=species_db)
	print()
	print(f"{species}  number of genes {len(gene_ids)}")
	gene_questions(cursor, species_db, gene_id)


def sampling_test(cursor, ensembl_db_name,all_species):
	for species in ['homo_sapiens'] + sample(all_species, 10):
		species_db = ensembl_db_name[species]
		gene_ids = get_gene_ids(cursor, biotype='protein_coding', db_name=species_db)
		print()
		print(f"{species}  number of genes {len(gene_ids)}")
		for gene_id in sample(gene_ids, 10):
			gene_questions(cursor, species_db, gene_id)


def exhaustive_test(cursor, ensembl_db_name,all_species):
	# I know some transcript are missing exon in mus caroli, paharus, and spretus
	for species in all_species:
		species_db = ensembl_db_name[species]
		gene_ids = get_gene_ids(cursor, biotype='protein_coding', db_name=species_db)
		print(f"{species}  number of genes {len(gene_ids)}")
		for gene_id in sample(gene_ids, 1000): # still not exhaustive - too boring
			gene_questions(cursor, species_db, gene_id)


def main():

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)

	#specific_test(cursor, 'monodelphis_domestica', ensembl_db_name, 8979)
	sampling_test(cursor, ensembl_db_name,all_species)
	#exhaustive_test(cursor, ensembl_db_name,all_species)

	cursor.close()
	db    .close()

	return


#########################################
if __name__ == '__main__':
	main()
