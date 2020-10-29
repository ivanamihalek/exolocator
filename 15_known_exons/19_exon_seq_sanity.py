#!/usr/bin/python3 -u

from random import sample
from el_utils.el_specific import  *
from el_utils.processes import parallelize
from config import Config


def test_rough_number_correpondence(cursor, ensembl_db_name, all_species):
	for species in ['homo_sapiens'] + sample(all_species, 10):
	#for species in ['monodelphis_domestica']:
		species_db = ensembl_db_name[species]
		gene_ids = get_gene_ids(cursor, biotype='protein_coding', db_name=species_db)
		print()
		for gene_id in sample(gene_ids,10):
		#for gene_id in [8979]:
			exons = get_sorted_canonical_exons(cursor, species_db, gene_id)
			exon_ids_str = ",".join([str(exon.id) for exon in exons])
			qry = f"select count(*) from exon_seq where exon_id in ({exon_ids_str})"
			ret = error_intolerant_search(cursor, qry)
			number_of_exons_w_seq= ret[0][0] if ret else 0
			# warn = "" if number_of_exons==gene2exon_ct else "    <====="
			print(f"{gene_id}     {number_of_exons}    {number_of_exons_w_seq}    ")


#########################################
def main():


	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)

	test_rough_number_correpondence(cursor, ensembl_db_name, all_species)


	cursor.close()
	db    .close()


#########################################
if __name__ == '__main__':
	main()
