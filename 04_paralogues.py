#!/usr/bin/python3 -u
# -u flag forces the flushing ('u' for unboffered output)

from el_utils.ensembl   import *
from el_utils.processes import parallelize
from config import Config

#########################################
def store_paralogues (cursor_species, gene_id, paralogues):

	for para in paralogues:
		[para_stable_id, species, cognate_genome_db_id] = para
		para_gene_id = stable2gene (cursor_species, para_stable_id)

		fixed_fields  = {'gene_id': gene_id,
						 'cognate_genome_db_id': cognate_genome_db_id,
						 'cognate_gene_id': para_gene_id}

		update_fields = {'source': 'ensembl'}

		store_or_update(cursor_species, 'paralogue', fixed_fields, update_fields, primary_key='orth_pair_id')

#########################################
def collect_paralogues(species_list, db_info):

	[ensembl_db_name] = db_info

	db_species  = connect_to_mysql(Config.mysql_conf_file)
	cursor_species = db_species.cursor()
	search_db(cursor_species,"set autocommit=1")

	ensembl_compara_name = get_compara_name(cursor_species)
	db_compara     = connect_to_mysql(Config.mysql_conf_file)
	cursor_compara = db_compara.cursor()
	switch_to_db (cursor_compara, ensembl_compara_name)

	for species in species_list:
		print(species)
		switch_to_db (cursor_species,  ensembl_db_name[species])
		# it looks I cannot demand that the gene is known, because for many species
		# most of the genes still have 'predicted' status
		gene_list = get_gene_ids(cursor_species, biotype='protein_coding')
		ct = 0
		for gene_id in gene_list:
			ct += 1
			# find stable
			stable_id = gene2stable(cursor_species, gene_id=gene_id)
			# memebr id refers to entries in compara db
			member_id = stable2member(cursor_compara, stable_id)

			#print gene_id, stable_id, member_id
			if not ct%100:
				print("\t", species, ct , "out of ", len(gene_list))
			# find all paralogue pairs suggested for this gene
			ortho_type = 'within_species_paralog'
			paralogues = get_orthologues(cursor_compara, ortho_type, member_id)
			if not paralogues: continue
			store_paralogues (cursor_species, gene_id, paralogues)
		print(species, 'done')

	cursor_species.close()
	db_species.close()
	cursor_compara.close()
	db_compara.close()

#########################################
def main():

	number_of_chunks = 10

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species (cursor)
	cursor.close()
	db.close()

	parallelize (number_of_chunks, collect_paralogues, all_species, [ensembl_db_name])

	return True


#########################################
if __name__ == '__main__':
	main()
