#!/usr/bin/python3 -u
# -u flag forces the flushing ('u' for unboffered output)

from el_utils.ensembl   import *
from el_utils.processes import parallelize
from config import Config
from time import time


#########################################
def store_paralogues (cursor, ensembl_db_name, paralogues):

	id_string = ",".join(paralogues)
	print(id_string)
	exit()
	group_id = store_without_checking(cursor, 'paralogue_groups', {'stable_ids':id_string}, database=ensembl_db_name)
	if group_id<0: exit()
	for stable_id in paralogues:
		qry = "update gene set paralogue_group_id=%d where stable_id='%s'"%(group_id,stable_id)
		error_intolerant_search(cursor, qry)


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
		time_species_started = time()
		print(species)
		switch_to_db (cursor_species,  ensembl_db_name[species])
		# it looks I cannot demand that the gene is known, because for many species
		# most of the genes still have 'predicted' status
		gene_ids = get_gene_ids(cursor_species, biotype='protein_coding')
		gene_list = dict([(gene_id,False) for gene_id in gene_ids])
		ct = 0
		time0 = time()
		for gene_id, seen in gene_list.items():
			ct += 1
			if seen: continue
			# find stable
			stable_id = gene2stable(cursor_species, gene_id=gene_id)
			# memebr id refers to entries in compara db
			member_id = stable2member(cursor_compara, stable_id)

			#print gene_id, stable_id, member_id
			if not ct%100:
				time1 = time()
				print("\t", species, ct , "out of", len(gene_list), "in %.1f min"%((time1-time0)/60))
				time0 = time1
			# find all paralogue pairs suggested for this gene
			ortho_type = 'within_species_paralog'
			paralogues = get_orthologues(cursor_compara, ortho_type, member_id)
			if not paralogues: continue
			store_paralogues (cursor_species, gene_id, [p[0] for p in paralogues])
		print(species, "done in %.1f min"%((time()-time_species_started)/60))

	cursor_species.close()
	db_species.close()
	cursor_compara.close()
	db_compara.close()

#########################################
def main():

	number_of_chunks = 1 # this is not working (deadlock?)

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
