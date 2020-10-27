#!/usr/bin/python3

from config import Config

from el_utils.ensembl import get_species
from el_utils.ncbi    import *
from el_utils.tree import Tree, species_tree


#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)

	# do we have the tree?
	ret = error_intolerant_search(cursor, "select value from exolocator_meta.taxonomy where name = 'species_tree'")
	if not ret:
		print("building the species tree ...")
		tree = species_tree(cursor, all_species)
		print("                     tree done.")
		print("storing")
		qry = f"insert into  exolocator_meta.taxonomy  (name,value) values ('species_tree','{tree.nhx_string()}') "
		error_intolerant_search(cursor, qry)
	else:
		print("reading the species tree ...")
		tree = Tree(ret[0][0])
		print("                     tree done.")

	trivial_name = {'Archosauria' : "birds_and_crocs",
	                    'Testudines'  : "turtles",
	                    'Lepidosauria' : "lizards_and_snakes",
	                    'Eutheria' : "mammals",
	                    'Marsupialia': "marsupials",
	                    'ornithorhynchus_anatinus':  "platypus",
	                    'Anura': "frogs",
	                    'latimeria_chalumnae': "coelacanth",
	                    'Euteleosteomorpha':"euteleosts",
	                    'Otomorpha':"otomorpha",
	                    'Osteoglossiformes':"osteoglossiforms",  # ray-finned fish
	                    'lepisosteus_oculatus':"spotted_gar",
	                    'erpetoichthys_calabaricus':"snakefish",  # more ray-finned fish
	                    'callorhinchus_milii':"elephant_shark",  # Australian ghostshark or elephant shark
	                    'Cyclostomata':"nightmare stuff"
	                    }
	species_subtrees = trivial_name.keys()
	# which species in the subtree has the best annotation so far?

	switch_to_db(cursor, "exolocator_meta")
	for tax_group in species_subtrees:
		node = tree.get_node(tax_group)
		print()
		print(tax_group)
		number_of_genes = {}
		number_of_transcripts = {}

		group_species = [node.name] if node.is_leaf else node.subtree_leafs()
		for species in group_species:
			qry = f"select count(*) from {ensembl_db_name[species]}.transcript"
			number_of_transcripts[species] = hard_landing_search(cursor,qry)[0][0]
			qry = f"select count(*) from {ensembl_db_name[species]}.gene"
			number_of_genes[species] = hard_landing_search(cursor,qry)[0][0]
		# we are using the reported number of transcript as an ad hoc measure of reliability of the genome annotation
		sorted_species_in_the_group = sorted(group_species, key = lambda s:number_of_transcripts[s], reverse=True)
		for species in sorted_species_in_the_group:
			strformat = "%50s:      transcripts:  %6d      genes:  %6d"
			print(strformat % (species, number_of_transcripts[species],number_of_genes[species]))
		fixed_fields  = {'name':tax_group}
		update_fields = {'trivial_name':trivial_name[tax_group],
		                 'representative_species':sorted_species_in_the_group[0],
						 'members': ",".join(sorted_species_in_the_group[1:])}
		print(fixed_fields)
		print(update_fields)
		print()
		store_or_update(cursor, "taxonomy_groups", fixed_fields, update_fields)


	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
