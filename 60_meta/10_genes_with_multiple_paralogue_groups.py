#!/usr/bin/python3 -u

from config import Config
from el_utils.ensembl import *

# https://github.com/networkx/networkx
import networkx as nx

#########################################
def find_paralogue_cliques(cursor_compara, cursor_species,  ensembl_db_name, species):
	switch_to_db (cursor_species, ensembl_db_name[species])
	# get taxon_id for the species
	qry = "select tax_id from exolocator_meta.species_names where species='%s'" % species
	taxon_id = hard_landing_search(cursor_species,qry)[0][0]
	print(species, taxon_id)
	# get species_tree_node_id from taxon_id
	qry =  "select * from species_tree_node where taxon_id=%d" % taxon_id
	species_tree_node_id = hard_landing_search(cursor_compara,qry)
	for row in species_tree_node_id:
		node_id = row[0]
		qry = "select count(*) from homology where species_tree_node_id=%d " % node_id
		qry += "and description='within_species_paralog'"
		count = hard_landing_search(cursor_compara,qry)[0][0]
		if count==0: continue
		# get all paralogue pairs for that species_tree_node_id
		qry = "select homology_id from homology where species_tree_node_id=%d " % node_id
		qry += "and description='within_species_paralog'"
		paralogous_pairs = {}
		for homology_id in [row[0] for row in hard_landing_search(cursor_compara,qry)]:
			qry2 = "select gene_member_id from homology_member where homology_id=%d" % homology_id
			paralogous_pairs[homology_id] = [row[0] for row in hard_landing_search(cursor_compara,qry2)]
			if len(paralogous_pairs[homology_id])!=2:
				print("pair does not seem to be a pair")
				print(homology_id, paralogous_pairs[homology_id])
				exit()
		print("\t node id:", node_id, "number of pairs:", len(paralogous_pairs))

		# turn pairs into graph
		graph = nx.Graph()
		graph.add_edges_from(paralogous_pairs.values())
		# find cliques for all pairs
		# nx.find_cliques returns iterator
		# each clique is given as a list of gene__member_ids
		cliques = list(nx.find_cliques(graph))
		cliques.sort(key=len,reverse=True)
		print("\t number of cliques:",  len(list(cliques)))
		# print("\t largest cliques")
		# for clk in cliques[:10]:
		# 	print(clk)
		# 	for gene_id in clk:
		# 		stable_id = gene_member2stable(cursor_compara, gene_id)
		# 		description = get_description(cursor_species, stable_id)
		# 		print("\t", gene_id, stable_id, description)

		# are any cliques intersecting? it looks like sometimes they are
		# for i in range(len(cliques)-1):
		# 	s1 = set(cliques[i])
		# 	for j in range(i+1, len(cliques)):
		# 		s2 = set(cliques[j])
		# 		if len(s1.intersection(s2))>5:
		# 			#print("intersect:", s1, s2, "intersection", s1.intersection(s2))
		# 			print(species, "intersection size:", len(s1), len(s2), len(s1.intersection(s2)))
		# 			#exit()


#########################################
def main():

	db_compara     = connect_to_mysql(Config.mysql_conf_file)
	cursor_compara = db_compara.cursor()
	ensembl_compara_name = get_compara_name(cursor_compara)
	switch_to_db (cursor_compara, ensembl_compara_name)

	db_species  = connect_to_mysql(Config.mysql_conf_file)
	cursor_species = db_species.cursor()
	# cursor_species checks for release number in exolocator_meta
	[all_species, ensembl_db_name] = get_species(cursor_species)

	all_species = ['cricetulus_griseus_chok1gshd']
	for species in all_species:
		find_paralogue_cliques(cursor_compara, cursor_species, ensembl_db_name, species)

	cursor_species.close()
	db_species.close()

	cursor_compara.close()
	db_compara.close()



	return True


#########################################
if __name__ == '__main__':
	main()
