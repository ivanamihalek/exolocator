#!/usr/bin/python3 -u

from config import Config
from el_utils.ensembl import *

# https://github.com/networkx/networkx
import networkx as nx

#########################################
def fuse_cliques(cliques):

	done = False
	while not done:
		done = True
		for i in range(len(cliques)-1):
			s1 = set(cliques[i])
			for j in range(i+1, len(cliques)):
				s2 = set(cliques[j])
				intersection_size = len(s1.intersection(s2))
				if intersection_size/len(s1)>0.5 or intersection_size/len(s2)>0.5:
					new_clique = list(s1.union(s2))
					print("\tfusing clique; number of cliques:", len(cliques))
					print("\t\t", len(s1), len(s2), intersection_size, len(new_clique))
					del cliques[j]
					del cliques[i]
					cliques.append(new_clique)
					done = False
					break
			if not done: break # break out of both loops


#########################################
def find_genome_db_id(cursor_compara, cursor_species):

	# genome_db has assembly version of the genome
	# assembly version of the genome can be found in meta table of each assembly database

	# working backwards:
	qry = "select meta_value  from meta where meta_key='assembly.default'"
	assembly_version = hard_landing_search(cursor_species, qry)[0][0]
	print(assembly_version)

	qry = "select genome_db_id from genome_db where assembly='%s'" % assembly_version
	genome_db_id  = hard_landing_search(cursor_compara, qry)[0][0]
	print(genome_db_id)

	return genome_db_id

#########################################
def find_paralogue_cliques(cursor_compara, cursor_species,  ensembl_db_name, species):
	print("#"*20, "\n", species)
	switch_to_db (cursor_species, ensembl_db_name[species])
	genome_db_id = find_genome_db_id(cursor_compara, cursor_species)

	# get taxon_id for the species
	qry = "select tax_id from exolocator_meta.species_names where species='%s'" % species
	taxon_id = hard_landing_search(cursor_species,qry)[0][0]
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
		ret =  hard_landing_search(cursor_compara,qry)
		print("\t number of homology entries:", len(ret))
		for homology_id in [row[0] for row in ret]:
			qry2  = "select hom.gene_member_id from homology_member as hom, gene_member as gen "
			qry2 += "where hom.homology_id=%d " % homology_id
			# make sure we are in the same assembly because for some species there might be several
			qry2 += "and hom.gene_member_id=gen.gene_member_id and gen.genome_db_id=%d" % genome_db_id
			print("\r%s"%qry2, end='')
			ret = error_intolerant_search(cursor_compara,qry2)
			if not ret or len(ret)==0:continue
			if len(ret)!=2:
				print("pair does not seem to be a pair")
				print(homology_id, ret)
				exit()

			paralogous_pairs[homology_id] = [row[0] for row in ret]


		print("\t node id:", node_id, "number of pairs:", len(paralogous_pairs))
		# networx can crash the whole system if the number of pairs is too big
		if len(paralogous_pairs)>10000:
			#flag_database(ensembl_db_name[species])
			continue

		# turn pairs into graph
		graph = nx.Graph()
		graph.add_edges_from(paralogous_pairs.values())
		# find cliques for all pairs
		# nx.find_cliques returns iterator
		# each clique is given as a list of gene__member_ids
		cliques = list(nx.find_cliques(graph))
		cliques.sort(key=len,reverse=True)
		print("\t number of cliques:",  len(list(cliques)))
		if len(list(cliques))>2000:
			#flag_database(ensembl_db_name[species])
			continue
		# print("\t largest cliques")
		# for clk in cliques[:10]:
		# 	print(clk)
		# 	for gene_id in clk:
		# 		stable_id = gene_member2stable(cursor_compara, gene_id)
		# 		description = get_description(cursor_species, stable_id)
		# 		print("\t", gene_id, stable_id, description)

		# are any cliques intersecting? it looks like sometimes they are
		# in some pahtological cases (e.g. cricetulus_griseus_chok1gshd)
		# it looks rather that a sizeable number of links is missing
		# TODO: this is to be revisited when I actually see the sequences
		# but for now take that two cliques tha overlap more than 50% are actually
		# the same clique with missing links
		fuse_cliques(cliques)


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

	#all_species = ['cricetulus_griseus_chok1gshd']
	#all_species = ['homo_sapiens']
	#all_species = ['oryzias_latipes']
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
