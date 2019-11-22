#!/usr/bin/python3 -u
# -u flag forces the flushing ('u' for unboffered output)
import networkx as nx

from el_utils.ensembl   import *
from config import Config
from time import time

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
					#print("\tfusing clique; number of cliques:", len(cliques))
					#print("\t\t", len(s1), len(s2), intersection_size, len(new_clique))
					del cliques[j]
					del cliques[i]
					cliques.append(new_clique)
					done = False
					break
			if not done: break # break out of both loops
	return(cliques)

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

def flag_database(cursor_species, db_name, flag_name, msg):
	print("\t flagging database {}, message: {}".format(db_name, msg))
	source =  sys.argv[0].split("/")[-1]
	store_or_update(cursor_species, "exolocator_meta.flags",
					fixed_fields  = {"genome_db":db_name, "raised_by":source, "flag":flag_name},
					update_fields = {"comment":msg} )
	exit()

#########################################
def store_paralogue_cliques (cursor_species, cliques):
	table = "paralogue_groups"
	error_intolerant_search(cursor_species, "delete from %s"%table)
	error_intolerant_search(cursor_species, "alter table %s AUTO_INCREMENT = 1"%table)
	for clique in cliques:
		#store_without_checking(cursor_species, table, fields={"stable_ids":",".join(clique)})
		print("\n","#"*20)
		for ensid in clique:
			print(ensid, hard_landing_search(cursor_species, "select description from gene where stable_id='%s'"%ensid)[0][0])
	return

#########################################
def find_paralogue_cliques(cursor_compara, cursor_species, ensembl_db_name, species):
	cliques = {}
	
	return cliques

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
	all_species = ['homo_sapiens']
	#all_species = ['oryzias_latipes']
	for species in all_species:
		time0 = time()
		# how long would all-against-all blasting take?
		cliques = find_paralogue_cliques(cursor_compara, cursor_species, ensembl_db_name, species)
		print("%s done in %.1f mins"%(species, (time()-time0)/60))
		if not cliques: continue
		store_paralogue_cliques(cursor_species, cliques)

	cursor_species.close()
	db_species.close()

	cursor_compara.close()
	db_compara.close()


	return True


#########################################
if __name__ == '__main__':
	main()
