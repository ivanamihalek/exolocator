#!/usr/bin/python3 -u

from config import Config
from el_utils.ensembl import *
from el_utils.tree import *
from el_utils.processes import  parallelize
from el_utils.el_specific import  *
from el_utils.map import *
#########################################
verbose = True


#########################################
def gene_has_a_map (cursor, db_name, exons):
	has_a_map = False
	for exon in exons:
		if not exon['is_canonical'] or  not exon['is_coding']: continue
		maps = get_maps(cursor, db_name, exon['exon_id'])
		if maps:
			has_a_map = True
			break
	return has_a_map


def exons_sane(species, exons):
	if not exons:
		print(f"\t{species} no exon info")
		return False
	if type(exons)==str and 'error' in exons.lower():
		print(f"\t{species}  {exons}")
		return False
	return True



#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./el_utils/kernprof.py -l this.py
# followed by
# python3 -m line_profiler this.py.lprof
# @profile
def maps_for_rep_species(rep_species, other_args):

	# tax_group_members = members of the taxonomy group represented by representative_species
	ensembl_db_name = other_args[0]
	tax_group_members = other_args[1:]
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor,  ensembl_db_name[rep_species])
	qry = "select distinct(g.gene_id) from exon_map e left join gene2exon g on e.exon_id = g.exon_id"
	genes_with_maps = [ret[0] for ret in hard_landing_search(cursor, qry)]
	print(f"genes with maps: {len(genes_with_maps)}")

	cursor.close()
	db.close()

	return True


#########################################
def main():

	no_threads = 1
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()


	[all_species, ensembl_db_name] = get_species (cursor)
	# eate_index(cursor, db_name,                         table, index_name, columns, verbose=False)
	# check whther index exists
	print('checking indexing on orthologues table')
	create_index(cursor, ensembl_db_name['homo_sapiens'], "orthologues",  "gene_cognate_db_idx", ["gene_id", "cognate_genome_db_id"])
	print('done with the index')

	members = {}
	qry = "select representative_species, members from exolocator_meta.taxonomy_groups "
	for rep_species, members_str in hard_landing_search(cursor, qry):
		members[rep_species] = members_str.split(",")
		if len(members[rep_species])==1: continue
		print(f"{rep_species} {len(members[rep_species])} ")
	print("================================")
	# exit()

	for rep_species in ['monodelphis_domestica']:
	#for rep_species in ['anas_platyrhynchos_platyrhynchos']:
		maps_for_rep_species(rep_species, [ensembl_db_name]+members[rep_species])

	cursor.close()
	db.close()

	# parallelize (no_threads, maps_for_gene_list, gene_list, [representative_species, ensembl_db_name])

	return True

#########################################
if __name__ == '__main__':
	main()

