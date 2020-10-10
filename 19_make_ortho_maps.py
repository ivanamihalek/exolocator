#!/usr/bin/python3 -u

from config import Config
from el_utils.ensembl import *
from el_utils.tree import *
from el_utils.processes import  parallelize


#########################################
verbose = True

#########################################
def store (cursor, maps, ensembl_db_name):

	for map in maps:
		fixed_fields  = {}
		update_fields = {}
		fixed_fields ['exon_id']              = map.exon_id_1
		fixed_fields ['exon_known']           = map.exon_known_1
		fixed_fields ['cognate_genome_db_id'] = species2genome_db_id(cursor, map.species_2)
		fixed_fields ['cognate_exon_id']      = map.exon_id_2
		fixed_fields ['cognate_exon_known']   = map.exon_known_2
		update_fields['cigar_line']           = map.cigar_line
		update_fields['similarity']           = map.similarity
		update_fields['source']               = map.source
		#####
		switch_to_db(cursor,ensembl_db_name['homo_sapiens'])
		store_or_update (cursor, 'exon_map', fixed_fields, update_fields)

	return True

#########################################
def  map_cleanup (cursor, ensembl_db_name, human_exons):

	switch_to_db(cursor,ensembl_db_name['homo_sapiens'])
	for exon in human_exons:
		qry  = "delete from exon_map where exon_id = %d " % exon.exon_id
		qry += " and exon_known = %d " % exon.is_known
		rows = search_db (cursor, qry, verbose=False)

	return True

#########################################
def gene_has_a_map (cursor, ensembl_db_name, human_exons):

	has_a_map = False
	# for human_exon in human_exons:
	#     if ( not human_exon.is_canonical or  not human_exon.is_coding): continue
	#     maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
	#     if maps:
	#         has_a_map = True
	#         break

	return has_a_map


#########################################
def maps_for_gene_list(gene_list, db_info):

	[all_species, ensembl_db_name, tree] = db_info
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	missing_exon_info = 0
	missing_seq_info  = 0
	ct                = 0
	no_maps           = 0

	#######################################################
	#
	gene_list = [616358] # gnao1
	for gene_id in gene_list:

		ct += 1
		switch_to_db(cursor,  ensembl_db_name['homo_sapiens'])
		if verbose: print(gene_id, gene2stable(cursor, gene_id), get_description(cursor, gene_id))
		#
		# # get _all_ exons
		human_exons = gene2exon_list(cursor, gene_id)
		if not human_exons:
			print(f"no exons found for {gene_id}")
			continue
		#
		# # get rid of the old maps
		# map_cleanup (cursor, ensembl_db_name, human_exons)
		#
		# # human as its own orthologue - let's be systematic
		# maps = self_maps (cursor, ensembl_db_name, human_exons)
		# store (cursor, maps, ensembl_db_name)
		#
		# # get_reliable_orthos will return unresolved orthologues only if no unequivocal 1-to-1 ortho is found
		# orthologues = get_reliable_orthos(cursor, ensembl_db_name, gene_id)
		#
		# for [ortho_gene_id, ortho_species] in  orthologues:
		# 	ortho_exons = gene2exon_list(cursor, ortho_gene_id, db_name=ensembl_db_name[ortho_species] )
		# 	# in the first round of building the database there will be no 'novel' exons
		# 	# but we put this here just in case we want to rerun the script anyway
		# 	# (the script 21_make_novel_exon_maps.py rebuilds the maps only for the species for which
		# 	# there are some novel exons - this one builds the maps from scratch)
		# 	ortho_exons += get_novel_exons (cursor, ortho_gene_id, 'sw_exon')
		# 	ortho_exons += get_novel_exons (cursor, ortho_gene_id, 'usearch_exon')
		#
		# 	if not ortho_exons:
		# 		missing_exon_info += 1
		# 		print("\t", ortho_species, "no exon info")
		# 		continue
		# 	# maps are based on pairwise alignements of human to other species
		# 	# multiple seqence alignements on exon-by-exon basis are produced in 17_ortho_exon_map_to_msa.py
		# 	# reconstruction of full length multiple seqence alignments is  done only in
		# 	# 30_db_migration/06_reconstruct_ortho_alnmts.py
		# 	maps = make_maps (cursor, ensembl_db_name,  cfg, acg, ortho_species, human_exons, ortho_exons, verbose)
		# 	if not maps:
		# 		missing_seq_info += 1
		# 		print("\t", ortho_species, "no maps")
		# 		continue
		#
		# 	no_maps += len(maps)
		# 	store (cursor, maps, ensembl_db_name)
		#
		# if not ct%100:
		# 	datastring = io.StringIO()
		# 	print("processed ", ct, "genes,  out of ", len(gene_list), "  ", end=' ', file=datastring)
		# 	print(no_maps, " maps;  no_exon_info: ", missing_exon_info, "no_seq_info:", missing_seq_info, file=datastring)
		# 	print(datastring.getvalue())
	print("gene list done")
	cursor.close()
	db.close()

	return True


#########################################
def main():

	no_threads = 1
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()


	[all_species, ensembl_db_name] = get_species (cursor)
	tree = species_tree(cursor, all_species)
	# print(tree.nhx_string()) # https://phylo.io
	switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
	gene_list = get_gene_ids (cursor, biotype='protein_coding')


	cursor.close()
	db.close()

	parallelize (no_threads, maps_for_gene_list, gene_list, [all_species, ensembl_db_name, tree])

	return True

#########################################
if __name__ == '__main__':
	main()

