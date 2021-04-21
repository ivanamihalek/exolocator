#!/usr/bin/python3 -u

from el_utils.map import *
#########################################
verbose = False


#########################################
def store(cursor, maps, db_name):

	for map in maps:
		fixed_fields  = {}
		update_fields = {}
		fixed_fields ['exon_id']              = map.exon_id_1
		fixed_fields ['cognate_genome_db_id'] = species2genome_db_id(cursor, map.species_2)
		fixed_fields ['cognate_exon_id']      = map.exon_id_2
		update_fields['cigar_line']           = map.cigar_line
		update_fields['similarity']           = map.similarity
		update_fields['source']               = map.source
		#####
		switch_to_db(cursor, db_name)
		store_or_update(cursor, 'exon_map', fixed_fields, update_fields, primary_key='exon_map_id')

	return True


#########################################
def write(cursor, outfile, id_offset, maps):
	map_id = id_offset
	for map in maps:
		map_id += 1
		fields = [map_id, map.exon_id_1, map.exon_id_2, species2genome_db_id(cursor, map.species_2)]
		# the last two fields: msa_bitstring, warning
		fields += [map.cigar_line,  map.similarity, map.source, "\\N", "\\N" ]
		print("\t".join([str(f) for f in  fields]), file=outfile)


#########################################
def map_cleanup(cursor, rep_species_db, rep_species_exons):
	exon_id_str = ",".join([str(exon.exon_id) for exon in rep_species_exons])
	qry  = f"delete from {rep_species_db}.exon_map where exon_id = ({exon_id_str})"
	search_db(cursor, qry, verbose=False)
	return True


#########################################
def exons_sane(species, exons):
	if not exons:
		print(f"\t{species} no exon info")
		return False
	if type(exons) == str and 'error' in exons.lower():
		print(f"\t{species}  {exons}")
		return False
	return True


def gene_has_a_map (cursor, db_name, exons):
	has_a_map = False
	for exon in exons:
		if not exon.is_canonical or not exon.is_coding: continue
		maps = get_maps(cursor, db_name, exon.exon_id)
		if maps:
			has_a_map = True
			break
	return has_a_map


#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./el_utils/kernprof.py -l this.py
# followed by
# python3 -m line_profiler this.py.lprof
# @profile
def maps_for_gene(cursor, ensembl_db_name, representative_species, orthologue):
	if verbose: print(f"\tin maps_for_gene for {representative_species};\n\torthologues: {orthologue}")

	rep_species_db = ensembl_db_name[representative_species]
	rep_species_exons = get_sorted_canonical_exons(cursor, rep_species_db, orthologue[representative_species])
	if not exons_sane(representative_species, rep_species_exons):
		if verbose: print("\tfailed exons sanity test")
		return []
	# deleting the old stuff from the database -- do I need this if I am writing to a file?
	map_cleanup(cursor, ensembl_db_name[representative_species], rep_species_exons)

	relevant_exons_rep_species = find_relevant_exons_and_pepseqs(cursor, rep_species_exons, rep_species_db, ref_species=True)
	if not relevant_exons_rep_species:
		if verbose: print(f"\tno relevant exons found for {representative_species}")
		return []

	maps = []
	for ortho_species, ortho_gene_id in orthologue.items():
		if ortho_species == representative_species: continue
		if verbose: print(ortho_species, ortho_gene_id)

		ortho_exons = get_sorted_canonical_exons(cursor, ensembl_db_name[ortho_species], ortho_gene_id)
		if not exons_sane(ortho_species, ortho_exons): return []

		# maps are based on pairwise alignments of representative to other species in its group
		# multiple sequence alignements on exon-by-exon basis are produced in 17_ortho_exon_map_to_msa.py
		# reconstruction of full length multiple seqence alignments is  done only in
		# 30_db_migration/06_reconstruct_ortho_alnmts.py
		# make_maps(cursor, ensembl_db_name, rep_species, ortho_species, rep_exons, ortho_exons, verbose=False):
		maps.extend(make_maps(cursor, ensembl_db_name, representative_species, ortho_species,
							relevant_exons_rep_species, ortho_exons, verbose))

	return maps


#########################################
def human_genes_with_map_to_tax_group(cursor, ensembl_db_name, rep_species_db_id):
	# this is taking way too long
	# qry  = "select c1.cognate_gene_id, c1.cognate_genome_db_id, c2.cognate_gene_id, c2.cognate_genome_db_id "
	# qry += "from orthologues c1 left join orthologues c2 on c1.gene_id=c2.gene_id "
	# qry += f"where c1.cognate_genome_db_id='{cognate_genome_db_id}'"
	# TODO fix cases where the ortholgue is not presetn in the representative species
	# TODO  but might be present in other species
	qry  = f"select distinct(gene_id) from {ensembl_db_name['homo_sapiens']}.orthologues "
	qry += f"where cognate_genome_db_id='{rep_species_db_id}'"
	# the following can fail only if there are no orthologues at all
	ret  = error_intolerant_search(cursor, qry)
	if not ret:
		print(f"no orthologues found for cognate_genome_db_id {rep_species_db_id}")
		return []
	return [line[0] for line in ret]


#########################################
def maps_for_rep_species(representative_species, other_args):

	# tax_group_members = members of the taxonomy group represented by representative_species
	outdir  = other_args[0]
	os.makedirs(f"{outdir}/{representative_species}", exist_ok=True)
	outfile = open(f"{outdir}/{representative_species}/exon_map.tsv", "w")

	ensembl_db_name = other_args[1]
	tax_group_members = other_args[2:]
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	if verbose: print('checking indexing on orthologues table')
	create_index(cursor, ensembl_db_name['homo_sapiens'], "orthologues",  "gene_cognate_db_idx", ["gene_id", "cognate_genome_db_id"])
	if verbose: print('done with the index')

	# we are still chiefly interested in human - which genes from this tax_group are we considering
	cognate_genome_db_id = {}
	for species in [representative_species] + tax_group_members:
		cognate_genome_db_id[species] = species2genome_db_id(cursor, species)
	rep_species_db_id = cognate_genome_db_id[representative_species]
	if verbose: print(f"{representative_species} {cognate_genome_db_id}")

	if representative_species == "homo_sapiens":
		human_genes = get_gene_ids(cursor, biotype='protein_coding', db_name=ensembl_db_name[representative_species])
	else:  # we go through a bit of redirection
		human_genes = human_genes_with_map_to_tax_group(cursor, ensembl_db_name, rep_species_db_id)
	if verbose: print(f"number of human genes that map to this taxonomical group: {len(human_genes)}")

	orthologue = {}
	ct = 0
	time0 = time()
	id_offset = 0
	# NARS example - giving me some trouble in the beginning
	# for human_gene_id in [621232]:
	for human_gene_id in human_genes:
		switch_to_db(cursor,  ensembl_db_name['homo_sapiens'])
		if representative_species == "homo_sapiens":  # we are interested in all human genes
			qry  = "select cognate_gene_id, cognate_genome_db_id from orthologues  "
			qry += f"where gene_id={human_gene_id}"
		else:   # we are interested in the genes from reps species that map to human gene
			qry  = "select c2.cognate_gene_id, c2.cognate_genome_db_id "
			qry += "from orthologues c1 left join orthologues c2 on c1.gene_id=c2.gene_id "
			qry += f"where c1.gene_id={human_gene_id} and c1.cognate_genome_db_id='{rep_species_db_id}'"
		ret = error_intolerant_search(cursor, qry)
		if not ret: continue  # no human orthologues in the representative species
		if representative_species == "homo_sapiens": orthologue[representative_species] = human_gene_id
		for line in ret:
			[cognate_gene_id, belongs_to_genome_db_id] = line
			if belongs_to_genome_db_id not in cognate_genome_db_id.values(): continue # we're not interested today
			species = genome_db_id2species(cursor, belongs_to_genome_db_id)
			orthologue[species] = cognate_gene_id
		if verbose: print("\n", human_gene_id, gene2stable(cursor, human_gene_id), get_description(cursor, human_gene_id))
		maps = maps_for_gene(cursor, ensembl_db_name, representative_species, orthologue)
		if len(maps) == 0: continue
		write(cursor, outfile, id_offset, maps)
		id_offset += len(maps)

		ct += 1
		if ct % 100 == 0:
			time_for_last_ten = (time() - time0)/60
			print("maps for %6d genes written   --  last 100  in %5.1f min" % (ct, time_for_last_ten))
			time0 = time()

	cursor.close()
	db.close()

	outfile.close()
	return True


#########################################
def main():

	no_threads = 1
	outdir = "raw_tables"
	os.makedirs(outdir, exist_ok=True)

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	[all_species, ensembl_db_name] = get_species (cursor)
	members = {}
	qry = "select representative_species, members from exolocator_meta.taxonomy_groups "
	for rep_species, members_str in hard_landing_search(cursor, qry):
		members[rep_species] = members_str.split(",")
		if len(members[rep_species])==1: continue
		print(f"{rep_species} {len(members[rep_species])} ")
	print("================================")
	# TODO look into cases where the total coding length is not divisible by 3
	# TODO uncomment when I'm done with this
	# for species in all_species:
	# 	create_index(cursor, ensembl_db_name[species], "exon_seq", "exon_id_idx", ["exon_id"])

	#for rep_species in members.keys():
	for rep_species in ['carassius_auratus']:
		if len(members[rep_species])==1: continue
		print(f"{rep_species} {len(members[rep_species])} ")
		maps_for_rep_species(rep_species, [outdir, ensembl_db_name]+members[rep_species])

	cursor.close()
	db.close()

	# parallelize (no_threads, maps_for_gene_list, gene_list, [representative_species, ensembl_db_name])

	return True


#########################################
if __name__ == '__main__':
	main()

