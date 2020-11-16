#!/usr/bin/python3 -u

from el_utils.processes import *
from el_utils.special_gene_sets import *

from config import Config
from itertools import combinations
import networkx as nx


#########################################
def mark_canonical (cursor, species_db, gene_id, exons):

	canonical_transcript_id = get_canonical_transcript_id(cursor, gene_id, species_db)
	if not canonical_transcript_id:
		print("canonical_transcript_id  not retrived for ",  gene_id)
		return False

	canonical_exon_ids = get_canonical_exon_ids(cursor, canonical_transcript_id, species_db)
	for exon in exons:
		# this has to be 0/1 and not T/F because we will be reading this into MySQL that does not know about T/F
		setattr(exon, "is_canonical", 0)  # we might not have this object until this point
		if exon.exon_id in canonical_exon_ids: exon.is_canonical = 1

	ret = get_canonical_coordinates(cursor, canonical_transcript_id, species_db)
	if not ret:
		return f"Error: no canonical_coordinates found for transcript {canonical_transcript_id} gene {gene_id}."
	[canonical_start_in_exon, canonical_start_exon_id, canonical_end_in_exon, canonical_end_exon_id] = ret

	for exon in exons:
		setattr(exon, "canon_transl_start", -1)
		setattr(exon, "canon_transl_end",  -1)
	start_found = False
	end_found   = False
	for exon in exons:
		# in translation [n.b: translation, not transcript] table canonical_start_in_exon and canonical_end_in_exon
		# refer to distance from the exon start - these are small numbers, like 15 or 138, starting from 1
		# from ensembl schema page, translation columns description
		# seq_start: 1-based offset into the relative coordinate system of start_exon_id (which is given in the next column)
		# seq_end: 1-based offset into the relative coordinate system of end_exon_id
		# and all that after the reverse direction is taken into the account
		# I find it extremely difficult to work with
		# I am storing the coding start position measured from the gene start, 0 offset, irrespective of reading direction
		if exon.exon_id == canonical_start_exon_id:
			start_found = True
			if exon.seq_region_strand > 0:
				exon.canon_transl_start = exon.start_in_gene + canonical_start_in_exon - 1
			else:
				exon.canon_transl_end = exon.end_in_gene - (canonical_start_in_exon - 1)
		if exon.exon_id == canonical_end_exon_id:
			end_found = True
			if exon.seq_region_strand > 0:
				exon.canon_transl_end = exon.start_in_gene + canonical_end_in_exon - 1
			else:
				exon.canon_transl_start = exon.end_in_gene - (canonical_end_in_exon - 1)

	# can somehting be canonical if we do not know where it starts
	if not start_found:
		return f"canonical translation start not found for {gene_id}"
	if not end_found:
		return f"canonical translation end not found for {gene_id}"

	return "ok"


#########################################
def fill_in_annotation_info (cursor, species_db, exons):

	for exon in exons:
		qry  = f"select {species_db}.transcript.analysis_id "
		qry += f"from {species_db}.transcript, {species_db}.exon_transcript "
		qry += "where exon_transcript.transcript_id = transcript.transcript_id "
		qry += f"and  exon_transcript.exon_id = {exon.exon_id} "

		rows = error_intolerant_search(cursor, qry)

		if not rows:
			exon.analysis_id = 0
		else:
			exon.analysis_id = rows[0][0]


#########################################
def get_exon_start(cursor, species_db, exon_id):

	qry  = f"select seq_region_start from {species_db}.exon  where exon_id = {exon_id}"
	rows = error_intolerant_search (cursor, qry)
	if not rows:
		print("start not found for ", exon_id)
		return None

	return rows[0][0]


#########################################
def get_exon_end(cursor, species_db, exon_id):

	qry  = f"select seq_region_end from {species_db}.exon where exon_id = {exon_id} "
	rows = error_intolerant_search (cursor, qry)
	if not rows:
		print("start not found for ", exon_id)
		return None

	return rows[0][0]


#########################################
def get_translated_region(cursor, species_db, gene_id):

	# get the region on the gene
	ret = get_gene_coordinates(cursor, species_db, gene_id)
	if ret:
		[gene_seq_id, gene_region_start, gene_region_end, gene_region_strand] = ret
	else:
		# in the current implementation this is never going to happen - hard landing if no gene coords
		print(f"region not retrieved for  {species_db} {gene_id}")
		return []

	transcript_ids = get_transcript_ids(cursor, species_db, gene_id)

	transl_region_start = gene_region_end
	transl_region_end   = gene_region_start

	for [transcript_id, tr_stable_id] in transcript_ids:

		qry  = "select seq_start, start_exon_id, seq_end, end_exon_id "
		qry += f" from {species_db}.translation where transcript_id={transcript_id}"
		rows = error_intolerant_search(cursor, qry)
		if not rows:
			# if the transcript is not protein coding, there will be no associated translation
			# there should eb an annotation for that, but we choose not to trust it
			#rows = search_db (cursor, qry, verbose=True)
			continue
		[exon_seq_start, start_exon_id, exon_seq_end, end_exon_id] = rows[0]

		if gene_region_strand > 0:
			start = {}
			start[start_exon_id] = get_exon_start(cursor, species_db, start_exon_id)
			start[end_exon_id]   = get_exon_start(cursor, species_db, end_exon_id)

			this_translation_region_start = start[start_exon_id] + exon_seq_start - 1
			this_translation_region_end   = start[end_exon_id]   + exon_seq_end   - 1

		else:
			end = {}

			end[start_exon_id] = get_exon_end (cursor, species_db, start_exon_id)
			end[end_exon_id]   = get_exon_end (cursor, species_db, end_exon_id)

			this_translation_region_start = end[end_exon_id]   - exon_seq_end   + 1
			this_translation_region_end   = end[start_exon_id] - exon_seq_start + 1

		if this_translation_region_start <= transl_region_start:
			transl_region_start = this_translation_region_start

		if this_translation_region_end >= transl_region_end:
			transl_region_end = this_translation_region_end

	return [transl_region_start, transl_region_end, gene_region_strand]


#########################################################
# resolve which exons are "master" among the ones
# which are havana, ensembl, predicted ....
#


#########################################
def find_master (exon_1, exon_2, is_ensembl, is_havana):

	master_exon    = None

	havana_exon    = None
	ensembl_exon   = None
	canonical_exon = None
	superset_exon  = None
	novel_exon     = None

	exon_start_1 = exon_1.start_in_gene
	exon_end_1   = exon_1.end_in_gene
	exon_start_2 = exon_2.start_in_gene
	exon_end_2   = exon_2.end_in_gene

	if exon_start_1 > exon_end_2 or exon_start_2 > exon_end_1:
		return None, None # Fully disjoint exons

	if exon_start_1 <= exon_start_2 < exon_end_2 <= exon_end_1:
		superset_exon = exon_1
	elif exon_start_2 <= exon_start_1 < exon_end_1 <= exon_end_2:
		superset_exon = exon_2

	if exon_1.is_canonical and not exon_2.is_canonical:
		canonical_exon = exon_1
	elif exon_2.is_canonical and not exon_1.is_canonical:
		canonical_exon = exon_2

	if is_havana[exon_1] and not is_havana[exon_2]:
		havana_exon = exon_1
	elif is_havana[exon_2] and not is_havana[exon_1]:
		havana_exon = exon_2

	if is_ensembl[exon_1] and not is_ensembl[exon_2]:
		ensembl_exon = exon_1
	elif is_ensembl[exon_2] and not is_ensembl[exon_1]:
		ensembl_exon = exon_2


	if canonical_exon is not None:
		# tough decision, but otherwise I get seqeunces which do not match ensembl's output
		# how can canonical exon not be havana - nobody ever loked at it?
		master_exon = canonical_exon
	elif havana_exon is not None:
		master_exon = havana_exon
	elif ensembl_exon  is not None:
		master_exon = ensembl_exon
	elif superset_exon is not None:
		master_exon = superset_exon
	elif novel_exon is not None:
		master_exon = novel_exon # ?

	if master_exon == exon_1:
		covered_exon = exon_2
	else:
		covered_exon = exon_1

	return master_exon, covered_exon


#########################################
def sort_out_covering_exons(cursor, exons):

	# havana is manually curated and gets priority
	is_ensembl = {}
	is_havana  = {}
	for exon in exons:
		logic_name = get_logic_name(cursor, exon.analysis_id)
		is_ensembl[exon] = ('ensembl' in logic_name)
		is_havana [exon] = ('havana'  in logic_name)
	# directed graph dg
	dg = nx.DiGraph()
	dg.add_nodes_from(exons)
	for exon1, exon2 in combinations(dg.nodes(),2):
		master, covered = find_master(exon1, exon2,is_ensembl, is_havana)
		if master is not None and covered is not None:
			dg.add_edge(master, covered)
	try:
		nx.find_cycle(dg)
		print("exon covering scheme results in a cyclic graph")
		exit()
	except:
		pass

	subtree_roots = [node for node in dg if dg.in_degree[node]==0]
	clusters = dict([(node, nx.descendants(dg,node)) for node in subtree_roots])
	# print("=========================")
	# for root, nodes in clusters.items():
	# 	print(root.exon_id, root.start_in_gene, root.end_in_gene, is_ensembl[root])
	# 	for n in nodes:
	# 		print("\t", n.exon_id, n.start_in_gene, n.end_in_gene, is_ensembl[n])

	for master_exon, covered_list in clusters.items():
		master_exon.covering_exon       = -1 # nobody's covering this guy
		for covered_exon in covered_list:
			covered_exon.covering_exon       = master_exon.exon_id

	return


#########################################
# check which exons are coding
def mark_coding (cursor, species_db, gene_id, exons):

	ret = get_translated_region(cursor, species_db, gene_id)
	if not ret: return False

	[transl_region_start,transl_region_end, strand] = ret

	translated_length = 0
	for exon in exons:

		# inocent until proven guilty
		exon.is_coding = 0
		# find related transcripts
		# if there is a related trascript, the thing is coding
		# (the default is 'not coding')
		qry  = f"select count(1) from {species_db}.exon_transcript "
		qry += f"where  exon_id = {exon.exon_id}"
		if not hard_landing_search(cursor, qry)[0][0]: continue

		# now need to check that it is within the coding region
		exon_start = get_exon_start(cursor, species_db, exon.exon_id)
		exon_end   = get_exon_end(cursor, species_db, exon.exon_id)

		translated_length += exon_end-exon_start+1

		if exon_end < transl_region_start or transl_region_end < exon_start:
			exon.is_coding = 0
		else:
			# there is _a_ translation that covers this exon
			# otherwise I could have two disjunct translations from
			# the saem gene, and a couple of exons in the middle that
			# are never translated - is that possible?
			exon.is_coding = 1


	return True


#########################################
# store into gene2exon table
# def store_exon (cursor, exon):
def format_tsv_line(cursor, exon):

	# note we are not storing the stable id here
	# (the field stable_id  does not exist in gene2exon table)
	# also, exon_seq_id field exists in gene2exon, but is apparently not used(?)

	fixed_fields  = {}
	update_fields = {}

	fixed_fields['exon_id']     = exon.exon_id
	# todo: define provenance codes somewhere
	fixed_fields['provenance']  = 1  # this is to distinguish ensembl exons from the exons we find or correct
	update_fields['gene_id']  = exon.gene_id
	update_fields['start_in_gene']      = exon.start_in_gene
	update_fields['end_in_gene']        = exon.end_in_gene
	update_fields['exon_seq_id']        = -1 # we'll fill this after we extract the exon seqs
	update_fields['strand']             = exon.seq_region_strand
	update_fields['phase']              = exon.phase
	update_fields['canon_transl_start'] = exon.canon_transl_start
	update_fields['canon_transl_end']   = exon.canon_transl_end
	update_fields['is_coding']          = exon.is_coding
	update_fields['is_canonical']       = exon.is_canonical
	update_fields['is_constitutive']    = exon.is_constitutive
	update_fields['covering_exon']      = exon.covering_exon
	update_fields['covering_provenance']  = 1
	update_fields['analysis_id']        = exon.analysis_id

	# if not store_or_update(cursor, 'gene2exon',
	# 						fixed_fields, update_fields,
	# 						primary_key = "gene2exon_id"):
	# 	print("failed storing exon ", exon.exon_id, "from gene",  exon.gene_id)
	all_fields = fixed_fields
	all_fields.update(update_fields)
	db_field_names = ["gene_id", "exon_id", "start_in_gene", "end_in_gene",
	                  "canon_transl_start", "canon_transl_end", "exon_seq_id", "strand", "phase",
	                  "provenance", "is_coding", "is_canonical", "is_constitutive", "covering_exon",
	                  "covering_provenance", "analysis_id"]
	for k,v in all_fields.items():
		if v is None: all_fields[k] = "\\N"

	tabbed_line = "\t".join([str(all_fields[field_name]) for field_name in db_field_names])
	return tabbed_line
	# abandoned: store_without_checking(cursor, 'gene2exon', all_fields)


#########################################
# fetch all exons for a given gene,
# and figure out whether they are canonical and/or coding
# which ones cover others, and what is the source of the annotation
def find_exon_info (cursor, gene_id, species_db):

	# print (gene_id, hgnc_symbol(cursor, gene_id), get_description(cursor, gene_id))

	# get all exons from the 'exon' table
	# here we enforce that they are "known", that is,
	# that they have a corresponding entry in transcript table
	exons = get_exons(cursor, species_db, gene_id, 'exon')
	gene_start = get_gene_start(cursor, species_db, gene_id)
	for exon in exons:
		exon.set_gene_id(gene_id)
		exon.set_within_gene_coords(gene_start)
	# mark the exons belonging to canonical transcript
	ret = mark_canonical(cursor, species_db, gene_id, exons)
	if ret != "ok": return ret

	# get annotation info
	fill_in_annotation_info(cursor, species_db, exons)

	# find covering exons
	sort_out_covering_exons(cursor, exons)

	# mark coding exons
	mark_coding(cursor, species_db, gene_id, exons)

	return exons


#########################################
# profile decorator is for the use with kernprof (a line profiler):
#  ./el_utils/kernprof.py -l this.py
# followed by
# python3 -m line_profiler this.py.lprof
# see here https://github.com/rkern/line_profiler#line-profiler
# the reason I am using local kernprof.py is that I don't know where pip
# installed its version (if anywhere)
#@profile
#########################################
def exons_for_gene(cursor, gene_id, species_db, count, outfile, logf):
	# find all exons associated with the gene id
	exons = find_exon_info(cursor, gene_id, species_db)

	if type(exons)==str:
		print(f"\t{species_db} {gene2stable (cursor, gene_id=gene_id)} {exons}", file=logf)
		return count
	elif not exons:
		print(f"\t{species_db} {gene2stable (cursor, gene_id=gene_id)} no exons found", file=logf)
		return count  # if I got to here in the pipeline this shouldn't happen

	# store into gene2exon table
	for exon in exons:
		count += 1
		tabbed_line = format_tsv_line(cursor, exon)
		outfile.write("%d\t%s\n"%(count, tabbed_line))

	return count


#########################################
def exons_for_species(species_list, db_info):

	[ensembl_db_name, outdir] = db_info
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	logf = open("log.{}.txt".format(get_process_id()),"w")

	for species in species_list:
		# load this file later with
		# sudo mysqlimport  --local <species db>  outfir/species/gene2exon.tsv
		# mysqlimport strips any extension and uses what's left as a table name
		# before you begin, do
		# mysql> SET GLOBAL local_infile = 1;
		os.makedirs(f"{outdir}/{species}", exist_ok=True)
		outfile = open(f"{outdir}/{species}/gene2exon.tsv", "w")

		logf.write(species+" started\n")
		gene_ids = get_gene_ids(cursor, biotype='protein_coding', db_name= ensembl_db_name[species])
		# gene_ids = [8979]
		count = 0
		time0 = time()
		for gene_id in gene_ids:
			if gene_ids.index(gene_id)%1000==0:
				pct_of_genes_processed = float(int(gene_ids.index(gene_id)) + 1)/len(gene_ids)*100
				print("%50s:  %5.1f%%    %ds" % (species, pct_of_genes_processed, time()-time0))
				time0 = time()
			count = exons_for_gene(cursor, gene_id, ensembl_db_name[species], count, outfile, logf)
		outfile.close()
		print(f"{species} done", file=logf)
	logf.close()
	cursor.close()
	db.close()

	return True


def check_species_done( all_species,  outdir):
	unprocessed_species = []
	for species in all_species:
		gene2ex_file = f"{outdir}/{species}/gene2exon.tsv"
		if not os.path.exists(gene2ex_file):
			unprocessed_species.append(species)

	return unprocessed_species


#########################################
def main():
	outdir = "raw_tables"
	os.makedirs(outdir, exist_ok=True)

	no_threads = 32
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	#all_species = ["mus_musculus"]
	#all_species.remove('homo_sapiens')

	cursor.close()
	db    .close()

	unprocessed_species = check_species_done(all_species, "raw_tables")

	parallelize(no_threads, exons_for_species, unprocessed_species, [ensembl_db_name, outdir])


#########################################
if __name__ == '__main__':
	main()

'''
In v 101 a bunch of canonical transcript coordinates were missing for mus caroli
not sure if I should worry about that - there ar 75 cases like that
exmaple MGP_CAROLIEiJ_G0027698 Ppp1cc protein phosphatase 1 catalytic subunit gamma
The smae for mus pahari and mus spretus
'''