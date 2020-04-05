#!/usr/bin/python3 -u

from el_utils.processes import *
from el_utils.special_gene_sets import *

from config import Config
from itertools import combinations
import networkx as nx


#########################################
def mark_canonical (cursor, gene_id, exons):

	canonical_transcript_id = get_canonical_transcript_id(cursor, gene_id)
	if not canonical_transcript_id:
		print("canonical_transcript_id  not retrived for ",  gene_id)
		return False

	canonical_exon_ids = get_canonical_exon_ids(cursor, canonical_transcript_id)
	for exon in exons:
		if exon.is_known and (exon.exon_id in canonical_exon_ids):
			exon.is_canonical = 1
		else:
			exon.is_canonical = 0
	[canonical_start_in_exon, canonical_start_exon_id,
	 canonical_end_in_exon, canonical_end_exon_id] = get_canonical_coordinates (cursor, canonical_transcript_id)

	for exon in exons:
		exon.canon_transl_start = -1
		exon.canon_transl_end = -1
	start_found = False
	end_found   = False
	for exon in exons:
		# in translation table canonical_start_in_exon and canonical_end_in_exon
		# refer to distance from the exon start - these are small numbers, like 15 or 138, starting from 1
		# and all thet after the reverse direction is taken into the account
		# I find it extremely difficult to work with
		# I am stroing the coding start position measured from the gene start, irrespective of reading direction
		if exon.exon_id == canonical_start_exon_id:
			start_found = True
			if exon.strand > 0:
				exon.canon_transl_start = exon.start_in_gene + canonical_start_in_exon - 1
			else:
				exon.canon_transl_end = exon.end_in_gene - (canonical_start_in_exon - 1)
		if exon.exon_id == canonical_end_exon_id:
			end_found = True
			if exon.strand > 0:
				exon.canon_transl_end = exon.start_in_gene +  canonical_end_in_exon -1
			else:
				exon.canon_transl_start = exon.end_in_gene - (canonical_end_in_exon - 1)

	# can somehting be canonical if we do not know where it starts
	if not start_found:
		print("canonical translation start not found for ", gene_id)
		return
	if not end_found:
		print("canonical translation end not found for ", gene_id)
		return
	# for each exon in canonical
	qry = "select exon_id  from exon_transcript where transcript_id= %d" % canonical_transcript_id
	rows = search_db (cursor, qry)
	if not rows:
		rows = search_db (cursor, qry, verbose=True)
		return

	canonical_ids = []
	for row in rows:
		canonical_ids.append(row[0])

	for exon in exons:
		exon.is_canonical = 0  # default
		if not exon.is_known:
			continue
		if exon.exon_id in canonical_ids:
			exon.is_canonical = 1

#########################################
def fill_in_annotation_info (cursor, gene_id, exons):

	for exon in exons:
		if exon.is_known:
			qry  = "select transcript.analysis_id "
			qry += " from transcript, exon_transcript "
			qry += " where exon_transcript.transcript_id = transcript.transcript_id "
			qry += " and  exon_transcript.exon_id = %d " % exon.exon_id
		else:
			qry  = " select prediction_transcript.analysis_id "
			qry += " from prediction_transcript, prediction_exon "
			qry += " where prediction_exon.prediction_transcript_id "
			qry += " = prediction_transcript.prediction_transcript_id and  "
			qry += " prediction_exon.prediction_exon_id= %d " %  exon.exon_id

		rows = error_intolerant_search(cursor, qry)

		if not rows:
			exon.analysis_id = 0
		else:
			exon.analysis_id = rows[0][0]

#########################################
#
#
# Note: the original function went something like this:
#
#def get_transcript_ids(cursor, gene_id, species):
#
#    rows = []
#    if ( species=='homo_sapiens'):
#        qry    = "SELECT transcript_id  FROM transcript "
#        qry   += " WHERE gene_id=%d AND status = 'known' AND biotype='protein_coding' "  \
#            %  gene_id
#        rows   = search_db (cursor, qry, verbose=False)
#
#    return  transcript_ids
#
#
# however, for human almost 2000   protein coding genes have canonical transcript annotated as 'putative' or 'novel'
# check:
# select transcript.status  from transcript, gene where gene.biotype='protein_coding' \\
# and gene.canonical_transcript_id = transcript.transcript_id and not transcript.status = 'KNOWN';
# 
# just go by the version in ensembl.py,
# that simply returns all transcript ids for a gene 
# (presumably we have cheked elsehwere that the gene is protein coding)
#########################################
def get_exon_start(cursor, exon_id):

	qry  = "select seq_region_start from exon "
	qry += "where exon_id = %d " % exon_id
	rows = search_db (cursor, qry)
	if (not rows or 'Error' in rows[0]):
		print("start not found for ", exon_id)
		return None

	return rows[0][0]


#########################################
def get_exon_end(cursor, exon_id):
	qry  = "select seq_region_end from exon "
	qry += "where exon_id = %d " % exon_id
	rows = search_db (cursor, qry)
	if (not rows or 'Error' in rows[0]):
		print("start not found for ", exon_id)
		return None

	return rows[0][0]


#########################################
def get_translated_region(cursor, gene_id, species):

	# get the region on the gene
	ret = get_gene_region (cursor, gene_id)
	if ret:
		[gene_seq_id,gene_region_start, gene_region_end,
		 gene_region_strand] = ret
	else:
		print("region not retrieved for ", species, gene_id, species)
		return []

	transcript_ids = get_transcript_ids(cursor, gene_id)

	transl_region_start = gene_region_end
	transl_region_end   = gene_region_start

	for [transcript_id, tr_stable_id] in transcript_ids:

		qry  = "SELECT seq_start, start_exon_id, seq_end, end_exon_id "
		qry += " FROM translation WHERE transcript_id=%d"  %  transcript_id
		rows = search_db(cursor, qry)
		if (not rows):
			# if the transcript is not protein coding, there will be no associated translation
			# there should eb an annotation for that, but we choose not to trust it
			#rows = search_db (cursor, qry, verbose=True)
			continue
		exon_seq_start = rows[0][0]
		start_exon_id  = rows[0][1]
		exon_seq_end   = rows[0][2]
		end_exon_id    = rows[0][3]

		if (gene_region_strand > 0):
			start = {}
			start[start_exon_id] = get_exon_start(cursor, start_exon_id)
			start[end_exon_id]   = get_exon_start(cursor, end_exon_id)

			this_translation_region_start = start[start_exon_id] + exon_seq_start - 1
			this_translation_region_end   = start[end_exon_id]   + exon_seq_end   - 1

		else:
			end = {}

			end[start_exon_id] = get_exon_end (cursor, start_exon_id)
			end[end_exon_id]   = get_exon_end (cursor, end_exon_id)

			this_translation_region_start = end[end_exon_id]   - exon_seq_end   + 1
			this_translation_region_end   = end[start_exon_id] - exon_seq_start + 1

		if (this_translation_region_start <= transl_region_start):
			transl_region_start = this_translation_region_start

		if (this_translation_region_end >= transl_region_end):
			transl_region_end = this_translation_region_end

	return [transl_region_start, transl_region_end, gene_region_strand]


#########################################################
# resolve which exons are "master" among the ones
# which are havana, ensembl, predicted ....
#


#########################################
def find_master (cursor, exon_1, exon_2, is_ensembl, is_havana):

	master_exon    = None
	covered_exon   = None

	havana_exon    = None
	ensembl_exon   = None
	canonical_exon = None
	superset_exon  = None
	known_exon     = None
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

	if exon_1.is_known and not exon_2.is_known:
		known_exon = exon_1
	elif exon_2.is_known and not exon_1.is_known:
		known_exon = exon_2

	if exon_1.is_known<2 and not exon_2.is_known < 2:
		novel_exon = exon_2
	elif exon_2.is_known<2  and not exon_1.is_known < 2:
		novel_exon = exon_1


	if canonical_exon is not None:
		# tough decision, but otherwise I get seqeunces which do not match ensembl's output
		# how can canonical exon not be havana - nobody ever loked at it?
		master_exon = canonical_exon
	elif havana_exon is not None:
		master_exon = havana_exon
	elif ensembl_exon  is not None:
		master_exon = ensembl_exon
	elif known_exon is not None:
		master_exon = known_exon
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
def sort_out_covering_exons (cursor, exons):

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
		master, covered = find_master(cursor, exon1,exon2,is_ensembl,is_havana)
		if master is not None and covered is not None:
			dg.add_edge(master,covered)
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
		master_exon.covering_exon_known = -1 # formal
		for covered_exon in covered_list:
			covered_exon.covering_exon       = master_exon.exon_id
			covered_exon.covering_exon_known = master_exon.is_known

	return


#########################################
# check which exons are coding
def mark_coding (cursor, gene_id, species, exons):

	ret = get_translated_region(cursor, gene_id, species)
	if not ret:
		return False

	[transl_region_start,transl_region_end, strand] = ret

	translated_length = 0
	for exon in exons:

		exon.is_coding = 0
		if (exon.is_known):
			# inocent until proven guilty
			exon.is_coding = 0
			# find related transcripts
			# if there is a related trascript, the thing is coding
			# (the default is 'not coding')
			qry  = "select count(1) from exon_transcript "
			qry += " where  exon_id = %d " % exon.exon_id
			rows = search_db (cursor, qry)
			if (not rows[0][0]):
				continue

			# now need to check that it is within the coding region
			exon_start = get_exon_start(cursor, exon.exon_id)
			exon_end   = get_exon_end(cursor, exon.exon_id)

			translated_length += exon_end-exon_start+1

			if exon_end < transl_region_start or transl_region_end < exon_start:
				exon.is_coding = 0
			else:
				# there is _a_ translation that covers this exon
				# otherwise I could have two disjunct translations from
				# the saem gene, and a couple of exons in the middle that
				# are never translated - is that possible?
				exon.is_coding = 1

		else: # exons belongs to a  predicted transcript
			# == we don't know if it is coding or not
			exon.is_coding = 0
			# if it is covered by a coding exon, it is coding

	return True


#########################################
# fetch all exons for a given gene,
# and figure out whether they are canonical and/or coding
# which ones cover others, and what is the source of the annotation
def find_exon_info (cursor, gene_id, species):

	# get all exons from the 'exon' table
	# here we enforce that they are "known", that is,
	# that they have a corresponding entry in transcript table
	exons = get_exons(cursor, gene_id, species, 'exon')
	# assorted predicted exons
	if not species == 'homo_sapiens':
		exons += get_exons(cursor, gene_id, species, 'prediction_exon')
		exons += get_exons(cursor, gene_id, species, 'sw_exon')
		exons += get_exons(cursor, gene_id, species, 'usearch_exon')

	# mark the exons belonging to canonical transcript
	mark_canonical(cursor, gene_id, exons)

	# get annotation info
	fill_in_annotation_info(cursor, gene_id, exons)

	# find covering exons
	sort_out_covering_exons(cursor, exons)

	# mark coding exons
	mark_coding(cursor, gene_id, species, exons)

	return exons


#########################################
# store into gene2exon table
# def store_exon (cursor, exon):
def format_tsv_line(cursor, exon):

	# note we are not storing the stable id here
	# (the field stable_id  does not exist in gene2exon table)
	# also, exon_seq_id field exists in gene2exon, but is apparently not used(?)

	fixed_fields  = {}
	update_fields = {}

	fixed_fields['exon_id']  = exon.exon_id
	# I should get rid of this 'is known' thing - ensembl dropped it
	fixed_fields['is_known'] = exon.is_known

	update_fields['gene_id']  = exon.gene_id
	update_fields['start_in_gene']      = exon.start_in_gene
	update_fields['end_in_gene']        = exon.end_in_gene
	update_fields['exon_seq_id']        = exon.exon_seq_id
	update_fields['strand']             = exon.strand
	update_fields['phase']              = exon.phase
	update_fields['canon_transl_start'] = exon.canon_transl_start
	update_fields['canon_transl_end']   = exon.canon_transl_end
	update_fields['is_coding']          = exon.is_coding
	update_fields['is_canonical']       = exon.is_canonical
	update_fields['is_constitutive']    = exon.is_constitutive
	update_fields['covering_exon']      = exon.covering_exon
	update_fields['covering_is_known']  = exon.covering_exon_known
	update_fields['analysis_id']        = exon.analysis_id

	# if not store_or_update(cursor, 'gene2exon',
	# 						fixed_fields, update_fields,
	# 						primary_key = "gene2exon_id"):
	# 	print("failed storing exon ", exon.exon_id, "from gene",  exon.gene_id)
	all_fields = fixed_fields
	all_fields.update(update_fields)
	db_field_names = ["gene_id", "exon_id", "start_in_gene", "end_in_gene",
	                 "canon_transl_start", "canon_transl_end", "exon_seq_id", "strand", "phase",
	                  "is_known", "is_coding", "is_canonical", "is_constitutive", "covering_exon",
	                  "covering_is_known", "analysis_id"]
	for k,v in all_fields.items():
		if v is None: all_fields[k] = "\\N"

	tabbed_line = "\t".join([str(all_fields[field_name]) for field_name in db_field_names])
	return tabbed_line
	# abandoned: store_without_checking(cursor, 'gene2exon', all_fields)

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
def gene2exon_all(species_list, db_info):

	ensembl_db_name,  = db_info
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	logf = open("log.{}.txt".format(get_process_id()),"w")


	for species in species_list:
		# load this file later with
		# sudo mysqlimport  --local <species db>  <"gene2exon.%s"%species>
		# mysqlimport strips any extension and uses what's left as a table name
		# before you begin, do
		# mysql> SET GLOBAL local_infile = 1;
		outfile = open("gene2exon.%s"%species, "w")
		switch_to_db (cursor,  ensembl_db_name[species])
		gene_ids = get_gene_ids (cursor, biotype='protein_coding')
		count = 0
		time0 = time()
		for gene_id in gene_ids:
			if gene_ids.index(gene_id)%10==0:
				print("%50s:  %5.1f%%    %ds" %  \
					(species, float(int(gene_ids.index(gene_id)) +1 )/len(gene_ids)*100,
				      time()-time0))
				time0 = time()

			#print (gene_id, get_description(cursor, gene_id))
			# find all exons associated with the gene id
			exons = find_exon_info (cursor, gene_id, species)
			if not exons:
				print(gene2stable (cursor, gene_id = gene_id), " no exons found (%s)" % species)
				continue  # if I got to here in the pipeline this shouldn't happen
			# store into gene2exon table
			for exon in exons:
				count += 1
				tabbed_line = format_tsv_line(cursor, exon)
				outfile.write("%d\t%s\n"%(count, tabbed_line))
		logf.write(species+"\n")
		outfile.close()
	logf.close()
	cursor.close()
	db.close()

	return True

#########################################
def gene2exon_orthologues(gene_list, db_info):

	ensembl_db_name,  = db_info
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	for gene_id in gene_list:

		switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
		orthologues  = get_orthos (cursor, gene_id, 'orthologue') # get_orthos changes the db pointer
		switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
		orthologues += get_orthos (cursor, gene_id, 'unresolved_ortho')

		for [ortho_gene_id, ortho_species] in orthologues:

			switch_to_db (cursor, ensembl_db_name[ortho_species])

			# find all exons associated with the gene id
			exons = find_exon_info (cursor, ortho_gene_id, ortho_species)
			if (not exons):
				print(gene2stable (cursor, gene_id = ortho_gene_id), " no exons found for", ortho_species)
				print()
				continue  # if I got to here in the pipeline this shouldn't happen


			exons.sort(key=lambda exon: exon.start_in_gene)
			# store into gene2exon table
			for exon in exons:
				store_exon (cursor, exon)

		print("progress:  %8.3f " %  ( float( int(gene_list.index(gene_id)) +1 )/len(gene_list)))

	cursor.close()
	db.close()

	return True



#########################################
def main():

	no_threads = 1
	special    = None

	if len(sys.argv) > 1 and  len(sys.argv)<3  or len(sys.argv) >= 2 and sys.argv[1]=="-h":
		print("usage: %s <set name> <number of threads>" % sys.argv[0])
		exit(1) # after usage statment
	elif len(sys.argv)==3:
		special = sys.argv[1].lower()
		if special == 'none': special = None
		no_threads = int(sys.argv[2])

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	all_species = ['homo_sapiens']

	print('=======================================')
	print(sys.argv[0])
	gene_list = []
	if special:
		print("using", special, "set")
		# if special == 'complement':
		# 	gene_list = get_complement_ids(cursor, ensembl_db_name, cfg)
		# else:
		# 	gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special)

	cursor.close()
	db    .close()

	# two versions of the main loop:
	# 1) over all species, and all genes in each species
	if not special:
		parallelize(no_threads, gene2exon_all, all_species,  [ensembl_db_name])
	# 2) over orthologues for a given list of genes
	else:
		parallelize(no_threads, gene2exon_orthologues, gene_list,  [ensembl_db_name])


#########################################
if __name__ == '__main__':
	main()

