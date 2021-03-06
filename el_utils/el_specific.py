import MySQLdb, subprocess, re, os, time
import inspect
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna
from bitstring import Bits

from  .ensembl import  *
from  .utils   import  *

from config import Config

def store_problem(cursor, db_name, gene_id, description):
	ret = error_intolerant_search(cursor, f"select * from {db_name}.problems where gene_id={gene_id}")
	if ret:
		[gene_id, prev_description] = ret[0]
		if description not in prev_description:
			qry = f"update {db_name}.problems set description='{prev_description} | {description}' where gene_id={gene_id}"
			error_intolerant_search(cursor, qry)
	else:
		store_without_checking(cursor, 'problems', {'gene_id':gene_id, 'description':description}, database=db_name)


def get_file_name(cursor, seq_region_id):
	# for human I do not expect that a gene would be split across different regions,
	# but for other species not so sure
	qry =  f"select file_ids from seq_region2file where seq_region_id={seq_region_id}"
	file_ids = hard_landing_search(cursor,qry)[0][0]
	qry = f"select file_name from file_names where file_id in ({file_ids})"
	file_names = [ret[0] for ret in hard_landing_search(cursor, qry)]

	if len(file_names)>1:
		print(f"multiple file names found for seq_region_id {seq_region_id}")
		exit()
	return file_names[0]


#########################################
def extract_region_from_fasta(fasta_full_path, seq_region_name, seq_region_start, seq_region_end):

	tmpfile = f"tmp.{os.getpid()}.fasta"
	if os.path.exists(tmpfile): os.remove(tmpfile)
	# %s == "no header line"
	cmd  = f"{Config.blastdbcmd} -db {fasta_full_path} -dbtype nucl -entry {seq_region_name} "
	cmd += f"-range {seq_region_start}-{seq_region_end} -out {tmpfile} -outfmt %s "
	subprocess.call(["bash", "-c", cmd])
	if not os.path.exists(tmpfile):
		print(f"{tmpfile} not produced")
		exit()
	with open(tmpfile) as inf:
		inseq = inf.read().replace("\n", "")
	os.remove(tmpfile)
	return inseq


#########################################
def get_gene_dna(cursor, species, db_name, gene_id):
	switch_to_db(cursor, db_name)
	qry =  f"select seq_region_id, seq_region_start, seq_region_end from gene where gene_id={gene_id}"
	[gene_region_id, gene_region_start, gene_region_end] = hard_landing_search(cursor, qry)[0]

	qry =  f"select name from seq_region where seq_region_id={gene_region_id}"
	seq_region_name = hard_landing_search(cursor, qry)[0][0]

	file_name = get_file_name(cursor, gene_region_id)

	# convention over configuration :}
	fasta_full_path = f"{Config.fasta_repo}/{species}/dna/{file_name}"
	if not os.path.exists(fasta_full_path):
		print(f"{fasta_full_path} not found")
		exit()

	return extract_region_from_fasta(fasta_full_path, seq_region_name, gene_region_start, gene_region_end)


def get_sorted_canonical_exons(cursor, db_name, gene_id):

	column_names = get_column_names(cursor, db_name, "gene2exon")
	exons = []
	ret = error_intolerant_search(cursor, f"select * from {db_name}.gene2exon where gene_id={gene_id}")
	if not ret:
		return []
	for line in ret:
		exon = Exon()
		exon.fields2attributes(column_names, line)
		if not exon.is_canonical: continue
		exons.append(exon)
	sorted_exons = sorted(exons,key=lambda x: x.start_in_gene)
	# sanity check
	total_length = 0
	for exon in sorted_exons:
		#print(exon)
		reading = exon.is_coding
		if reading:
			total_length += exon.end() - exon.start()  + 1

	if total_length%3 != 0:
		# take care of that later - Esmebl has quite a number of those
		return f"Error: total coding length for canonical exons {total_length} not divisible by 3"

	return sorted_exons


#################
def exons2cdna(gene_region_dna, sorted_exons):
	reverse = None
	for exon in sorted_exons:
		if reverse is None:
			reverse = (exon.strand<0)
		elif reverse != (exon.strand<0):
			print(f"gene {exon.gene_id}, exon {exon.exon_id} strand mismatch")
			exit()

	cdna= ""
	for exon in sorted_exons:
		reading = exon.is_coding
		if reading:
			cdna += gene_region_dna[exon.start():exon.end()+1]

	if reverse: cdna = str(Seq(cdna, generic_dna).reverse_complement())

	return cdna



#
# #########################################
# def check_afa_age (cfg, directory, stable_id, max_days=3):
#
#     afa_age = "old"
#
#     afa_fnm  = "{0}/{1}.afa".format(directory, stable_id)
#
#     if (os.path.exists(afa_fnm) and os.path.getsize(afa_fnm) > 0 ):
#         time_modified = os.path.getmtime(afa_fnm)
#         number_of_days_since_modified = (time.time() - time_modified)/(60*60*24)
#         if number_of_days_since_modified < max_days:
#             #print "\t %s last modified %s. Moving on." % (stable_id, time.ctime(os.path.getmtime(afa_fnm) ))
#             afa_age  = "new"
#     return afa_age
#
#
#
# #########################################
# def strip_stop(pepseq):
#     if (not pepseq or len(pepseq)==0):
#         return pepseq
#     if ( pepseq[-1] == '*'):
#         pepseq = pepseq[:-1]
#     return pepseq
#
# #########################################
# def  transl_reconstruct (cursor,  gene_id, gene_seq, canonical_coding_exons,
#                          is_mitochondrial, verbose = False):
#
#     """
#     Given tha dna sequence, gene information an the list of its canonical exons, reconstruct canonical translation.
#
#     Pay attention to whether the gene is mitochondrial.
#     Here in particular we can catch the 'false stop codons' corresponding to selenocysteines.
#     (Because they are stored by their position
#     in the translation.)
#
#     """
#
#     canonical_exon_pepseq = {}
#     translated_seq = ""
#
#     [can_transl_start_exon, can_transl_start_position,
#      can_transl_end_exon, can_transl_end_position] = canonical_transl_info (cursor, gene_id)
#
#
#     # do we have any selenocysteines by any chance
#     selenoC_pos = get_selenocysteines (cursor,  gene_id)
#
#     carry = ""
#     ok_so_far = True
#     # sanity checking
#     for exon in canonical_coding_exons:
#         #find exon sequence within the gene
#         start = exon.start_in_gene
#         if (exon is canonical_coding_exons[0]):
#             if ( not exon.exon_id == can_transl_start_exon ):
#                 print(" error start cantransl:  gene_id ",  gene_id, end=' ')
#                 print(" exon_id ", exon.exon_id, " canon: ", can_transl_start_exon)
#                 return [{}, ""]
#                 #exit (1)
#             start +=  exon.canon_transl_start
#
#         if ( exon is canonical_coding_exons[-1]):
#             if ( not exon.exon_id == can_transl_end_exon ):
#                 print(" error end cantransl:  gene_id ",  gene_id, end=' ')
#                 print(" exon_id ", exon.exon_id, " canon: ", can_transl_end_exon)
#                 return [{}, ""]
#                 #exit (1)
#             end = exon.start_in_gene + exon.canon_transl_end
#         else:
#             end = exon.end_in_gene
#
#         if (not exon.phase == -1 and not exon.phase == len(carry)):
#             #print "Houston we have a problem: exon phase =", exon.phase,
#             #print " the length of carry =", len(carry),
#             #print " (gene_id %d, exon_id %d) " % (gene_id, exon.exon_id)
#             if (  exon.phase ):
#                 start += 3-exon.phase
#             carry = ""
#
#
#         exon_seq     =  gene_seq[ start: end+1]
#         exon_seq_for_transl_purposes = carry + exon_seq
#
#         remainder    = len(exon_seq_for_transl_purposes)%3
#         if ( remainder == 0 ):
#             carry = ""
#         elif (remainder == 1 ):
#             carry    = exon_seq_for_transl_purposes[-1:]
#             exon_seq_for_transl_purposes  = exon_seq_for_transl_purposes[:-1]
#         else:
#             carry    = exon_seq_for_transl_purposes[-2:]
#             exon_seq_for_transl_purposes = exon_seq_for_transl_purposes[:-2]
#
#         dnaseq = Seq (exon_seq_for_transl_purposes, generic_dna)
#         pepseq = dnaseq.translate()
#
#         # turn to the corresponding BioPython object
#         if ( is_mitochondrial ):
#             pepseq = dnaseq.translate(table="Vertebrate Mitochondrial")
#         else:
#             pepseq = dnaseq.translate()
#
#         # replace stop codons from selenoC positions, if there are such
#         if (selenoC_pos):
#             for pos_in_full_length_translation in selenoC_pos:
#                 pos_in_exon_translation = pos_in_full_length_translation-len(translated_seq)
#                 if pos_in_exon_translation<0 or pos_in_exon_translation>len(pepseq):
#                     continue
#                 tempseq = pepseq.tomutable()
#                 tempseq[pos_in_exon_translation] = 'U'
#                 pepseq  = tempseq.toseq()
#
#         # strip the last stop codon only
#         if ( exon is canonical_coding_exons[-1]):
#             pepseq = strip_stop(pepseq)
#
#         pepseq0 =  pepseq
#         if verbose:
#             print("phase 0", pepseq0)
#
#         # if there are stil stop codons we'll give another shot
#         # to the possibility that it is mitochondrial, (and we ddin't know it)
#         # after that we cry foul
#         if ( '*' in pepseq):
#             pepseq = dnaseq.translate(table="Vertebrate Mitochondrial")
#
#         # strip the last stop codon only
#         if ( exon is canonical_coding_exons[-1]):
#             pepseq = strip_stop(pepseq)
#
#             # some further  desperate measures
#         ok_so_far = True
#         if ( '*' in pepseq):
#             ok_so_far = False
#             dnaseq = Seq (exon_seq_for_transl_purposes[1:], generic_dna)
#             pepseq = dnaseq.translate()
#             pepseq = strip_stop(pepseq)
#             pepseq1 =  pepseq
#             if verbose:
#                 print("phase 1", pepseq1)
#         else:
#             ok_so_far = True
#
#         if (not ok_so_far and '*' in pepseq):
#             dnaseq = Seq (exon_seq_for_transl_purposes[2:], generic_dna)
#             pepseq = dnaseq.translate()
#             pepseq = strip_stop(pepseq)
#             pepseq2 =  pepseq
#             if verbose:
#                 print("phase 2", pepseq2)
#         else:
#             ok_so_far = True
#
#         if (not ok_so_far and  '*' in pepseq):
#             if verbose:
#                 print("Error: stop codon ")
#         else:
#             ok_so_far = True
#
#         if ( not ok_so_far):
#             return [{}, ""]
#         translated_seq += pepseq # I need the seq in selenoC - to decide
#         # where the position of  U should be
#         canonical_exon_pepseq[exon.exon_id] = pepseq.tostring()
#
#     return [canonical_exon_pepseq,translated_seq]
#
# #########################################
# def compare_seqs (canonical_translation, translated_seq, verbose=False):
#
#     """
#     Return true if the two input sequences differ in no more than two positions.
#     """
#     comparison_ok = True
#
#     while (len(translated_seq) and translated_seq[0] == 'X'):
#         translated_seq = translated_seq[1:]
#
#     difference = len(translated_seq) - len(canonical_translation)
#     if ( abs(difference) > 3):
#         comparison_ok = False
#         if verbose:
#             print()
#             print(">canon")
#             print(canonical_translation)
#             print(">exons")
#             print(translated_seq)
#             print()
#     else:
#         diff  =  0
#         start = -1
#         for i in range(len(translated_seq)):
#             if ( i >= len(canonical_translation)):
#                 break
#             if (not translated_seq[i] ==  canonical_translation[i]):
#                 diff += 1
#                 if start < 0:
#                     start = i
#         if (diff > 2):
#             comparison_ok = False
#             if verbose:
#                 print()
#                 print(">canon")
#                 print(canonical_translation)
#                 print(">exons")
#                 print(translated_seq)
#                 print(translated_seq[start], canonical_translation[start])
#                 print("nuber of  diff sites: ", diff, " starting from ", start)
#                 print()
#
#     return comparison_ok
#
#
#
# #########################################
# def  get_gene_seq (acg, cursor, gene_id, species, verbose = False):
#
#     """
#     Given gene_id, return dna region which reproduces the correct canonical translation.
#     """
#     null = ["",{}, "", "", None, None]
#
#     # i'm not quite clear why Ensembl is doing this, but sometimes we need the alternative
#     # region - ("PATCH" deposited as tte right sequence, but its missing most of the gene)
#     # so first establish whether it is the case: find canonical translation.
#     canonical_translation  = get_canonical_transl (acg, cursor, gene_id, species)
#     if not canonical_translation: return null
#
#     #########################################
#     # which file should we be looking in, which sequence, from where to where
#     primary_seq_info = get_primary_seq_info (cursor, gene_id, species)
#     if not primary_seq_info:
#         return null
#     [seq_name, file_names, seq_region_start, seq_region_end,
#      seq_region_strand, is_mitochondrial] = primary_seq_info
#     # find all canonical exons associated with the gene id
#     canonical_coding_exons = get_canonical_exons (cursor, gene_id)
#     # extract raw gene  region TODO - store the information about which
#     # file_name we ended up using, and the start and end in that region
#     [gene_seq, file_name] = extract_gene_seq (acg, species, seq_name, file_names, seq_region_strand,
#                                               seq_region_start, seq_region_end)
#     # reconstruct the translation from the raw gene_seq and exon boundaries
#     [canonical_exon_pepseq,translated_seq] = transl_reconstruct (cursor, gene_id, gene_seq,
#                                                                  canonical_coding_exons, is_mitochondrial)
#     if (translated_seq):
#         # compare the two sequences and cry foul if they are not the same:
#         comparison_ok = compare_seqs (canonical_translation, translated_seq)
#     else:
#         comparison_ok = False
#     # if we succefully translated the exons, and came up with the same answer
#     # as the canonical translation, we are done here
#     if (comparison_ok):
#         return [gene_seq, canonical_exon_pepseq, file_name, seq_name, seq_region_start, seq_region_end]
#     if verbose:
#         print("Using primary seq info: failed comparison with canonical sequence.")
#         print("primary_seq_info: ")
#         print(primary_seq_info)
#         print("canonical translation:")
#         print(canonical_translation)
#         print("translated:")
#         print(translated_seq)
#         print("canonical exons:")
#         print(canonical_exon_pepseq)
#
#     #########################################
#     # otherwise repeat the procedure with the alternative seq info:
#     alt_seq_info = get_alt_seq_info (cursor, gene_id, species)
#     if (not alt_seq_info):
#         return null
#     [seq_name, file_names, seq_region_start, seq_region_end,
#      seq_region_strand, is_mitochondrial] = alt_seq_info
#     # find all canonical exons associated with the gene id
#     canonical_coding_exons = get_canonical_exons (cursor, gene_id)
#     # extract raw gene  region
#     [gene_seq, file_name] = extract_gene_seq(acg, species, seq_name, file_names, seq_region_strand,
#                                              seq_region_start, seq_region_end)
#     # reconstruct the translation from the raw gene_seq and exon boundaries
#     [canonical_exon_pepseq,translated_seq] = transl_reconstruct (cursor, gene_id, gene_seq, canonical_coding_exons,
#                                                                  is_mitochondrial)
#     if (translated_seq):
#         # compare the two sequences and cry foul if they are not the same:
#         comparison_ok = compare_seqs (canonical_translation, translated_seq)
#     else:
#         comparison_ok = False
#     # if we succefully translated the exons, and came up with the same answer
#     # as the canonical translation, we are done here
#     if (comparison_ok):
#         return [gene_seq, canonical_exon_pepseq, file_name, seq_name, seq_region_start, seq_region_end]
#
#     if verbose:
#         print()
#         print("Using alt seq info: failed comparison with canonical sequence.")
#         print("alt seq info:")
#         print(alt_seq_info)
#         print("canonical translation:")
#         print(canonical_translation)
#         print("translated:")
#         print(translated_seq)
#         print("canonical exons:")
#         print(canonical_exon_pepseq)
#
#     return null
#
#

# #########################################
def get_orthos(cursor, ref_species, other_species, ensembl_db_name, gene_id):
	ortho_gene_id = {}
	qry  = f"select cognate_gene_id, cognate_genome_db_id from {ensembl_db_name[ref_species]}.orthologues "
	qry += f"where gene_id={gene_id}"
	for line in error_intolerant_search(cursor, qry):
		[cognate_gene_id, cognate_genome_db_id] = line
		qry = f"select db_name from exolocator_meta.db_names where genome_db_id={cognate_genome_db_id}"
		db_name = hard_landing_search(cursor, qry)[0][0]
		#stable_transl_id = gene2stable_canon_transl_id(cursor, cognate_gene_id, db_name)
		species = db_name.split("core")[0].rstrip("_")
		if species not in other_species: continue
		ortho_gene_id[species]=cognate_gene_id
	return ortho_gene_id


# #########################################
# def get_reliable_orthos(cursor, ensembl_db_name, gene_id):
#
#     all_orthologues = []
#
#     # one2one   orthologues
#     switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
#     known_orthologues      = get_orthos (cursor, gene_id, 'orthologue')
#     # not-clear orthologues
#     switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
#     unresolved_orthologues = get_orthos (cursor, gene_id, 'unresolved_ortho')
#
#     # get rid of the unresolved orthologues if a resolved orthologue for the species already exists
#     species_with_known_orthologues = []
#     for  [ortho_gene_id, ortho_species] in known_orthologues:
#         species_with_known_orthologues.append(ortho_species)
#     all_orthologues = known_orthologues
#
#     for  [ortho_gene_id, ortho_species] in unresolved_orthologues:
#         if ortho_species in species_with_known_orthologues: continue
#         all_orthologues.append( [ortho_gene_id, ortho_species] )
#
#     # for each ortho check that it is not a pseudogene
#     all_protein_coding_orthologues = []
#     for  [ortho_gene_id, ortho_species] in all_orthologues:
#         biotype = get_gene_biotype(cursor, ortho_gene_id, ensembl_db_name[ortho_species]);
#         if ( biotype == 'protein_coding'):
#             all_protein_coding_orthologues.append( [ortho_gene_id, ortho_species] )
#         #else:
#         #print ortho_gene_id, ortho_species, biotype
#
#     return all_protein_coding_orthologues
#
#
#
# #########################################
# def sort_names (sorted_species, alignment):
#
#     sorted_names = []
#     for species in sorted_species:
#         for seq_name  in list(alignment.keys()):
#
#             if seq_name[-1].isdigit():
#                 aux = seq_name.split("_")
#                 base_name = "_".join(aux[:-1])
#
#                 if base_name [-1].isdigit():
#                     aux = base_name.split("_")
#                     base_name = "_".join(aux[:-1])
#             else:
#                 base_name = seq_name
#             if (species == base_name):
#                 sorted_names.append(seq_name)
#     return sorted_names
#
# #########################################
# def align_nucseq_by_pepseq (aligned_pepseq, nucseq):
#
#     # coding dna sequence, by assumption:
#     cds = nucseq
#     translated_cds = Seq(cds).translate().tostring()
#     if not len(aligned_pepseq.replace('-','')) == len(translated_cds):
#         print(aligned_pepseq.replace('-',''))
#         print("in align_nucseq_by_pepseq():  length mismatch: ", end=' ')
#         print(len(aligned_pepseq.replace('-','')), len(translated_cds))
#         return ""
#     codon = iter(map(''.join, list(zip(*[iter(nucseq)]*3))))
#     #aligned_nucseq = ''.join(('---' if c=='-' else next(codon) for c in aligned_pepseq))
#     aligned_nucseq = ''
#     for c in aligned_pepseq:
#         if c == '-': aligned_nucseq += '---'
#         else:        aligned_nucseq += next(codon)
#     return aligned_nucseq
#
# #########################################
# def expand_pepseq (aligned_pep_sequence, exon_seqs, flank_length):
#
#     dna_aln_expanded = ""
#
#     # check if this is a padding seqeunce:
#     if not exon_seqs or  not aligned_pep_sequence.replace('-',''):
#         dna_aln_expanded = '-'*(2*flank_length+3*len(aligned_pep_sequence) )
#         return dna_aln_expanded
#
#     [pepseq, pepseq_transl_start,
#      pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs
#     if not pepseq or pepseq_transl_start is None or pepseq_transl_end is None:
#         return ""
#     # coding dna sequence:
#     cds = dna_seq[pepseq_transl_start:pepseq_transl_end]
#     if not cds:
#         return ""
#
#     aligned_nucseq  = align_nucseq_by_pepseq(aligned_pep_sequence, cds)
#     if not aligned_nucseq:
#         print(aligned_pep_sequence)
#         print(Seq(cds).translate().tostring())
#         print(cds)
#         print(pepseq_transl_start, pepseq_transl_end," reconstruction failed")
#         return ""
#
#     effective_left_flank  = ""
#     effective_right_flank = ""
#     #######
#     effective_left_flank  = left_flank
#     if pepseq_transl_start>0:
#         effective_left_flank += dna_seq[:pepseq_transl_start]
#     if len(effective_left_flank) > flank_length:
#         effective_left_flank = effective_left_flank[-flank_length:]
#     effective_left_flank = effective_left_flank.lower()
#
#     if 1:
#         #######
#         effective_right_flank = right_flank
#         delta = len(dna_seq)-pepseq_transl_end
#         if delta>0:
#             effective_right_flank = dna_seq[-delta:]+effective_right_flank
#         if len(effective_right_flank) > flank_length:
#             effective_right_flank = effective_right_flank[:flank_length]
#         effective_right_flank = effective_right_flank.lower()
#
#     #######
#     # pad the flanking seqs to the needed length
#     effective_left_flank  = effective_left_flank.rjust (flank_length, '-')
#     effective_right_flank = effective_right_flank.ljust(flank_length, '-')
#
#     dna_aln_expanded = effective_left_flank + aligned_nucseq + effective_right_flank
#
#     return dna_aln_expanded
#
#
#
# #########################################
# def check_has_sw_exons (cursor, ensembl_db_name, human_exon_id, human_exon_known, minsim):
#
#     has_sw_exons = False
#
#     # find all other exons that map to the human exon
#     maps    = get_maps(cursor, ensembl_db_name, human_exon_id, human_exon_known)
#     maps    = [m for m in maps if not m.exon_id_2 is None]
#     maps_sw = [m for m in maps if m.source=='sw_sharp' and m.similarity >minsim]
#
#     if maps_sw:
#         has_sw_exons = True
#
#     return has_sw_exons
#
# #########################################
# def make_exon_alignment(cursor, ensembl_db_name, human_exon_id, human_exon_known, mitochondrial,
#                         min_similarity,  flank_length, first_human_exon = True):
#
#     sequence_pep = {}
#     sequence_dna = {}
#     shortest_l = -1 # Uninitialized  leading padding length
#     shortest_r = -1 # Uninitialized trailing padding length
#
#     pep_aln_length = 0
#     dna_aln_length = 0
#     # find all other exons that map to the human exon
#     maps    = get_maps(cursor, ensembl_db_name, human_exon_id, human_exon_known)
#     maps    = [m for m in maps if not m.exon_id_2 is None]
#     maps_sw = [m for m in maps if m.source=='sw_sharp' or m.source=='usearch']
#
#     for map in maps:
#
#         if map.similarity < min_similarity: continue
#         # get the raw (unaligned) sequence for the exon that maps onto human
#         exon_seqs = get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
#         if (not exon_seqs):
#             #print " exon_seqs for" , map.source
#             continue
#         [pepseq, pepseq_transl_start,
#          pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs[1:]
#
#         # rpl11 starts with an exon that translates into 2 aa's,
#         # rpl10A has a single methionine (or so they say) followed by a split codon
#         # *supposedly there is evidence at the protein level
#         # but will this give me tons of junk elsewhere? ...
#         pepseq_noX = pepseq.replace ('X','')
#         if  len(pepseq_noX)<3:
#             # if this is the first exon, and if it starts with M, we'll let it off the hook
#             # abd then if it's human, we'll also salvage it at any price
#             if first_human_exon and pepseq_noX[0] == 'M' or map.species_2=='homo_sapiens':
#                 pass
#             else:
#                 continue
#
#                 # check
#         dnaseq  = Seq (dna_seq[pepseq_transl_start:pepseq_transl_end], generic_dna)
#         if (mitochondrial):
#             pepseq2 = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
#         else:
#             pepseq2 = dnaseq.translate().tostring()
#
#
#         if (not pepseq == pepseq2):
#             continue
#
#         # inflate the compressed sequence
#         if not map.bitmap:
#             continue
#
#         bs = Bits(bytes=map.bitmap)
#         if (not bs.count(1) == len(pepseq)): continue # check bitmap has correct number of 1s
#         usi = iter(pepseq)
#         #reconst_pepseq = "".join(('-' if c=='0' else next(usi) for c in bs.bin))
#         reconst_pepseq = ''
#         for c in bs.bin:
#             if c == '0': reconst_pepseq += '-'
#             else:        reconst_pepseq += next(usi)
#
#         # come up with a unique name for this sequence
#         species       = map.species_2
#         # let's also have the start in gene here - might make our lives easier later
#         exon2 = get_exon (cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[species])
#         sequence_name = species + "_" + str(map.exon_id_2)+"_"+str(map.exon_known_2)+"_"+str(exon2.start_in_gene)
#
#
#         if reconst_pepseq:
#             sequence_pep[sequence_name] = reconst_pepseq
#             pep_aln_length = len(reconst_pepseq)
#
#             reconst_ntseq = expand_pepseq (reconst_pepseq, exon_seqs[1:], flank_length)
#             if reconst_ntseq:
#                 sequence_dna[sequence_name] = reconst_ntseq
#                 dna_aln_length = len(reconst_ntseq)
#
#     # strip common gaps
#     sequence_stripped_pep = strip_gaps (sequence_pep)
#     if not sequence_stripped_pep:
#         c=inspect.currentframe()
#         #print " in %s:%d" % ( c.f_code.co_filename, c.f_lineno)
#         return ['','']
#     # strip common gaps
#     sequence_stripped_dna = strip_gaps (sequence_dna)
#     if not sequence_stripped_dna:
#         c=inspect.currentframe()
#         #print " in %s:%d" % ( c.f_code.co_filename, c.f_lineno)
#         return ['', '']
#
#     return [sequence_stripped_pep, sequence_stripped_dna]
#
# #########################################
# def get_canonical_transl (acg, cursor, gene_id, species, strip_X = True):
#
#     canonical_translation = ""
#
#     canonical_transl_id = gene2stable_canon_transl(cursor, gene_id)
#     if ( not canonical_transl_id):
#         print("no canonical transl id found for ", gene_id)
#         return ""
#
#     cmd = acg.generate_fastacmd_protein_command (canonical_transl_id, species,
#                                                  "all", None)
#     fasta = subprocess.getoutput(cmd)
#     if (not fasta):
#         print(gene2stable (cursor, gene_id = gene_id), end=' ')
#         print("fasta not found for ", canonical_transl_id)
#         return ""
#
#     canonical_translation = ""
#     for line in fasta.split("\n"):
#         if ('>' in line):
#             continue
#         line.rstrip()
#         canonical_translation += line
#
#     if strip_X:
#         while (len(canonical_translation) and canonical_translation[0] == 'X'):
#             canonical_translation = canonical_translation[1:]
#
#     return canonical_translation
