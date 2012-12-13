#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re, os

from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update

from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  Map, get_maps
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.translation        import phase2offset, translation_bounds, crop_dna, translate
from   el_utils.threads import parallelize
from bitstring import Bits

# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

#########################################
def translate_to_trivial(cursor, all_species):

    trivial_name = {}
    for species in all_species:
        taxid                 = species2taxid (cursor, species)
        trivial_name[species] = taxid2trivial(cursor, taxid)

    return trivial_name

#########################################
def fract_identity (seq1, seq2):
    
    fract_identity = 0.0
    if ( not len(seq1)):
        return fract_identity

    for i in range(len(seq1)):
        if (seq1[i] == '-'): continue
        if seq1[i] == seq2[i]: fract_identity += 1.0
    
    fract_identity /= float(len(seq1))
    return fract_identity

#########################################
def merged_sequence (template_seq, sequence_pieces, nucseq_pieces):
    
    template_length = len(template_seq)
    merged = '-'*template_length

    for piece in sequence_pieces:
        if (not len(piece) == template_length):
            print "length mismatch for aligned exons (?)"
            return merged

    # check whether any two pieces overlap
    overlap = []
    for pos in range(template_length): 

        for i in range(len(sequence_pieces)):
            if (sequence_pieces[i][pos] == '-'): continue

            for j in range (i+1, len(sequence_pieces)):
                if (sequence_pieces[j][pos] == '-'): continue
                index = str(i) + " " + str(j)
                if not index in overlap:
                    overlap.append(index)
                break

    # if yes, get rid of the less similar one
    to_delete = []
    for index in overlap:
        tmp = index.split()
        i = int(tmp[0])
        j = int(tmp[1])
        #print "overlap: ", i, j 
        if ( fract_identity (template_seq, sequence_pieces[i]) < 
             fract_identity (template_seq, sequence_pieces[j]) ):
            to_delete.append(i)
        else:
            to_delete.append(j)

    if (to_delete):
        to_delete.sort()
        to_delete.reverse()
        for i in range(len(to_delete)):
            deletable = to_delete[i]
            # not sure what is this:
            if deletable >= len(sequence_pieces): continue
            del sequence_pieces[deletable]
            del nucseq_pieces[deletable]

    # a piece of fudge - will I evr come backt to clean it up ...
    # remove flanking regions that ended up inside - this is a multiple seq alignment now
    # that neds to be fixed
    # replace the match(es) in the middle with empty strings
    for piece_ct in range(len(sequence_pieces)):
        dna_piece = nucseq_pieces  [piece_ct]
        if ( piece_ct < len(sequence_pieces) -1):
            # delete right flank"
            nucseq_pieces  [piece_ct] = dna_piece[:-flank_length]
        elif (piece_ct > 0 ):
            # delete left flank
            nucseq_pieces  [piece_ct] = dna_piece[flank_length:]

    # if not, go ahead and merge
    merged     = ""
    merged_dna = ""
    for pos in range(template_length):
        new_char  = '-'
        new_codon = '---'
        
        for piece_ct in range(len(sequence_pieces)):
            pep_piece = sequence_pieces[piece_ct]
            dna_piece = nucseq_pieces  [piece_ct]
            if pep_piece[pos] == '-': 
                continue
            else:
                new_char  =  pep_piece[pos]
                new_codon =  dna_piece[pos*3:pos*3+3]
                break
        merged     += new_char
        merged_dna += new_codon


    return [merged, merged_dna]   

#########################################
def align_nucseq_by_pepseq(aligned_pepseq, nucseq):
    if (not len(aligned_pepseq.replace('-',''))*3 == len(nucseq)):
        print aligned_pepseq.replace('-','')
        print "length mismatch: ", len(aligned_pepseq.replace('-',''))*3, len(nucseq)
        return ""
    codons = iter(map(''.join, zip(*[iter(nucseq)]*3)))
    aligned_nucseq = ''.join(('---' if c=='-' else next(codons) for c in aligned_pepseq))
    return aligned_nucseq

#########################################
def expand_pepseq (reconstructed_pep_sequence, exon_seqs):
    
    dna_aln_expanded = ""
    
    flank_length     = 10 # where should this come from?
    
    [pepseq, pepseq_transl_start, 
     pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs
    # coding dna sequence:
    cds = dna_seq[pepseq_transl_start:pepseq_transl_end]
    if not cds: 
        return ""

    aligned_nucseq  = align_nucseq_by_pepseq(reconstructed_pep_sequence, cds)
    if not aligned_nucseq: 
        print reconstructed_pep_sequence
        print cds
        print pepseq_transl_start, pepseq_transl_end," reconstruction failed"
        return ""

    effective_left_flank  = ""
    effective_right_flank = "" 
    #######
    effective_left_flank  = left_flank
    if pepseq_transl_start>0:
        effective_left_flank += dna_seq[:pepseq_transl_start]
    if len(effective_left_flank) > flank_length: 
        effective_left_flank = effective_left_flank[-flank_length:]
    effective_left_flank = effective_left_flank.lower()

    if 1:
        #######
        effective_right_flank = right_flank
        delta = len(dna_seq)-pepseq_transl_end
        if delta>0:
            effective_right_flank = dna_seq[-delta:]+effective_right_flank
        if len(effective_right_flank) > flank_length: 
            effective_right_flank = effective_right_flank[:flank_length]
        effective_right_flank = effective_right_flank.lower()

    #######
    # pad the flanking seqs to the needed length
    effective_left_flank  = effective_left_flank.rjust (flank_length, '-')
    effective_right_flank = effective_right_flank.ljust(flank_length, '-')

    dna_aln_expanded = effective_left_flank + aligned_nucseq + effective_right_flank

    return dna_aln_expanded

#########################################
def make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial):

    sequence_pep = {}
    sequence_dna = {}
    shortest_l = -1 # Uninitialized leading padding length
    shortest_r = -1 # Uninitialized trailing padding length

    # where is human sequence?

    # find all other exons that map to the human exon
    maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)

    for map in maps:
        # get the raw (unaligned) sequence for the exon that maps onto human
        exon_seqs = get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
        if (not exon_seqs):
            print map
            exit (1)
        [pepseq, pepseq_transl_start, 
         pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs[1:]

        # check
        dnaseq  = Seq (dna_seq[pepseq_transl_start:pepseq_transl_end], generic_dna)
        if (mitochondrial):
            pepseq2 = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
        else:
            pepseq2 = dnaseq.translate().tostring()
        

        if (not pepseq == pepseq2):
            print " ! ", pepseq
            print " ! ", pepseq2
            exit (1)
            
        # inflate the compressed sequence
        bs = Bits(bytes=map.bitmap)
        if (not bs.count(1) == len(pepseq)): continue # check bitmap has correct number of 1s
        usi = iter(pepseq)
        reconst_pepseq = "".join(('-' if c=='0' else next(usi) for c in bs.bin))
  
        # come up with a unique name for this sequence
        species       = map.species_2
        sequence_name = species + "_" + str(map.exon_id_2)+"_"+str(map.exon_known_2)

        if reconst_pepseq: 
            sequence_pep[sequence_name] = reconst_pepseq
            pep_aln_length = len(reconst_pepseq)

            reconst_ntseq = expand_pepseq (reconst_pepseq, exon_seqs[1:])
            if reconst_ntseq: 
                sequence_dna[sequence_name] = reconst_ntseq
                dna_aln_length = len(reconst_ntseq)


    # strip common gaps
    all_gaps = {}  
    for pos in range(pep_aln_length):
        all_gaps[pos] = True
        for name, seq in sequence_pep.iteritems():
            if (not seq[pos]=='-'):
                all_gaps[pos] = False
                break

    sequence_stripped_pep = {}
    for name, seq in sequence_pep.iteritems():
        sequence_stripped_pep[name] = ""
        for pos in range(pep_aln_length):
            if all_gaps[pos]: continue
            sequence_stripped_pep[name] += seq[pos]

    # strip common gaps
    all_gaps = {}  
    for pos in range(dna_aln_length):
        all_gaps[pos] = True
        for name, seq in sequence_dna.iteritems():
            if (not seq[pos]=='-'):
                all_gaps[pos] = False
                break

    sequence_stripped_dna = {}
    for name, seq in sequence_dna.iteritems():
        sequence_stripped_dna[name] = ""
        for pos in range(dna_aln_length):
            if all_gaps[pos]: continue
            sequence_stripped_dna[name] += seq[pos]

    return [sequence_stripped_pep, sequence_stripped_dna]



#########################################
def print_notes (notes_fnm, orthologues, exons, sorted_species, specid2name, human_stable_id, source):

    # write to string
    out_string  = "% Notes to accompany the alignment of (tentative) orthologues\n"
    out_string += "%% for the canonical transcript of the human gene %s," %  human_stable_id
    out_string += "\n" 
    out_string += "% The alignment shows the exons that correspond to human transcript,\n" 
    out_string += "% from the following genes: \n" 
    out_string += "%% %7s  %25s  %30s \n" % ('name_short', 'name_long', 'stable_id')

    for species in sorted_species:
        [orth_id, spec_id] =  used_orthologue[species]
        spec_long = specid2name[spec_id]
        stable_id = get_stable_id (orth_id, spec_long)
        out_string += "%7s %30s  %30s \n" % (species, specid2name[spec_id], stable_id)

    out_string += "\n" 
    out_string += "% The following exons were used in the alignment\n" 

    for i in range(len(used_exons['HOM_SAP'])):
        out_string += "%% exon %3d\n" % i
        out_string += "%% %7s  %15s %15s  %40s\n" % ('name_short', 'gene_from', 'gene_to', 'source')
        for species in sorted_species:
            if ( used_exons.has_key(species) and i < len(used_exons[species])):
                ret =  used_exons[species][i]
                if ret:
                    [exon_id, is_known,  is_canonical, start_in_gene, end_in_gene, protein_seq, analysis_id] = ret
                
                    out_string += "%7s  %15d %15d  %40s \n" % (species, start_in_gene,
                                                           end_in_gene, source[species][analysis_id])
                else:
                    out_string += "%7s  %15s %15s  %40s \n" % (species, "-",  "-",  "-")
            else:
                out_string += "%7s  %15s %15s  %40s \n" % (species, "-",  "-",  "-")
                
                

    of = open (notes_fnm, "w")
    print >> of, out_string
    of.close()

    return True

#########################################
def parse_aln_name (name):
    fields     = name.split("_")
    exon_id    = int(fields[-2])
    exon_known = int(fields[-1])
    species    =  "_".join(fields[:-2])
    return [species, exon_id, exon_known]

#########################################
def get_name (seq_name, species, ortho_gene_id):
    sequence_name = ""

    if ( not seq_name.has_key(species)): # we see this species for the first time
        sequence_name = species
        seq_name[species] = {ortho_gene_id:sequence_name}
    else:
        # we have seen an exon from this species and it came from the same gene
        if seq_name[species].has_key(ortho_gene_id):
            sequence_name = seq_name[species][ortho_gene_id]
        else:
            # if sequence for the species exists,
            # but not for this gene, mangle the name
            number_of_genes = len(seq_name[species].keys()) + 1
            sequence_name = species + "_" + str(number_of_genes)
            seq_name[species][ortho_gene_id] = sequence_name
    return sequence_name

#########################################

def make_alignments ( gene_list, db_info):

    [local_db, ensembl_db_name] = db_info

    verbose      = False
    flank_length = 10

    if local_db:
        db     = connect_to_mysql()
        cfg    = ConfigurationReader()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    # find db ids adn common names for each species db
    [all_species, ensembl_db_name] = get_species (cursor)
    # walk the taxonomical tree, and sort the species according to
    # the (distance of) the last common ancestor
    sorted_species = species_sort(cursor, all_species, 'homo_sapiens')
    # in the afa headers use 'trivial' names for the species: cow, dog, pig, ...
    trivial_name   = translate_to_trivial(cursor, all_species)

    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
    #gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    # for each human gene
    gene_ct = 0
    #for gene_id in gene_list:
    for gene_id in [412667]: #  wls       
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)
        afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        #if (os.path.exists(afa_fnm) and os.path.getsize(afa_fnm) > 0):
        #    continue
        gene_ct += 1
        if (not gene_ct%100): print gene_ct, "out of ", len(gene_list)
        if verbose: print gene_id, stable_id, get_description (cursor, gene_id)

        # find all exons we are tracking in the database
        human_exons = gene2exon_list(cursor, gene_id)
        canonical_exons = []
        for human_exon in human_exons:
            if not human_exon.is_canonical:
                continue
            canonical_exons.append(human_exon)

        # the exons are not guaranteed to be in order
        canonical_exons.sort(key=lambda exon: exon.start_in_gene)

        # reconstruct the alignment with orthologues
        mitochondrial        = is_mitochondrial(cursor, gene_id)
 
        alnmt_pep    = {}
        alnmt_dna    = {}
        has_a_map    = False
        for human_exon in canonical_exons:
            [alnmt_pep[human_exon], alnmt_dna[human_exon]]  = \
                make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial)   
            if alnmt_pep[human_exon]: has_a_map=True

        # >>>>>>>>>>>>>>>>>>
        if not has_a_map: continue

        # find which species we have in this story
        dna_sequence = {}
        pep_sequence = {}
        seq_name     = {}
        for human_exon in canonical_exons:
            for name,seq in alnmt_pep[human_exon].iteritems():
                (species, exon_id, exon_known) = parse_aln_name(name)
                ortho_gene_id = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
                # gene --> name -- retrieve old name, or construct new one
                sequence_name = get_name (seq_name, trivial_name[species], ortho_gene_id) 
                if not pep_sequence.has_key(sequence_name): 
                    pep_sequence[sequence_name] = {}
                if not pep_sequence[sequence_name].has_key(human_exon): 
                    pep_sequence[sequence_name][human_exon] = []
                pep_sequence[sequence_name][human_exon].append(seq)

        for sequence_name in pep_sequence.keys():
            print sequence_name
            for human_exon in canonical_exons:
                if pep_sequence[sequence_name].has_key(human_exon):
                    if ( len(pep_sequence[sequence_name][human_exon]) == 1):
                        print human_exon.start_in_gene, pep_sequence[sequence_name][human_exon][0]
                    else:
                        # if two sequences map to the same human exon, merge
                        template = alnmt_pep[human_exon]['homo_sapiens']
                        [pep, dna] = merged_sequence (template, alnmt_pep[sequence_name][human_exon], 
                                                       alnmt_dna[sequence_name][human_exon])
                        print human_exon.start_in_gene, pep_sequence[sequence_name][human_exon]
        exit(1)
        # stich the exons together, this time taking into account that the alignment
        # doesn't have to be one to one
                       
        headers        = []
        output_pep     = {}
        output_dna     = {}
        ortholog_count = {}
        exit(1)
        # find place for best_afa -- for the moment can put it to the scratch space:
        afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        output_fasta (afa_fnm, headers, output_pep)
        print afa_fnm
        afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        output_fasta (afa_fnm, headers, output_dna)
        print afa_fnm

        exit(1)

        # notes to accompany the alignment:
        notes_fnm  = "{0}/notes/{1}.txt".format(cfg.dir_path['afs_dumps'], stable_id)
        print notes_fnm
        print_notes (notes_fnm, orthologues, exons, sorted_species, specid2name, human_stable_id, source)
        
        

#########################################
def main():
    
    no_threads = 1

    local_db = False

    if local_db:
        db = connect_to_mysql()
    else:
        db = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)


    species                        = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_list                      = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
    cursor.close()
    db.close()

    parallelize (no_threads, make_alignments, gene_list, [local_db, ensembl_db_name])
    
    return True


#########################################
if __name__ == '__main__':
    main()

'''
    #for gene_id in [412667]: #  wls
    #for gene_id in [378768]: #  p53
    #for gene_id in [378766]: #  dynein
'''
