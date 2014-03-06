import MySQLdb, commands, re, os
import inspect
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna
from bitstring import Bits

from  ensembl import  *
from  utils   import  *
from  map     import  *



#########################################
def strip_stop(pepseq):
    if (not pepseq or len(pepseq)==0):
        return pepseq
    if ( pepseq[-1] == '*'):
        pepseq = pepseq[:-1]
    return pepseq

#########################################
def  transl_reconstruct (cursor,  gene_id, gene_seq, canonical_coding_exons, 
                         is_mitochondrial, verbose = False):

    """
    Given tha dna sequence, gene information an the list of its canonical exons, reconstruct canonical translation.

    Pay attention to whether the gene is mitochondrial.
    Here in particular we can catch the 'false stop codons' corresponding to selenocysteines.
    (Because they are stored by their position
    in the translation.)

    """

    canonical_exon_pepseq = {}
    translated_seq = "" 

    [can_transl_start_exon, can_transl_start_position,
     can_transl_end_exon, can_transl_end_position] = canonical_transl_info (cursor, gene_id)


    # do we have any selenocysteines by any chance
    selenoC_pos = get_selenocysteines (cursor,  gene_id)

    carry = ""
    ok_so_far = True
    # sanity checking
    for exon in canonical_coding_exons:
        #print
        #print "exon", exon.exon_id
        #find exon sequence within the gene
        start = exon.start_in_gene
        if (exon is canonical_coding_exons[0]):
            if ( not exon.exon_id == can_transl_start_exon ):
                print " error start cantransl:  gene_id ",  gene_id,
                print  " exon_id ", exon.exon_id, " canon: ", can_transl_start_exon
                return [{}, ""]
                #exit (1)
            start +=  exon.canon_transl_start

        if ( exon is canonical_coding_exons[-1]):
            if ( not exon.exon_id == can_transl_end_exon ):
                print " error end cantransl:  gene_id ",  gene_id,
                print  " exon_id ", exon.exon_id, " canon: ", can_transl_end_exon
                return [{}, ""]
                #exit (1)
            end = exon.start_in_gene + exon.canon_transl_end
        else:
            end = exon.end_in_gene

        if (not exon.phase == -1 and not exon.phase == len(carry)):
            #print "Houston we have a problem: exon phase =", exon.phase,
            #print " the length of carry =", len(carry), 
            #print " (gene_id %d, exon_id %d) " % (gene_id, exon.exon_id)
            if (  exon.phase ):
                start += 3-exon.phase
            carry = ""


        exon_seq     =  gene_seq[ start: end+1]
        exon_seq_for_transl_purposes = carry + exon_seq

        remainder    = len(exon_seq_for_transl_purposes)%3
        if ( remainder == 0 ):
            carry = ""
        elif (remainder == 1 ):
            carry    = exon_seq_for_transl_purposes[-1:]
            exon_seq_for_transl_purposes  = exon_seq_for_transl_purposes[:-1]
        else:
            carry    = exon_seq_for_transl_purposes[-2:]
            exon_seq_for_transl_purposes = exon_seq_for_transl_purposes[:-2]

        dnaseq = Seq (exon_seq_for_transl_purposes, generic_dna)
        pepseq = dnaseq.translate()

        # turn to the corresponding BioPython object
        if ( is_mitochondrial ):
            pepseq = dnaseq.translate(table="Vertebrate Mitochondrial")
        else:
            pepseq = dnaseq.translate()
            
        # replace stop codons from selenoC positions, if there are such
        if (selenoC_pos):
            for pos_in_full_length_translation in selenoC_pos:
                pos_in_exon_translation = pos_in_full_length_translation-len(translated_seq)
                if pos_in_exon_translation<0 or pos_in_exon_translation>len(pepseq):
                    continue
                tempseq = pepseq.tomutable()
                tempseq[pos_in_exon_translation] = 'U'
                pepseq  = tempseq.toseq()
 
        # strip the last stop codon only
        if ( exon is canonical_coding_exons[-1]):
            pepseq = strip_stop(pepseq)  

        pepseq0 =  pepseq
        if verbose:
            print "phase 0", pepseq0

        # if there are stil stop codons we'll give another shot 
        # to the possibility that it is mitochondrial, (and we ddin't know it)
        # after that we cry foul
        if ( '*' in pepseq):
            pepseq = dnaseq.translate(table="Vertebrate Mitochondrial")

        # strip the last stop codon only
        if ( exon is canonical_coding_exons[-1]):
            pepseq = strip_stop(pepseq)  

        # some further  desperate measures 
        ok_so_far = True
        if ( '*' in pepseq):
            ok_so_far = False
            dnaseq = Seq (exon_seq_for_transl_purposes[1:], generic_dna)
            pepseq = dnaseq.translate()
            pepseq = strip_stop(pepseq)
            pepseq1 =  pepseq
            if verbose:
                print "phase 1", pepseq1
        else:
            ok_so_far = True

        if (not ok_so_far and '*' in pepseq):
            dnaseq = Seq (exon_seq_for_transl_purposes[2:], generic_dna)
            pepseq = dnaseq.translate()
            pepseq = strip_stop(pepseq)
            pepseq2 =  pepseq
            if verbose:
                print "phase 2", pepseq2
        else:
            ok_so_far = True

        if (not ok_so_far and  '*' in pepseq):
            if verbose:
                print "Error: stop codon "
        else:
            ok_so_far = True

        if ( not ok_so_far):
            return [{}, ""]
        translated_seq += pepseq # I need the seq in selenoC - to decide
                                 # where the position of  U should be
        canonical_exon_pepseq[exon.exon_id] = pepseq.tostring()

    return [canonical_exon_pepseq,translated_seq] 

#########################################
def compare_seqs (canonical_translation, translated_seq, verbose=False):

    """
    Return true if the two input sequences differ in no more than two positions.
    """
    comparison_ok = True

    while (len(translated_seq) and translated_seq[0] == 'X'):
        translated_seq = translated_seq[1:]

    difference = len(translated_seq) - len(canonical_translation)
    if ( abs(difference) > 3):
        comparison_ok = False
        if verbose:
            print
            print ">canon"
            print canonical_translation
            print ">exons"
            print translated_seq
            print
    else:
        diff  =  0
        start = -1
        for i in range(len(translated_seq)):
            if ( i >= len(canonical_translation)):
                break
            if (not translated_seq[i] ==  canonical_translation[i]):
                diff += 1
                if start < 0:
                    start = i
        if (diff > 2):
            comparison_ok = False
            if verbose:
                print
                print ">canon"
                print canonical_translation
                print ">exons"
                print translated_seq
                print translated_seq[start], canonical_translation[start]
                print "nuber of  diff sites: ", diff, " starting from ", start
                print

    return comparison_ok



#########################################
def  get_gene_seq (acg, cursor, gene_id, species, verbose = False):

    """
    Given gene_id, return dna region which reproduces the correct canonical translation.
    """
    null = ["",{}, "", "", None, None]

    #########################################
    # which file should we be looking in, which sequence, from where to where
    ret = get_primary_seq_info (cursor, gene_id, species)
    if (not ret):
        return null
    [seq_name, file_names, seq_region_start, seq_region_end, 
     seq_region_strand, is_mitochondrial] = ret
    # i'm not quite clear why Ensembl is doing this, but sometimes we need the alternative
    # region - ("PATCH" deposited as tte right sequence, but its missing most of the gene)
    # so first establish whether it is the case: find canonical translation.
    canonical_translation  = get_canonical_transl (acg, cursor, gene_id, species)
    if not canonical_translation: return null
    # find all canonical exons associated with the gene id
    canonical_coding_exons = get_canonical_exons (cursor, gene_id)
    # extract raw gene  region TODO - store the information about which 
    # file_name we ended up using, and the start and end in that region
    [gene_seq, file_name] = extract_gene_seq (acg, species, seq_name, file_names, seq_region_strand,  
                                             seq_region_start, seq_region_end)
    # reconstruct the translation from the raw gene_seq and exon boundaries
    [canonical_exon_pepseq,translated_seq] = transl_reconstruct (cursor, gene_id, gene_seq, 
                                                                 canonical_coding_exons, is_mitochondrial)
    if (translated_seq):
        # compare the two sequences and cry foul if they are not the same:
        comparison_ok = compare_seqs (canonical_translation, translated_seq)
    else:
        comparison_ok = False
    # if we succefully translated the exons, and came up with the same answer 
    # as the canonical translation, we are done here
    if (comparison_ok):
        return [gene_seq, canonical_exon_pepseq, file_name, seq_name, seq_region_start, seq_region_end]
    if verbose:
        print "Using primary seq info: failed comparison with canonical sequence."
        print "canonical:"
        print canonical_exon_pepseq
        print "translated:"
        print translated_seq
 
    #########################################
    # otherwise repeat the procedure with the alternative seq info:
    ret = get_alt_seq_info (cursor, gene_id, species)
    if (not ret):
        return null
    [seq_name, file_names, seq_region_start, seq_region_end, 
     seq_region_strand, is_mitochondrial] = ret
    # check  canonical translation
    canonical_translation  = get_canonical_transl (acg, cursor, gene_id, species)
    # find all canonical exons associated with the gene id
    canonical_coding_exons = get_canonical_exons (cursor, gene_id)
    # extract raw gene  region
    [gene_seq, file_name] = extract_gene_seq(acg, species, seq_name, file_names, seq_region_strand,  
                                seq_region_start, seq_region_end)
    # reconstruct the translation from the raw gene_seq and exon boundaries
    [canonical_exon_pepseq,translated_seq] = transl_reconstruct (cursor, gene_id, gene_seq, canonical_coding_exons, 
                                                                  is_mitochondrial)
    if (translated_seq):
        # compare the two sequences and cry foul if they are not the same:
        comparison_ok = compare_seqs (canonical_translation, translated_seq)
    else:
        comparison_ok = False
    # if we succefully translated the exons, and came up with the same answer 
    # as the canonical translation, we are done here
    if (comparison_ok):
        return [gene_seq, canonical_exon_pepseq, file_name, seq_name, seq_region_start, seq_region_end]

    if verbose:
        print "Using alt seq info: failed comparison with canonical sequence."

    return null 


#########################################
def get_reliable_orthos(cursor, ensembl_db_name, gene_id):

    all_orthologues = []

    # one2one   orthologues
    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    known_orthologues      = get_orthos (cursor, gene_id, 'orthologue')
    # not-clear orthologues
    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    unresolved_orthologues = get_orthos (cursor, gene_id, 'unresolved_ortho')

    # get rid of the unresolved orthologues if a resolved orthologue for the species already exists
    species_with_known_orthologues = []
    for  [ortho_gene_id, ortho_species] in known_orthologues:
        species_with_known_orthologues.append(ortho_species)
    all_orthologues = known_orthologues
    for  [ortho_gene_id, ortho_species] in unresolved_orthologues:
        if ortho_species in species_with_known_orthologues: continue
        all_orthologues.append( [ortho_gene_id, ortho_species] )

    return all_orthologues



#########################################
def sort_names (sorted_species, alignment):

    sorted_names = []
    for species in sorted_species:
        for seq_name  in alignment.keys():

            if seq_name[-1].isdigit():
                aux = seq_name.split("_")
                base_name = "_".join(aux[:-1])

                if base_name [-1].isdigit():
                    aux = base_name.split("_")
                    base_name = "_".join(aux[:-1])
            else:
                base_name = seq_name
            if (species == base_name):
                sorted_names.append(seq_name)
    return sorted_names

#########################################
def align_nucseq_by_pepseq (aligned_pepseq, nucseq):

    # coding dna sequence, by assumption:
    cds = nucseq
    translated_cds = Seq(cds).translate().tostring()
    if not len(aligned_pepseq.replace('-','')) == len(translated_cds):
        print aligned_pepseq.replace('-','')
        print "in align_nucseq_by_pepseq():  length mismatch: ", 
        print len(aligned_pepseq.replace('-','')), len(translated_cds)
        return ""
    codon = iter(map(''.join, zip(*[iter(nucseq)]*3)))
    #aligned_nucseq = ''.join(('---' if c=='-' else next(codon) for c in aligned_pepseq))
    aligned_nucseq = ''
    for c in aligned_pepseq:
        if c == '-': aligned_nucseq += '---'
        else:        aligned_nucseq += next(codon)
    return aligned_nucseq

#########################################
def expand_pepseq (aligned_pep_sequence, exon_seqs, flank_length):
    
    dna_aln_expanded = ""

    # check if this is a padding seqeunce:
    if not exon_seqs or  not aligned_pep_sequence.replace('-',''):
        dna_aln_expanded = '-'*(2*flank_length+3*len(aligned_pep_sequence) )
        return dna_aln_expanded

    [pepseq, pepseq_transl_start, 
     pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs
    if not pepseq or pepseq_transl_start is None or pepseq_transl_end is None:
        return ""
    # coding dna sequence:
    cds = dna_seq[pepseq_transl_start:pepseq_transl_end]
    if not cds: 
        return ""

    aligned_nucseq  = align_nucseq_by_pepseq(aligned_pep_sequence, cds)
    if not aligned_nucseq: 
        print aligned_pep_sequence
        print Seq(cds).translate().tostring()
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
def check_has_sw_exons (cursor, ensembl_db_name, human_exon_id, human_exon_known, minsim):

    has_sw_exons = False

    # find all other exons that map to the human exon
    maps    = get_maps(cursor, ensembl_db_name, human_exon_id, human_exon_known)
    maps    = filter (lambda m: not m.exon_id_2 is None, maps)
    maps_sw = filter (lambda m: m.source=='sw_sharp' and m.similarity >minsim, maps)

    if maps_sw:  
        has_sw_exons = True
 
    return has_sw_exons

#########################################
def make_exon_alignment(cursor, ensembl_db_name, human_exon_id, human_exon_known, mitochondrial, 
                        min_similarity,  flank_length):

    sequence_pep = {}
    sequence_dna = {}
    shortest_l = -1 # Uninitialized  leading padding length
    shortest_r = -1 # Uninitialized trailing padding length

    pep_aln_length = 0
    dna_aln_length = 0
    # find all other exons that map to the human exon
    maps    = get_maps(cursor, ensembl_db_name, human_exon_id, human_exon_known)
    maps    = filter (lambda m: not m.exon_id_2 is None, maps)
    maps_sw = filter (lambda m: m.source=='sw_sharp' or m.source=='usearch', maps)

    for map in maps:

        if map.similarity < min_similarity: continue
        # get the raw (unaligned) sequence for the exon that maps onto human
        exon_seqs = get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
        if (not exon_seqs):
            print " exon_seqs for" , map.source
            continue
        [pepseq, pepseq_transl_start, 
         pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs[1:]

        if     len(pepseq)<3: continue
        pepseq_noX = pepseq.replace ('X','')
        if len(pepseq_noX)<3: continue
       

        # check
        dnaseq  = Seq (dna_seq[pepseq_transl_start:pepseq_transl_end], generic_dna)
        if (mitochondrial):
            pepseq2 = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
        else:
            pepseq2 = dnaseq.translate().tostring()
        

        if (not pepseq == pepseq2):
            continue
            
        # inflate the compressed sequence
        if not map.bitmap:
            continue

        bs = Bits(bytes=map.bitmap)
        if (not bs.count(1) == len(pepseq)): continue # check bitmap has correct number of 1s
        usi = iter(pepseq)
        #reconst_pepseq = "".join(('-' if c=='0' else next(usi) for c in bs.bin))
        reconst_pepseq = ''
        for c in bs.bin:
            if c == '0': reconst_pepseq += '-'
            else:        reconst_pepseq += next(usi)

        # come up with a unique name for this sequence
        species       = map.species_2
        sequence_name = species + "_" + str(map.exon_id_2)+"_"+str(map.exon_known_2)


        if reconst_pepseq: 
            sequence_pep[sequence_name] = reconst_pepseq
            pep_aln_length = len(reconst_pepseq)

            reconst_ntseq = expand_pepseq (reconst_pepseq, exon_seqs[1:], flank_length)
            if reconst_ntseq: 
                sequence_dna[sequence_name] = reconst_ntseq
                dna_aln_length = len(reconst_ntseq)

    # strip common gaps
    sequence_stripped_pep = strip_gaps (sequence_pep)
    if not sequence_stripped_pep:  
        c=inspect.currentframe()
        print " in %s:%d" % ( c.f_code.co_filename, c.f_lineno)
        return ['','']
    # strip common gaps
    sequence_stripped_dna = strip_gaps (sequence_dna)
    if not sequence_stripped_dna:  
        c=inspect.currentframe()
        print " in %s:%d" % ( c.f_code.co_filename, c.f_lineno)
        return ['', '']

    return [sequence_stripped_pep, sequence_stripped_dna]

#########################################
def get_canonical_transl (acg, cursor, gene_id, species, strip_X = True):

    canonical_translation = ""

    canonical_transl_id = gene2stable_canon_transl(cursor, gene_id)
    if ( not canonical_transl_id):
        print "no canonical transl id found for ", gene_id
        return ""

    cmd = acg.generate_fastacmd_protein_command (canonical_transl_id, species, 
                                                 "all", None)
    fasta = commands.getoutput(cmd)
    if (not fasta):
        print gene2stable (cursor, gene_id = gene_id), 
        print "fasta not found for ", canonical_transl_id
        return ""

    canonical_translation = ""
    for line in fasta.split("\n"):
        if ('>' in line):
            continue
        line.rstrip()
        canonical_translation += line

    if strip_X:
        while (len(canonical_translation) and canonical_translation[0] == 'X'):
            canonical_translation = canonical_translation[1:]

    return canonical_translation
