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
        exit (1)
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
            exit(1)
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
        exit(1)
    # strip common gaps
    sequence_stripped_dna = strip_gaps (sequence_dna)
    if not sequence_stripped_dna:  
        c=inspect.currentframe()
        print " in %s:%d" % ( c.f_code.co_filename, c.f_lineno)
        exit(1)

    return [sequence_stripped_pep, sequence_stripped_dna]

