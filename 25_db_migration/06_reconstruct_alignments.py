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
def merged_sequence (template_seq, sequence_pieces, nucseq_pieces, flank_length=10):
    
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

    # a piece of fudge - will I ever come backt to clean it up ...
    # remove flanking regions that ended up inside - this is a multiple seq alignment now
    # this needs to be fixed
    # replace the match(es) in the middle with empty strings
    left_flank  = nucseq_pieces [0][0:flank_length]
    right_flank = nucseq_pieces[-1][-flank_length:]
    for piece_ct in range(len(sequence_pieces)):
        dna_piece = nucseq_pieces  [piece_ct]
        nucseq_pieces  [piece_ct] = dna_piece[flank_length:-flank_length]

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

    merged_dna = left_flank+merged_dna+right_flank
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
def strip_gaps (sequence):
    
    seq_stripped = {}

    all_gaps = {}  
    
    for name, seq in sequence.iteritems():
        aln_length = len(seq)
        break

    for pos in range(aln_length):
        all_gaps[pos] = True
        for name, seq in sequence.iteritems():
            if (not seq[pos]=='-'):
                all_gaps[pos] = False
                break

    for name, seq in sequence.iteritems():
        seq_stripped[name] = ""
        for pos in range(aln_length):
            if all_gaps[pos]: continue
            seq_stripped[name] += seq[pos]

    return seq_stripped

#########################################
def make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial):

    sequence_pep = {}
    sequence_dna = {}
    shortest_l = -1 # Uninitialized leading padding length
    shortest_r = -1 # Uninitialized trailing padding length

    pep_aln_length = 0
    dna_aln_length = 0
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

        if len(pepseq)<2: continue

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
    sequence_stripped_pep = strip_gaps (sequence_pep)
    # strip common gaps
    sequence_stripped_dna = strip_gaps (sequence_dna)

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
def  find_human_template(exon_alignment):

    template_name = filter (lambda spec: 'homo_sapiens' in spec, exon_alignment.keys())[0]
    template_seq  = exon_alignment[template_name]

    return [template_name, template_seq]

#########################################
def sort_names (sorted_species, alignment):

    sorted_names = []
    for species in sorted_species:
        for seq_name  in alignment.keys():
            if (species in seq_name):
                sorted_names.append(seq_name)
    return sorted_names

#########################################
def gapped_exon (exon_alignment_pep, exon_almt_dna):
    [pep, dna] = ["",""]
    [template_name, template_pep] = find_human_template(exon_alignment_pep)
    pep = '-'*len(template_pep)
    template_dna = exon_almt_dna[template_name]
    dna = '-'*len(template_dna)
    return [pep, dna]
                     
#########################################
def flags_init( flags, seq_name, flag_name):
    if not flags.has_key(seq_name): 
        flags[seq_name] = {}
    if not flags[seq_name].has_key(flag_name):
        flags[seq_name][flag_name] = []
    return

#########################################
def fix_one2many (cfg, acg, exon_seq_name, concat_name, list_of_human_exons, seqid, alnmt_pep, output_pep ):
    
    new_alignment = {}

    exon_numbers = map (lambda ex: seqid[ex], list_of_human_exons)
    smallest_id  = exon_numbers[0]
    largest_id   = exon_numbers[-1]
    # sanity
    prev = ""
    for human_exon in list_of_human_exons:
        current = alnmt_pep[human_exon][exon_seq_name].replace ("-", "")
        if prev and not prev == current: # should be all one and the same
            print "oink? "
            exit (1)
        prev = current
    seq_to_fix = current
    # pull  the slice out of the alignment
    slice_seq = {}
    # exon boundaries:
    delimiter = re.compile("Z")
    exon_aln_start = []
    exon_aln_end   = []
    start          =  0
    prev_end       =  0
    # use human as the reference - in other species the boundaries might
    # be at different positions
    for match in delimiter.finditer(output_pep['human']):
        
        start = prev_end 
        end   = match.start()
        exon_aln_start.append(start)
        exon_aln_end.append(end)
        prev_end   =  match.end()

    start = prev_end + 1
    end   = len(output_pep['human'])
    exon_aln_start.append(start)
    exon_aln_end.append(end)
    ####################################
    slice_start = exon_aln_start[smallest_id]
    slice_end   = exon_aln_end  [largest_id]
    for name, seq in output_pep.iteritems():
        if (name== concat_name): continue
        slice_seq[name] = seq[slice_start:slice_end]
        
        
    tmp_name = exon_seq_name+"_"+"_".join(map (lambda x: str(x), exon_numbers))
    afa_fnm  = "{0}/{1}.afa".format(cfg.dir_path['scratch'], tmp_name)
    output_fasta (afa_fnm, slice_seq.keys(), slice_seq)

    tmp_name  = exon_seq_name
    fasta_fnm = "{0}/{1}.fa".format(cfg.dir_path['scratch'], tmp_name)
    output_fasta (fasta_fnm, [concat_name], {concat_name:seq_to_fix})

    tmp_name = exon_seq_name+"_"+"_".join(map (lambda x: str(x), exon_numbers))
    out_fnm  = "{0}/{1}.out.afa".format(cfg.dir_path['scratch'], tmp_name)

    mafftcmd = acg.generate_mafft_profile (afa_fnm, fasta_fnm, out_fnm)
    ret      = commands.getoutput(mafftcmd)
    ret      = commands.getoutput("cat "+out_fnm)

    realigned = {}
    seq = ""
    for line in ret.split('\n'):
        if '>' in line:
            if ( seq ):
                realigned[name] = seq
            name = line.replace (">", "")
            name = name.replace (" ", "")
            seq = ""
        else:
            seq += line
    if seq: realigned[name] = seq

    # replace the slice with the re-aligned one
    new_alignment = {}
    for name, seq in output_pep.iteritems():
        new_alignment[name]  = ""
        new_alignment[name] += seq[:slice_start] # this should presumable include "Z"
        new_alignment[name] += realigned[name]
        new_alignment[name] += seq[slice_end:]
 
    tmp_name  = exon_seq_name
    fasta_fnm = "{0}/{1}.ctrl.afa".format(cfg.dir_path['scratch'], tmp_name)
    output_fasta (fasta_fnm, new_alignment.keys(), new_alignment)
    
    return new_alignment


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
    #for gene_id in [412667]: #  wls   
    for gene_id in [378768]: #  p53

        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)
        afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        #if (os.path.exists(afa_fnm) and os.path.getsize(afa_fnm) > 0):
        #    continue
        gene_ct += 1
        if (not gene_ct%100): print gene_ct, "out of ", len(gene_list)
        if verbose: print gene_id, stable_id, get_description (cursor, gene_id)

        # find all exons we are tracking in the database
        human_exons     = gene2exon_list(cursor, gene_id)
        canonical_exons = []
        for human_exon in human_exons:
            if not human_exon.is_canonical or  not human_exon.is_coding:
                continue
            canonical_exons.append(human_exon)

        # the exons are not guaranteed to be in order
        canonical_exons.sort(key=lambda exon: exon.start_in_gene)

        # >>>>>>>>>>>>>>>>>>
        # reconstruct the per-exon alignment with orthologues
        mitochondrial = is_mitochondrial(cursor, gene_id)
 
        alnmt_pep = {}
        alnmt_dna = {}
        has_a_map = False
        ct        = 0
        seqid     = {}
        for human_exon in canonical_exons:
            [alnmt_pep[human_exon], alnmt_dna[human_exon]]  = \
                make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial)   
            if alnmt_pep[human_exon]: has_a_map=True
            seqid[human_exon] = ct
            ct += 1

        # do we have a sequence mapping to multiple human exons?
        ortho_exon_to_human_exon = {}
        for human_exon in canonical_exons:
            for exon_seq_name in alnmt_pep[human_exon].keys():
                if not ortho_exon_to_human_exon.has_key(exon_seq_name):
                    ortho_exon_to_human_exon[exon_seq_name] = [human_exon]
                else:
                    ortho_exon_to_human_exon[exon_seq_name].append(human_exon)

        # >>>>>>>>>>>>>>>>>>
        if not has_a_map: continue

        # >>>>>>>>>>>>>>>>>>
        # find which species we have, and for how many exons
        # we may have two orthologues for the same species
        dna_sequence = {}
        pep_sequence = {}
        seq_name     = {}
        flags        = {}
        parent_seq_name = {}
        for human_exon in canonical_exons:
            for exon_seq_name, exon_seq in alnmt_pep[human_exon].iteritems():
                (species, exon_id, exon_known) = parse_aln_name(exon_seq_name)
                ortho_gene_id = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
                # gene --> name -- retrieve old name, or construct new one
                parent_seq_name[exon_seq_name] = get_name (seq_name, trivial_name[species], ortho_gene_id) 
                

        # >>>>>>>>>>>>>>>>>>
        # flag the cases when one ortho exon maps to many human for later
        for exon_seq_name, concat_seq_name  in parent_seq_name.iteritems():
            if ( len(ortho_exon_to_human_exon[exon_seq_name]) > 1):
                flags_init(flags, concat_seq_name, 'one_ortho2many_human')
                to_fix = [exon_seq_name, ortho_exon_to_human_exon[exon_seq_name]]
                flags[concat_seq_name]['one_ortho2many_human'].append(to_fix)

        # >>>>>>>>>>>>>>>>>>
        for human_exon in canonical_exons:
            for exon_seq_name, exon_seq in alnmt_pep[human_exon].iteritems():
                concat_seq_name = parent_seq_name[exon_seq_name]
                if not pep_sequence.has_key(concat_seq_name): 
                    pep_sequence[concat_seq_name] = {}
                    dna_sequence[concat_seq_name] = {}
                if not pep_sequence[concat_seq_name].has_key(human_exon): 
                    pep_sequence[concat_seq_name][human_exon] = []
                    dna_sequence[concat_seq_name][human_exon] = []
                pep_sequence[concat_seq_name][human_exon].append(exon_seq)
                dna_sequence[concat_seq_name][human_exon].append(alnmt_dna[human_exon][exon_seq_name])
  
        # >>>>>>>>>>>>>>>>>>
        # concatenate the aligned exons for each species, taking into account that the alignment
        # doesn't have to be one to one
        headers     = []
        output_pep  = {}
        output_dna  = {}
        for concat_seq_name in pep_sequence.keys():
            output_pep[concat_seq_name] = ""
            output_dna[concat_seq_name] = ""
            flagged_exons = []
            if concat_seq_name in flags.keys() and flags[concat_seq_name].has_key('one_ortho2many_human'):
                for to_fix in flags[concat_seq_name]['one_ortho2many_human']:
                    flagged_exons = flagged_exons+to_fix[1]
            for human_exon in canonical_exons:
                if human_exon in flagged_exons:
                    # one ortho seq maps to multiple human exons
                    [pep, dna] = gapped_exon (alnmt_pep[human_exon], alnmt_dna[human_exon])
                elif pep_sequence[concat_seq_name].has_key(human_exon):
                    if ( len(pep_sequence[concat_seq_name][human_exon]) == 1):
                        # we have a neat one-to-one mapping
                        pep = pep_sequence[concat_seq_name][human_exon][0]
                        dna = dna_sequence[concat_seq_name][human_exon][0]
                    else:
                        # if two sequences map to the same human exon, merge
                        [pep, dna] = gapped_exon (alnmt_pep[human_exon], alnmt_dna[human_exon])
                        # flag for later
                        flags_init(flags, concat_seq_name, 'many_ortho2one_human')
                        to_fix = [human_exon, pep_sequence[concat_seq_name][human_exon]]
                        flags[concat_seq_name]['many_ortho2one_human'].append(to_fix)
                else: 
                    # no exon in this species
                    [pep, dna] = gapped_exon (alnmt_pep[human_exon], alnmt_dna[human_exon])

                if output_pep[concat_seq_name]: output_pep[concat_seq_name] += '-Z-'
                output_pep[concat_seq_name] += pep
                if output_dna[concat_seq_name]: output_dna[concat_seq_name] += '---Z---'
                output_dna[concat_seq_name] += dna
                   
            headers.append(concat_seq_name)

        # clean all gaps
        output_pep = strip_gaps(output_pep)
        output_dna = strip_gaps(output_dna)

        # >>>>>>>>>>>>>>>>>>
        # cleanup for the case when the match is not one-to-one
        # take out the slice of the alignment
        # re-align the slice
        # put the slice back in the alignment
        # figure out how to expand DNA too
        if 1:
            for seq_to_fix in flags.keys():
                 if flags[seq_to_fix].has_key('one_ortho2many_human'): 
                    for to_fix in flags[seq_to_fix]['one_ortho2many_human']:
                        exon_seq_name       = to_fix[0]
                        list_of_human_exons = to_fix[1]
                        output_pep = fix_one2many (cfg, acg, exon_seq_name, seq_to_fix, 
                                                   list_of_human_exons, seqid, alnmt_pep, output_pep)

        #exit (1)

 
        # >>>>>>>>>>>>>>>>>>
        # find place for best_afa -- for the moment can put it to the scratch space:
        #afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        afa_fnm  = 'test.afa'
        sorted_trivial_names = map(lambda species: trivial_name[species], sorted_species)
        sorted_seq_names = sort_names (sorted_trivial_names, output_pep)
        output_fasta (afa_fnm, sorted_seq_names, output_pep)
        print afa_fnm
        exit(1)
        # afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        afa_fnm  = 'test.nt.afa'
        output_fasta (afa_fnm, sorted_seq_names, output_dna)
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
    #[template_name, template_seq] = find_human_template(alnmt_pep[human_exon])
    #[pep, dna] = merged_sequence (template_seq, pep_sequence[concat_seq_name][human_exon], 
    #                               dna_sequence[concat_seq_name][human_exon])
 '''
