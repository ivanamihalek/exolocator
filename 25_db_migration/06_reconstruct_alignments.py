#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import pdb
#pdb.set_trace()

import MySQLdb, commands, re, os

from el_utils.mysql   import  connect_to_mysql, connect_to_db
from el_utils.mysql   import  switch_to_db,  search_db, store_or_update

from el_utils.ensembl import  *
from el_utils.utils   import  erropen, output_fasta, input_fasta
from el_utils.map     import  Map, get_maps
from el_utils.tree    import  species_sort
from el_utils.ncbi    import  taxid2trivial
from el_utils.almt_cmd_generator import AlignmentCommandGenerator
from el_utils.config_reader      import ConfigurationReader
from el_utils.translation        import phase2offset, translation_bounds, crop_dna, translate
from el_utils.threads import parallelize
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
def check_seq_length(sequence, msg):

    if not sequence.values():
        return False
    aln_length = len(sequence.values()[0])
    if not aln_length:
        return False
    for name, seq in sequence.iteritems():
        if not len(seq) == aln_length:
            print msg
            print "seq length check failure: ",  name, len(seq),  aln_length
            afa_fnm  = 'test.afa'
            output_fasta (afa_fnm, sequence.keys(), sequence)
            print afa_fnm
            exit(1)
            return False
    return True

#########################################
def check_seq_overlap (template_seq, pep_seq_pieces, pep_seq_names):
    
    seq_names_to_remove = []

    template_length = len(template_seq)

    for piece in pep_seq_pieces:
        if (not len(piece) == template_length):
            print "length mismatch for aligned exons (?)"
            return []

    # check whether any two pieces overlap
    overlap = []
    for pos in range(template_length): 

        for i in range(len(pep_seq_pieces)):
            if (pep_seq_pieces[i][pos] == '-'): continue

            for j in range (i+1, len(pep_seq_pieces)):
                if (pep_seq_pieces[j][pos] == '-'): continue
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
        if ( fract_identity (template_seq, pep_seq_pieces[i]) < 
             fract_identity (template_seq, pep_seq_pieces[j]) ):
            to_delete.append(i)
        else:
            to_delete.append(j)

    for i in to_delete:
        seq_names_to_remove.append(pep_seq_names[i])

    return seq_names_to_remove


#########################################
def align_nucseq_by_pepseq(aligned_pepseq, nucseq):
    if (not len(aligned_pepseq.replace('-',''))*3 == len(nucseq)):
        print aligned_pepseq.replace('-','')
        print "length mismatch: ", len(aligned_pepseq.replace('-',''))*3, len(nucseq)
        return ""
    codon = iter(map(''.join, zip(*[iter(nucseq)]*3)))
    aligned_nucseq = ''.join(('---' if c=='-' else next(codon) for c in aligned_pepseq))
    return aligned_nucseq

#########################################
def aligned_nt_length ( aligned_pepseq_template,  aligned_pepseq, flank_length):

    length = 0 
    for pos in range(len(aligned_pepseq) ):

        if aligned_pepseq_template[pos]=='Z':
            length += (2*flank_length+1)
        else:
            length += 3
       
    return  length

#########################################
def align_nucseq_by_pepseq_w_template( aligned_pepseq_template,  aligned_pepseq, nucseq, flank_length):

    if (not len(aligned_pepseq.replace('-',''))*3 == len(nucseq)):
        print aligned_pepseq.replace('-','')
        print "length mismatch: ", len(aligned_pepseq.replace('-',''))*3, len(nucseq)
        return ""

    codon         = iter(map(''.join, zip(*[iter(nucseq)]*3)))

    aligned_nucseq = ""
    for pos in range(len(aligned_pepseq) ):

        if aligned_pepseq_template[pos]=='Z':
            if (aligned_pepseq[pos] == '-'):
                aligned_nucseq += "-"*(2*flank_length+1)
            else:
                aligned_nucseq += "-"*(2*flank_length-2)
                aligned_nucseq += next(codon)
        elif (aligned_pepseq[pos] == '-'):
            aligned_nucseq += "---"
        else:
            aligned_nucseq += next(codon)
       
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
def strip_gaps (sequence):

    
    seq_stripped = {}

    all_gaps = {}  

    if not check_seq_length(sequence, 'strip gaps:'): 
        return sequence
    
    aln_length = len(sequence.itervalues().next())

    if aln_length is None or aln_length==0:
        return sequence
    
    for name, seq in sequence.iteritems():
        if not len(seq): 
            continue
        sequence[name] = seq.replace("-Z-", "BZB")


    for pos in range(aln_length):
        all_gaps[pos] = True
        for name, seq in sequence.iteritems():
            if not len(seq): 
                continue
            if (not seq[pos]=='-'):
                all_gaps[pos] = False
                break


    for name, seq in sequence.iteritems():
        if not len(seq): 
            continue
        seq_stripped[name] = ""
        for pos in range(aln_length):
            if all_gaps[pos]: continue
            seq_stripped[name] += seq[pos]

    for name, seq in seq_stripped.iteritems():
        if not len(seq): 
            continue
        seq_stripped[name] = seq.replace("BZB", "-Z-")

    return seq_stripped


#########################################
def expand_delimiter(sequence, delim, delim_expanded):
    
    seq_new = {}
    for name, seq in sequence.iteritems():
        seq_new[name] = seq.replace(delim, delim_expanded)

    return seq_new

#########################################
def make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial):

    sequence_pep = {}
    sequence_dna = {}
    shortest_l = -1 # Uninitialized  leading padding length
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
            #print " ! ", pepseq
            #print " ! ", pepseq2
            #exit (1)
            continue
            
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
    
    ret = filter (lambda spec: 'homo_sapiens' in spec, exon_alignment.keys())
    if not ret:
        return ["", ""]
    template_name = ret[0]
    template_seq  = exon_alignment[template_name]

    return [template_name, template_seq]

#########################################
def sort_names (sorted_species, alignment):

    sorted_names = []
    for species in sorted_species:
        for seq_name  in alignment.keys():

            if seq_name[-1].isdigit():
                aux = seq_name.split("_")
                base_name = "_".join(aux[:-1])

                if base_name [-1].isdigit():
                    aux = base.split("_")
                    base_name = "_".join(aux[:-1])
            else:
                base_name = seq_name
            if (species == base_name):
                sorted_names.append(seq_name)
    return sorted_names
                     
            
#########################################
def move_Z (peptide_alnmt, global_bdry_position, template_name):

    # see which boundaries can be re-aligned with human
    for name, seq in peptide_alnmt.iteritems():
        if name==template_name: continue
        pep_exons = peptide_alnmt[name].split ('-Z-')
        
        bdry_position  = 0
        for exon_ct in range(1,len(pep_exons)):
            bdry_position += len(pep_exons[exon_ct-1])
            if not bdry_position in global_bdry_position:
                for pos in global_bdry_position:
                    if bdry_position < pos:
                        gaps ='-'*(pos-bdry_position)
                        if peptide_alnmt[name][bdry_position+3:pos+3] ==gaps:
                            pepseq  = peptide_alnmt[name][:bdry_position]
                            pepseq += gaps+'-Z-'
                            pepseq += peptide_alnmt[name][pos+3:]
                            peptide_alnmt[name] = pepseq
                            bdry_new = False
                            break
                    else:
                       gaps ='-'*(bdry_position - pos) 
                       if peptide_alnmt[name][pos:bdry_position] ==gaps:
                            pepseq  = peptide_alnmt[name][:pos]
                            pepseq += '-Z-'+gaps
                            pepseq += peptide_alnmt[name][bdry_position+3:]
                            peptide_alnmt[name] = pepseq
                            bdry_new = False
                            break
            bdry_position += 3
    

#########################################
def  correct_overlap_Z(peptide_alnmt, name, pos, global_bdry_position):

    length       = len(peptide_alnmt[name])
    triple       = not '-' in peptide_alnmt[name][pos:pos+3]
    insert_left  = peptide_alnmt[name][pos] == '-' 
    insert_left  = insert_left or  triple and (not pos or  peptide_alnmt[name][pos-1]=='-')
    insert_right = peptide_alnmt[name][pos+2] == '-'
    insert_right = insert_right or  triple and (not pos>=length-3 or peptide_alnmt[name][pos+3]=='-')
    #pdb.set_trace()
    if insert_left:
        # insert gaps on the left
        if triple:
            insert = '---'
            insert_start = pos

        elif peptide_alnmt[name][pos:pos+2] == '--':
            insert = '-'
            insert_start = pos+2
        else:
            insert = '--'
            insert_start = pos+1

        for name2, seq2 in peptide_alnmt.iteritems():
            if peptide_alnmt[name2][pos:pos+3] == '-Z-':
                pepseq  = seq2[:pos+3]
                pepseq += insert
                pepseq += seq2[pos+3:]
                peptide_alnmt[name2] = pepseq
            else:
                pepseq  = seq2[:insert_start]
                pepseq += insert
                pepseq += seq2[insert_start:]
                peptide_alnmt[name2] = pepseq
        # also: shift the boundary positions
        for i in range(len(global_bdry_position) ):
            if global_bdry_position[i] <= pos: continue
            global_bdry_position[i] += len(insert)

    elif insert_right:

        # insert gaps on the right
        if triple:
            insert = '---'
            insert_start = pos+3
        elif peptide_alnmt[name][pos+1:pos+3] == '--':
            insert = '-'
            insert_start = pos+1
        else:
            insert = '--'
            insert_start = pos+2


        for name2, seq2 in peptide_alnmt.iteritems():
            if peptide_alnmt[name2][pos:pos+3] == '-Z-':
                pepseq  = seq2[:pos]
                pepseq += insert
                pepseq += seq2[pos:]
                peptide_alnmt[name2] = pepseq
            else:
                pepseq  = seq2[:insert_start]
                pepseq += insert
                pepseq += seq2[insert_start:]
                peptide_alnmt[name2] = pepseq

        # also: shift the boundary positions
        for i in range(len(global_bdry_position) ):
            if global_bdry_position[i] < pos: continue
            global_bdry_position[i] += len(insert)


#########################################
# what I really need is the alignment program that enforces boundaries ...
def boundary_cleanup(peptide_alnmt, sorted_seq_names): 

    global_bdry_position = []
    # now use remaining species as anchors
    for sorted_seq_name in sorted_seq_names:
        pep_exons = peptide_alnmt[sorted_seq_name].split ('-Z-')
        bdry_position  = 0
        for exon_ct in range(1,len(pep_exons)):
            bdry_position += len(pep_exons[exon_ct-1])
            global_bdry_position.append(bdry_position)
            bdry_position += 3
 
        global_bdry_position = sorted(global_bdry_position)
        # see which boundaries can be re-aligned with human
        move_Z (peptide_alnmt, global_bdry_position, sorted_seq_name)

        for name, seq in peptide_alnmt.iteritems():
            if name==sorted_seq_name: continue

            for pos in global_bdry_position:

                if 'Z' in peptide_alnmt[name][pos:pos+3]:
                    continue

                if peptide_alnmt[name][pos:pos+3] == '---':
                    continue

                correct_overlap_Z(peptide_alnmt, name, pos, global_bdry_position)

        # see which boundaries can be re-aligned with human
        move_Z (peptide_alnmt, global_bdry_position, sorted_seq_name)
            
#########################################
def expand (aligned_peptide, old_aligned_dna, flank_length):   
    output_dna = {}

    aligned_exons          = aligned_peptide.split('Z')
    old_aligned_exons_dna  = old_aligned_dna.split('Z')

    if  not len(aligned_exons) ==  len(old_aligned_exons_dna):
        print
        print aligned_peptide
        print
        print old_aligned_dna
        print 
        print aligned_exons
        print
        print old_aligned_exons_dna
        print
        exit(1)

    expected_lenth = 0
    aligned_dna = ""
    for ct in range(len(aligned_exons)):

        aligned_peptide  = aligned_exons[ct]
        peptide_stripped = aligned_peptide.replace("-","")

        
        nt_length = 3*len(aligned_peptide)+2*flank_length

        if (not peptide_stripped):
            aligned_exon = "-"*nt_length
    
        else:
            dna                = old_aligned_exons_dna[ct].replace("-","")

            left_flank_pattern  = re.compile('[atcgn]+[ACTGN]')
            right_flank_pattern = re.compile('[ACTGN][actgn]+')


            left_flank = ""
            for match in left_flank_pattern.finditer(dna):
                start = match.start()
                end   = match.end() - 1 # for the capital letter in the regex
                left_flank =  dna[start:end]
            right_flank = ""
            for match in right_flank_pattern.finditer(dna):
                start = match.start()+ 1
                end   = match.end()
                right_flank = dna[start:end]

            if ( left_flank):
                len_left = len(left_flank)
            else:
                len_left = 0

            if ( right_flank):
                len_right = -len(right_flank)
            else:
                len_right = len(dna)

            left_flank   = left_flank.rjust (flank_length, '-')
            right_flank  = right_flank.ljust(flank_length, '-')
            aligned_exon = left_flank + align_nucseq_by_pepseq(aligned_peptide, dna[len_left:len_right]) + right_flank
             
        if aligned_dna: aligned_dna += 'Z'

        if not len(aligned_exon) == nt_length:
            print aligned_peptide
            print aligned_exon
            exit(1)

        aligned_dna +=  aligned_exon
       
 
    return  aligned_dna

#########################################
def find_bounds (sequence, delim):
    delimiter = re.compile(delim)

    # slice position in the peptide alignment
    exon_aln_start_pep = []
    exon_aln_end_pep   = []
    start              = 0
    prev_end           = 0
    for match in delimiter.finditer(sequence):
        start = prev_end 
        end   = match.start()
        exon_aln_start_pep.append(start)
        exon_aln_end_pep.append(end)
        prev_end   =  match.end()

    start = prev_end
    end   = len(sequence)
    exon_aln_start_pep.append(start)
    exon_aln_end_pep.append(end)

    return     [exon_aln_start_pep, exon_aln_end_pep]    

#########################################
def fill_hack(pep_slice, protected_name, exon_aln_start, exon_aln_end):

    filled_w_human_seq = {}
    offset = exon_aln_start[0]
    alpha  = re.compile('[A-Z]')

    for name, seq in pep_slice.iteritems():

        if name=='human':        continue
        if name==protected_name: continue

        for ct in range (len(exon_aln_start)):
            start = exon_aln_start[ct] - offset
            end   = exon_aln_end[ct] - offset
            if re.search(alpha, seq[start:end]):
                 continue
            pep_slice[name]  = seq[:start]
            pep_slice[name] += pep_slice['human'][start:end]
            pep_slice[name] += seq[end:]
            if not  filled_w_human_seq.has_key(name): 
                filled_w_human_seq[name] = []
            filled_w_human_seq[name].append(ct)
    
    return filled_w_human_seq


#########################################
def expand_protein_to_dna_alnmt (cursor, ensembl_db_name, cfg, acg, sorted_trivial_names, 
                                 concatenated_exons, alnmt_pep, output_pep, flank_length):

    output_dna         = {}

    dna_wo_exon_bdries = {}
    left_flank         = {}
    right_flank        = {}
    local_exon_position  = {}
    global_exon_position = []
    for concat_name, seq in output_pep.iteritems():

        pep_exons = seq.split ('-Z-')

        print len(pep_exons), len(concatenated_exons[concat_name])
        continue

        conc_exons_names = []
        dna     = {}
        peptide = {}

        for human_exon, list_of_exon_names in concatenated_exons[concat_name].iteritems():
            print human_exon.exon_id
            print list_of_exon_names
            for exon_name in list_of_exon_names:
                
                tokens     = exon_name.split('_')
                species    = "_".join(tokens[:-2])
                exon_id    = int(tokens[-2])
                exon_known = int(tokens[-1])
                exon_seqs  = get_exon_seqs(cursor, exon_id, exon_known, ensembl_db_name[species])[1:]
                print exon_name
                print exon_seqs
                print
            print

        print len(pep_exons), len(concatenated_exons[concat_name])

        exit(1)

        for human_exon, list_of_exon_names in concatenated_exons[concat_name].iteritems():
            for exon_name in list_of_exon_names:
                if not exon_name in conc_exons_names:
                    conc_exons_names.append(exon_name)
                    peptide[exon_name]  =  alnmt_pep[human_exon][exon_name].replace ('-', '')

        local_exon_position[concat_name] = []
        left_flank [concat_name] = {}
        right_flank[concat_name] = {}
        dna_wo_exon_bdries[concat_name] = ""
        for pep_exon in pep_exons:
            pepe =  pep_exon.replace ('-', '')
            if not pepe: 
                dna_seq= '-'*(2*flank_length+3*len(pep_exon))
            else:
                matching_exon_names = filter (lambda  exon_name: 
                                              peptide[exon_name] == pepe, 
                                              conc_exons_names)
                if not matching_exon_names: 
                    dna_seq= '-'*(2*flank_length+3*len(pep_exon))
                else:

                    # what if I have two matching seqs (?)
                    tokens     = matching_exon_names[0].split('_')
                    species    = "_".join(tokens[:-2])
                    exon_id    = int(tokens[-2])
                    exon_known = int(tokens[-1])
                    exon_seqs  = get_exon_seqs(cursor, exon_id, exon_known, ensembl_db_name[species])[1:]
                    dna_seq    = expand_pepseq (pep_exon, exon_seqs)

            pos        = len(dna_wo_exon_bdries[concat_name])

            dna_wo_exon_bdries[concat_name] += dna_seq[flank_length:-flank_length]
            left_flank [concat_name][pos]    = dna_seq[:flank_length]
            right_flank[concat_name][pos]    = dna_seq[-flank_length:]

            if not pos: continue
            local_exon_position[concat_name].append(pos) 
            if not pos  in global_exon_position:
                global_exon_position.append(pos)
        
    exit(1)
    
    global_exon_position = sorted(global_exon_position)


    # comb one more time to insert the flanking regions
    for concat_name, seq in output_pep.iteritems():
        prev_pos = 0
        if  left_flank[concat_name].has_key(prev_pos):
            output_dna[concat_name] = left_flank[concat_name][prev_pos]
        else:
            output_dna[concat_name] = '-'*flank_length 

        for pos in global_exon_position:

            if pos in local_exon_position[concat_name]:
                boundary = '---'+'-Z-'+'---'
            else:
                boundary = ''
            output_dna[concat_name] += dna_wo_exon_bdries[concat_name][prev_pos:pos]
            output_dna[concat_name] += boundary
            prev_pos = pos

 

        output_dna[concat_name] += dna_wo_exon_bdries[concat_name][prev_pos:]
        if  right_flank[concat_name].has_key(prev_pos):
            output_dna[concat_name] += right_flank[concat_name][prev_pos]
        else:
            output_dna[concat_name] += '-'*flank_length 

    return output_dna


#########################################
def undo_fill_hack(pep_slice, filled_w_human_seq):

    # new bdry positions:
    [exon_start, exon_end] = find_bounds(pep_slice['human'], 'Z')

    for name in filled_w_human_seq.keys():
        for ct in filled_w_human_seq[name]:
            
            start = exon_start[ct]
            end   = exon_end  [ct]
            seq   = pep_slice [name]
            pep_slice[name]  = seq[:start]
            pep_slice[name] += '-'*(end-start)
            pep_slice[name] += seq[end:]

    return

  
#########################################
#########################################
#########################################
def check_duplicates (list):

    seen = []
    for a in list:
        if a == 'none': continue
        if a in seen:
            return True
        seen.append(a)

    return False
    
#########################################
def check_Z_right(pepseq):
    
    if not pepseq:
        return pepseq

    if pepseq[-1] == 'Z': 
        return pepseq[:-1]

    pattern = re.compile('Z\-*$')
    match   = re.search(pattern, pepseq)
    if (match):
        new_pepseq  = pepseq[:match.start()]
        new_pepseq += '-'
        new_pepseq += pepseq[match.start()+1:match.end()-1] + 'Z'
        new_pepseq += pepseq[match.end():]
        return new_pepseq[:-1]
    else:  
        #new_pepseq  = pepseq[:-1] + 'Z'
        return pepseq[:-1]
              
    
#########################################
def check_Z_left(pepseq):

    if not pepseq:
        return pepseq
    if  pepseq[0] == 'Z': 
        return pepseq[1:]

    pattern = re.compile('^\-*Z')
    match   = re.search(pattern, pepseq)

    if (match):
        new_pepseq  = pepseq[:match.start()]
        new_pepseq += 'Z'
        new_pepseq += pepseq[match.start()+1:match.end()-1] + '-'
        new_pepseq += pepseq[match.end():]        
        return new_pepseq[1:]
    else:  
        #new_pepseq  = 'Z' + pepseq[1:]
        return pepseq[1:]
    

#########################################
def fix_one2many (cfg, acg, sorted_trivial_names, canonical_human_exons, human_exon_map,
                  concatenated_exon_names, seq_to_fix, overlapping_maps, alnmt_pep, output_pep):


    count = 0

    if not overlapping_maps: return output_pep

    if not output_pep.has_key('human'): return output_pep

    verbose = seq_to_fix == 'wallaby' or  seq_to_fix == 'elephant'
    #verbose = False

    new_alignment_pep = {}

    names_of_concatenated_exons = []
    for human_exon in canonical_human_exons:
        has_map = False
        for ortho_exon in concatenated_exon_names[seq_to_fix]:
            if alnmt_pep[human_exon].has_key(ortho_exon):
                if not ortho_exon in names_of_concatenated_exons:
                    names_of_concatenated_exons.append(ortho_exon)
                has_map = True

    print
    print
    print seq_to_fix

    if verbose:
        print "\n".join( map (lambda seq: seq.replace('-','')+" *** ", output_pep[seq_to_fix].split('Z') ))
        print names_of_concatenated_exons
        print "lengths before:", len(output_pep[seq_to_fix].split('Z')),
        print len(concatenated_exon_names[seq_to_fix]), len(names_of_concatenated_exons)
        print 
        print "############################"
        print "flagged: "
        print
        for [human_exons, ortho_exons] in overlapping_maps:
            for human_exon in human_exons:
                print human_exon.exon_id
                for exon_seq_name in ortho_exons:
                    print "\t", exon_seq_name
                    for he2 in human_exons:
                        if alnmt_pep[he2].has_key(exon_seq_name):
                            pepseq = alnmt_pep[he2][exon_seq_name]
                            print "\t", pepseq.replace('-', '')
                            break
                print
        print "############################"


    # find sequential numbers of exons that we have in this story
    seqid = {}
    ct    = 0
    for human_exon in canonical_human_exons:
        seqid[human_exon] = ct
        ct += 1
    number_of_human_exons = ct

    current_pep = output_pep


    # for each unresolved "map"  cut out the slice and re-align
    for  [human_exons, ortho_exons] in overlapping_maps:

        exon_numbers = sorted(map (lambda ex: seqid[ex], human_exons))
        
        smallest_id  = exon_numbers[0]
        largest_id   = exon_numbers[-1]

        # remove duplicates ? can this happen
        ortho_exons = list(set(ortho_exons))

        # check sequence overlap, if several map to the same  human exon
        ortho_seq_to_remove = []
        for human_exon in human_exons:
            [template_name, template_seq]  = find_human_template(alnmt_pep[human_exon])
            sequence_pieces = []
            seq_piece_names = []
            for exon_seq_name in ortho_exons:
                if not alnmt_pep[human_exon].has_key(exon_seq_name): continue
                sequence_pieces.append(alnmt_pep[human_exon][exon_seq_name])
                seq_piece_names.append(exon_seq_name)
                ortho_seq_to_remove = check_seq_overlap(template_seq, sequence_pieces, seq_piece_names)
                if verbose: print "to remove: ", ortho_seq_to_remove

        names_of_concatenated_exons = filter (lambda exon: exon not in ortho_seq_to_remove, names_of_concatenated_exons)

        # join sequences that are deemed to be ok
        pep_seq_pieces = [] 
        for ortho_exon in ortho_exons:
            if ortho_exon in ortho_seq_to_remove: continue
            for human_exon in human_exons:
                if alnmt_pep[human_exon].has_key(ortho_exon):
                    pep_seq_pieces.append( alnmt_pep[human_exon][ortho_exon].replace("-", "") )
                    break
        pep_seq_to_fix = "JJJJJJJ".join(pep_seq_pieces)
     

        # pull  the slice out of the alignment
        # use human as the reference - in other species the boundaries might
        # be at different positions
        # 
        delimiter = re.compile("Z")

        # slice position in the peptide alignment
        exon_aln_start_pep = []
        exon_aln_end_pep   = []
        start              = 0
        prev_end           = 0
        for match in delimiter.finditer(current_pep['human']):
            start = prev_end 
            end   = match.start()
            exon_aln_start_pep.append(start)
            exon_aln_end_pep.append(end)
            prev_end   =  match.end()

        start = prev_end + 1
        end   = len(current_pep['human'])
        exon_aln_start_pep.append(start)
        exon_aln_end_pep.append(end)

        if len(exon_aln_start_pep) <= smallest_id or  len(exon_aln_start_pep) <= largest_id:
            print  len(exon_aln_start_pep), smallest_id, len(exon_aln_start_pep), largest_id
            return  current_pep

        ####################################
        pep_slice_start = exon_aln_start_pep[smallest_id]
        pep_slice_end   = exon_aln_end_pep  [largest_id]
        pep_slice = {}
        for name in current_pep.keys():            
            if (name== seq_to_fix):
                pep_slice[seq_to_fix] = pep_seq_to_fix
            else:
                pep_slice[name] = current_pep[name][pep_slice_start:pep_slice_end]


        ####################################
        exon_seq_name = ortho_exons[0]
        tmp_name  = exon_seq_name+"_"+"_".join(map (lambda x: str(x), exon_numbers))
        if len(tmp_name) > 50: tmp_name= tmp_name[0:50]
        afa_fnm = "{0}/{1}.afa".format(cfg.dir_path['scratch'], tmp_name)

        tmp_name = exon_seq_name+"_"+"_".join(map (lambda x: str(x), exon_numbers))
        if len(tmp_name) > 50: tmp_name= tmp_name[0:50]
        out_fnm = "{0}/{1}.out.afa".format(cfg.dir_path['scratch'], tmp_name)


        filled_w_human_seq = fill_hack(pep_slice, seq_to_fix, 
                                       exon_aln_start_pep[smallest_id:largest_id+1],
                                       exon_aln_end_pep[smallest_id:largest_id+1])
        B_hack = True
        for name in pep_slice.keys():
            if (len (pep_slice[name]) > 3000): # mafft chokes otherwise
                B_hack = False
                break
        for name in pep_slice.keys():
            if B_hack:
                pep_slice[name]= pep_slice[name].replace(  "-", "B"  )
            else:
                pep_slice[name]= pep_slice[name].replace("-Z-", "BZB")
                

        output_fasta (afa_fnm, pep_slice.keys(), pep_slice) 
        mafftcmd = acg.generate_mafft_command (afa_fnm)
        ret      = commands.getoutput(mafftcmd)

        pep_slice_realigned = {}
        pepseq = ""
        for line in ret.split('\n'):
            if '>' in line:
                if ( pepseq ):
                    pep_slice_realigned[name] = pepseq
                name = line.replace (">", "")
                name = name.replace (" ", "")
                pepseq = ""
            else:
                pepseq += line
        if pepseq: pep_slice_realigned[name] = pepseq

        for name in pep_slice_realigned.keys():
            pep_slice_realigned[name] = pep_slice_realigned[name].replace("B", "-")

        undo_fill_hack (pep_slice_realigned, filled_w_human_seq)

        # cleaning the JJJJJ thing
        pepseq = pep_slice_realigned[seq_to_fix]
        delimiter = re.compile("J+-*J+")
        for match in delimiter.finditer(pepseq):
            start = match.start()
            end   = match.end()
            length = end-start-1
            half   = length/2
            pepseq = pepseq[:start]+'-'*half + 'Z'+'-'*(length-half)+pepseq[end:]
        pep_slice_realigned[seq_to_fix] = pepseq

        old_aln_length = len(current_pep['human'])
        new_aln_length = pep_slice_start + len(pepseq) + (old_aln_length-pep_slice_end+1)

        # replace the slice with the re-aligned one
        new_alignment_pep = {}
        for name, pepseq in current_pep.iteritems():

            if name == seq_to_fix: #what in the hell is thich cecking for?
                new_alignment_pep[name]  = check_Z_right(pepseq[:pep_slice_start]) # will take off one char on the right
                if pep_slice_start>0:
                    new_alignment_pep[name] += 'Z'
                new_alignment_pep[name] += pep_slice_realigned[name]
                if largest_id < number_of_human_exons-1:
                    new_alignment_pep[name] += 'Z'
                new_alignment_pep[name] += check_Z_left( pepseq[pep_slice_end:]) # will take off one char on the left

            else:
                new_alignment_pep[name]  = pepseq[:pep_slice_start]
                new_alignment_pep[name] += pep_slice_realigned[name]
                new_alignment_pep[name] += pepseq[pep_slice_end:] 

                pep1 = pepseq.replace('Z','').replace('-','')
                pep2 = new_alignment_pep[name].replace('Z','').replace('-','')

                if ( not pep1 == pep2 ):
                    print name
                    print pep1
                    print "\n".join( map (lambda seq: seq.replace('-','') + " *** ", pepseq.split('Z')))
                    print "777777777777777"
                    print pep2
                    print "\n".join( map (lambda seq: seq.replace('-','') + " *** ", new_alignment_pep[name].split('Z')))
                    print "777777777777777"
                    print pepseq[pep_slice_end-1:pep_slice_end+2]
                    print mafftcmd
                    exit(1)
 
        if not check_seq_length(new_alignment_pep, "new_alignment_pep"):
            return current_pep

    


        #boundary_cleanup(new_alignment_pep, new_alignment_pep.keys())
        
        current_pep = new_alignment_pep

    pep_exons = new_alignment_pep[seq_to_fix].split ('Z') 
    exon_ct   = 0
    print len(names_of_concatenated_exons)
    for pe in pep_exons:
        pe = pe.replace('-','')
        if pe: 
            exon_ct += 1
            print exon_ct, pe

    if not exon_ct == len(names_of_concatenated_exons):
        print seq_to_fix, 'length mismatch'
        print "lengths:", exon_ct, len(names_of_concatenated_exons)
        print "\n".join( map (lambda seq: seq.replace('-','') + " *** ", pep_exons) )
        print names_of_concatenated_exons
        return output_pep

    #if verbose: exit(1)

    return new_alignment_pep


#########################################
def find_overlapping_maps (ortho_exon_to_human_exon, exon_seq_names, alnmt_pep):

    overlapping_maps = []

    # groups of exons that map onto each other in a non-trivial way
    join_groups = []
      
    for i in range (len(exon_seq_names)):

        exon_seq_name         = exon_seq_names[i]
        human_exons           = ortho_exon_to_human_exon[exon_seq_name]
        overlapping_human_set = False

        for j in range (i+1, len(exon_seq_names)):

            exon_seq_name_2  = exon_seq_names[j]  
            human_exons_2    = ortho_exon_to_human_exon[exon_seq_name_2]

            if set( human_exons) & set (human_exons_2) :
                overlapping_human_set = True
                group_found           = False

                for join_group in join_groups:
                    if exon_seq_name in join_group:
                        join_group.append(exon_seq_name_2)
                        group_found = True
                    elif  exon_seq_name_2 in join_group:
                        join_group.append(exon_seq_name)
                        group_found = True
                    if group_found: break
                if not group_found:
                    join_groups.append([exon_seq_name, exon_seq_name_2])

        if not overlapping_human_set and len(human_exons) > 1:
            overlapping_maps.append([human_exons, [exon_seq_name]])

    if (join_groups):
        for join_group in join_groups:
            human_exons = []
            for exon_seq_name in join_group:
                for hu_ex in ortho_exon_to_human_exon[exon_seq_name]:
                    if not hu_ex in human_exons:
                        human_exons.append(hu_ex)
            ortho_exons = list(set(join_group) )
            overlapping_maps.append([human_exons, ortho_exons])


    return overlapping_maps

#########################################
#########################################
#########################################
#########################################
def make_alignments ( gene_list, db_info):

    [local_db, ensembl_db_name] = db_info

    verbose      = True
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
    # in the afa headers use 'trivial' names for the species: cow, dog, pig, ...
    trivial_name   = translate_to_trivial(cursor, all_species)

    # walk the taxonomical tree, and sort the species according to
    # the (distance of) the last common ancestor
    sorted_species = {}
    sorted_trivial_names = {}
    for qry_species in ['homo_sapiens']:
        sorted_species[qry_species] = species_sort(cursor, all_species, qry_species)
        trivial = trivial_name[qry_species]
        sorted_trivial_names[trivial] = map(lambda species: trivial_name[species], sorted_species[qry_species])
        
    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
 

    # for each human gene
    gene_ct = 0
    #for gene_id in gene_list:
    #for gene_id in [412667]: #  wls   
    for gene_id in [378768]: #  p53
    #for gene_id in [389337]: #inositol polyphosphate-4-phosphatase
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)

        gene_ct += 1
        if verbose: 
            print gene_id, stable_id, get_description (cursor, gene_id)
        elif (not gene_ct%100): 
            print gene_ct, "out of ", len(gene_list)

        #afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        #if (os.path.exists(afa_fnm) and os.path.getsize(afa_fnm) > 0):
        #    continue
        # find all exons we are tracking in the database
        human_exons     = gene2exon_list(cursor, gene_id)
        canonical_human_exons = []
        for human_exon in human_exons:
            if not human_exon.is_canonical or  not human_exon.is_coding:
                continue
            canonical_human_exons.append(human_exon)

        # the exons are not guaranteed to be in order
        canonical_human_exons.sort(key=lambda exon: exon.start_in_gene)

        # >>>>>>>>>>>>>>>>>>
        # reconstruct the per-exon alignment with orthologues
        mitochondrial = is_mitochondrial(cursor, gene_id)
 
        alnmt_pep = {}
        alnmt_dna = {}
        has_a_map = True
        for human_exon in canonical_human_exons:
            [alnmt_pep[human_exon], alnmt_dna[human_exon]]  = \
                make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial)   
            if not alnmt_pep[human_exon]: 
                has_a_map = False
                break

        # >>>>>>>>>>>>>>>>>>
        # bail out if there is a problem
        if not has_a_map: continue

        # >>>>>>>>>>>>>>>>>>
        # do we have a sequence mapping to multiple human exons?
        ortho_exon_to_human_exon = {}
        for human_exon in canonical_human_exons:
            for exon_seq_name in alnmt_pep[human_exon].keys():
                if not ortho_exon_to_human_exon.has_key(exon_seq_name):
                    ortho_exon_to_human_exon[exon_seq_name] = [human_exon]
                else:
                    ortho_exon_to_human_exon[exon_seq_name].append(human_exon)
 
        # >>>>>>>>>>>>>>>>>>
        # find which species we have, and for how many exons
        # we may have two orthologues for the same species
        seq_name = {}
        overlapping_maps   = {}
        parent_seq_name = {}
        for human_exon in canonical_human_exons:
            for exon_seq_name, exon_seq in alnmt_pep[human_exon].iteritems():
                (species, exon_id, exon_known) = parse_aln_name(exon_seq_name)
                ortho_gene_id                  = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
                # gene --> name -- retrieve old name, or construct new one
                parent_seq_name[exon_seq_name] = get_name (seq_name, trivial_name[species], ortho_gene_id) 
                
        concatenated_exon_names = {}
        human_exon_map     = {}
        for human_exon in canonical_human_exons:
            for exon_seq_name in alnmt_pep[human_exon].keys():
                concat_seq_name = parent_seq_name[exon_seq_name]

                if not concatenated_exon_names.has_key(concat_seq_name): concatenated_exon_names[concat_seq_name] = []
                if not exon_seq_name in  concatenated_exon_names[concat_seq_name]:
                    concatenated_exon_names[concat_seq_name].append(exon_seq_name)
                    
                if not human_exon_map.has_key(concat_seq_name): human_exon_map[concat_seq_name] = {}
                if not human_exon_map[concat_seq_name].has_key(human_exon): 
                    human_exon_map[concat_seq_name][human_exon] = []
                human_exon_map[concat_seq_name][human_exon].append(exon_seq_name)


        # >>>>>>>>>>>>>>>>>>
        # flag the cases when one orthologue exon maps to many human (and vice versa) for later
        for concat_seq_name, concat_exons in concatenated_exon_names.iteritems():
            overlapping_maps[concat_seq_name] = find_overlapping_maps (ortho_exon_to_human_exon, concat_exons, alnmt_pep)


        # >>>>>>>>>>>>>>>>>>
        # concatenate the aligned exons for each species, taking into account that the alignment
        # doesn't have to be one to one
        headers     = []
        output_pep  = {}
        output_dna  = {}

        for concat_seq_name in concatenated_exon_names.keys():

            if not human_exon_map.has_key(concat_seq_name):  continue # this shouldn't happen but oh well

            output_pep[concat_seq_name] = ""

            # single out one ortho to mult human cases

            flagged_human_exons = []
            if overlapping_maps.has_key(concat_seq_name):
                for  [human_exons, ortho_exons] in overlapping_maps[concat_seq_name]:
                    if len(human_exons) == 1 and len(ortho_exons) == 1: continue 
                    flagged_human_exons = flagged_human_exons + human_exons

            for human_exon in canonical_human_exons:

                aln_length = len(alnmt_pep[human_exon].itervalues().next())
                pep = '-'*aln_length
                if human_exon in flagged_human_exons:
                    # one ortho seq maps to multiple human exons
                    pep = '-'*aln_length

                elif human_exon_map[concat_seq_name].has_key(human_exon):
                    # we have a neat one-to-one mapping
                    exon_seq_name = human_exon_map[concat_seq_name][human_exon][0]
                    pep = alnmt_pep[human_exon][exon_seq_name]
                    if not pep:  # what's this?
                        pep = '-'*aln_length
                else: 
                    # no exon in this species
                    pep =  '-'*aln_length
                   
                if output_pep[concat_seq_name]: output_pep[concat_seq_name] += '-Z-'
                output_pep[concat_seq_name] += pep

                if not pep:
                    print human_exon.exon_id
             
            headers.append(concat_seq_name)

        #########################################################
        # >>>>>>>>>>>>>>>>>>
        sorted_seq_names = sort_names (sorted_trivial_names['human'], output_pep)

        if ( 0 and os.path.exists(afa_fnm) and os.path.getsize(afa_fnm) > 0):
            output_pep = input_fasta(afa_fnm)
        else:
            for seq_to_fix in overlapping_maps.keys():
                output_pep = fix_one2many (cfg, acg, sorted_trivial_names, 
                                           canonical_human_exons, human_exon_map, concatenated_exon_names,
                                           seq_to_fix, overlapping_maps[seq_to_fix], alnmt_pep, output_pep)
                print
                print "wallaby after ", seq_to_fix
                print "\n".join( map (lambda seq: seq.replace('-','')+" *** ", output_pep['wallaby'].split('Z') ))
                print 
            exit(1)
            if not check_seq_length (output_pep, "ouput_pep"): 
                print "length check failure"
                continue

            # >>>>>>>>>>>>>>>>>>
            # find place for best_afa -- for the moment can put it to the scratch space:
            boundary_cleanup(output_pep, sorted_seq_names)
            output_pep = strip_gaps(output_pep)
            #afa_fnm  = 'test.afa'
            output_fasta (afa_fnm, sorted_seq_names, output_pep)
            print afa_fnm
            #continue
 
        # >>>>>>>>>>>>>>>>>>
        output_dna = expand_protein_to_dna_alnmt (cursor, ensembl_db_name, cfg, acg, 
                                                  sorted_trivial_names, concatenated_exon_names,  
                                                  alnmt_pep, output_pep, flank_length)
        #output_dna = strip_gaps(output_dna)

        afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        output_fasta (afa_fnm, sorted_seq_names, output_dna)
        print afa_fnm

        #continue
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
