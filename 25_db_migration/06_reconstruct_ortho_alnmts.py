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
from alignment import * # C implementation of smith waterman
from   random  import  choice
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
            print msg, 
            print "seq length check failure ",  name, len(seq),  aln_length
            afa_fnm = msg+'.afa'
            output_fasta (afa_fnm, sequence.keys(), sequence)
            print afa_fnm
            #exit(1)
            return False
    return True

#########################################
def find_initial_pos (pepseq_pieces, remaining_indices):

    initial_pos = []
    for i in range(len(pepseq_pieces)):
        if not i in remaining_indices: continue
        for pos in range(len(pepseq_pieces[i])):
            if not pepseq_pieces[i][pos] == '-':
                initial_pos.append(pos)
                break
    return initial_pos
    
#########################################
def check_seq_overlap (template_seq, pep_seq_pieces, pep_seq_names, new_names_of_exons):
    
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
        
    new_names_of_exons = filter (lambda exon: exon not in seq_names_to_remove, new_names_of_exons)

    remaining_names    = filter (lambda exon: exon not in seq_names_to_remove, pep_seq_names)
    # what is the order in which these sequences map to human?
    if len(remaining_names) > 1:
        # sort the names of the remaininig pieces by their initial position in the map to human exon(s)
        remaining_indices = filter (lambda idx: not idx in to_delete, range(len(pep_seq_names)))
        initial_pos       = find_initial_pos(pep_seq_pieces, remaining_indices) # for ex [57, 25, 12]
        a                 = range(len(initial_pos))    # for example    [0, 1, 2]
        a.sort (lambda i,j: cmp(initial_pos[i], initial_pos[j]) ) # for example  [2, 1, 0]
        to_reorder = []
        for idx in  range(len(new_names_of_exons)):
            name = new_names_of_exons[idx]
            if name in remaining_names:
                to_reorder.append(idx)

        for i in range(len(remaining_names)):
            if i >= len(to_reorder): return new_names_of_exons

        for i in range(len(remaining_names)):
            idx = to_reorder[i]
            reordered_name = remaining_names[a[i]]
            new_names_of_exons[idx] = reordered_name
            
    return new_names_of_exons


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
        seq_stripped[name] = seq_stripped[name].replace("BZB", "-Z-")

    return seq_stripped

#########################################
def make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial, flank_length):

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

            reconst_ntseq = expand_pepseq (reconst_pepseq, exon_seqs[1:], flank_length)
            if reconst_ntseq: 
                sequence_dna[sequence_name] = reconst_ntseq
                dna_aln_length = len(reconst_ntseq)

    # strip common gaps
    sequence_stripped_pep = strip_gaps (sequence_pep)
    # strip common gaps
    sequence_stripped_dna = strip_gaps (sequence_dna)

    # TODO
    #if not sequence_stripped_pep:
    #    print "blah"
    #    exit(1)

    return [sequence_stripped_pep, sequence_stripped_dna]
#########################################
def parse_aln_name (name):
    fields     = name.split("_")
    exon_id    = int(fields[-2])
    exon_known = int(fields[-1])
    species    =  "_".join(fields[:-2])
    return [species, exon_id, exon_known]


#########################################
def check_notes_directory (cfg):
    
    directory = "{0}/notes".format(cfg.dir_path['afs_dumps'])
    if not os.path.exists(directory):
        try:
            os.makedirs(directory) 
        except:
            print "error making", directory
            exit(1)

    return directory

#########################################
def print_notes (cursor, cfg,  ensembl_db_name, output_pep, names_of_exons, sorted_seq_names, stable_id):

    # write to string
    out_string  = "% Notes to accompany the alignment of (tentative) orthologues\n"
    out_string += "%% for the canonical transcript of the human gene %s," % stable_id
    out_string += "\n" 
    out_string += "% The alignment shows the exons that correspond to human transcript,\n" 
    out_string += "% from the following genes: \n" 
    out_string += "%% %-30s  %-30s  %-30s \n" % ('species', 'common_name', 'gene_id')

    gene_id        = {}
    stable_gene_id = {}
    sci_name       = {}
    for name in output_pep.keys():
        for exon_name in names_of_exons[name]:
            [species, exon_id, exon_known] = parse_aln_name(exon_name)
            ortho_gene_id                  = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
            if gene_id.has_key(name):
                if not ortho_gene_id == gene_id[name]:
                    print "error concatenating exons: ", name, gene_id.has_key(name), ortho_gene_id
            else:
                gene_id[name]        = ortho_gene_id
                stable_gene_id[name] = gene2stable(cursor, ortho_gene_id, ensembl_db_name[species])
                sci_name[name]       = species

    sorted_seq_names = filter (lambda name: name in output_pep.keys(), sorted_seq_names)
    for name in sorted_seq_names:
        out_string += "%-30s %-30s  %-30s \n" % ( name, sci_name[name],  stable_gene_id[name])

    out_string += "\n" 
    out_string += "% The following exons were used in the alignment\n" 

    for name in sorted_seq_names:
        out_string += "%% %s   %s \n" % (name, stable_gene_id[name])
        out_string += "%% %15s  %15s  %15s     %6s  %10s  %s \n" % \
           ('exon_id', 'gene_from', 'gene_to', 'coding', 'cannonical', 'source')
        for exon_name in names_of_exons[name]:
            [species, exon_id, exon_known] = parse_aln_name(exon_name)
            exon = get_exon (cursor, exon_id, exon_known, ensembl_db_name[species])
            if exon.is_known:
                exon_stable_id = exon2stable(cursor, exon_id, ensembl_db_name[species])
            else:
                exon_stable_id = 'anon'

            source = get_logic_name (cursor, exon.analysis_id,  ensembl_db_name[species])
            out_string += "  %15s  %15s  %15s     %-6s  %-10s  %s \n" % \
                (exon_stable_id, exon.start_in_gene, exon.end_in_gene,
                 exon.is_coding, exon.is_canonical, source)
               
    
    directory = check_notes_directory (cfg)
    notes_fnm = directory + '/'+stable_id+'.txt'
    print notes_fnm
    of = erropen (notes_fnm, "w")
    print >> of, out_string
    of.close()

    return True


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
def find_human_template(exon_alignment):
    
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
def correct_overlap_Z(peptide_alnmt, name, pos, global_bdry_position):
    #pdb.set_trace()
    #print name
    
    length       = len(peptide_alnmt[name])
    if pos+3 >= length: return
    triple       = not '-' in peptide_alnmt[name][pos:pos+3]
    insert_left  = peptide_alnmt[name][pos] == '-' 
    insert_left  = insert_left or  triple and (not pos or  peptide_alnmt[name][pos-1]=='-')
    insert_right = peptide_alnmt[name][pos+2] == '-'
    insert_right = insert_right or  triple and (not pos>=length-3 or peptide_alnmt[name][pos+3]=='-')
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
def splice_gaps(alnmt, to_splice):

    number_of_previous_splices = 0
    for pos in to_splice:
        pos_corrected = pos + number_of_previous_splices

        number_of_previous_splices += 1
        for name in alnmt.keys():
            temp = alnmt[name]
            alnmt[name] = temp[:pos_corrected]+'-'+temp[pos_corrected:]


#########################################
# what I really need is the alignment program that enforces boundaries ...
def boundary_cleanup (peptide_alnmt, sorted_seq_names): 

    global_bdry_position = []

    # check that no boundary marker is missing one '-'
    delimiter = re.compile("Z")
    
    for name, seq in peptide_alnmt.iteritems():
        to_splice = []
        for match in delimiter.finditer(seq):
            if not match: continue
            if match.start() and not seq[match.start()-1] == '-':
                to_splice.append(match.start())
            if  match.start()+1 < len(seq) and not seq[match.start()+1] == '-':
                to_splice.append(match.start()+1)
        splice_gaps(peptide_alnmt, to_splice)

    # use sorted species as anchors
    for sorted_seq_name in sorted_seq_names:
        pep_exons = peptide_alnmt[sorted_seq_name].split ('-Z-')
        bdry_position  = 0
        for exon_ct in range(1,len(pep_exons)):
            bdry_position += len(pep_exons[exon_ct-1])
            global_bdry_position.append(bdry_position)
            bdry_position += 3
 
        global_bdry_position = sorted(global_bdry_position)

        # see which boundaries can be re-aligned with anchor's Z
        move_Z (peptide_alnmt, global_bdry_position, sorted_seq_name)

        # find sequences that have amino acid overlapping with the 
        # exon boundary marker for our "anchor" sequence
        for name, seq in peptide_alnmt.iteritems():
            if name==sorted_seq_name: continue  # this is the anchor

            for pos in global_bdry_position:
                # this one has the bounday marker itself
                if 'Z' in peptide_alnmt[name][pos:pos+3]:
                    continue
                # this one is not a problem
                if peptide_alnmt[name][pos:pos+3] == '---':
                    continue
                # this needs correction
                correct_overlap_Z(peptide_alnmt, name, pos, global_bdry_position)

        # see which boundaries can be re-aligned with anchor's Z
        move_Z (peptide_alnmt, global_bdry_position, sorted_seq_name)
            
#########################################
def find_exon_boundaries (peptide_alignment):

    global_bdry_position = []
    local_bdry_position  = {}

    delimiter = re.compile("-Z-")
    for name, seq in peptide_alignment.iteritems():
        local_bdry_position[name] = []
        for match in delimiter.finditer(seq):
            if not match: continue
            pos = match.start()
            if not pos in global_bdry_position:
                global_bdry_position.append(pos)
            if not pos in local_bdry_position[name]:
                local_bdry_position[name].append(pos)

    global_bdry_position = sorted(global_bdry_position)

    return global_bdry_position, local_bdry_position

#########################################
def name2count (output_pep, names_of_exons):

    name_ct2exon_ct = {}
    exon_ct2name_ct = {}

    for name, seq in output_pep.iteritems():
        pep_exons = seq.split ('-Z-')
        exon_ct            = 0
        nonempty_exon_ct   = 0
        name_ct2exon_ct[name] = []
        exon_ct2name_ct[name] = []
        for pe in pep_exons:
            exon_ct2name_ct[name].append(-1)
            pe = pe.replace('-','')
            if pe:
                name_ct2exon_ct[name].append(exon_ct)
                exon_ct2name_ct[name][exon_ct] = nonempty_exon_ct
                nonempty_exon_ct += 1
            exon_ct += 1

        # sanity checking
        if not nonempty_exon_ct == len( names_of_exons[name]):
            print " foul:    nonempty_exon_ct={0}".format(nonempty_exon_ct),
            print " len(names_of_exons[name])={0}".format(len(names_of_exons[name]))
            print name
            return [{},{}]

    return [name_ct2exon_ct, exon_ct2name_ct]

#########################################
def exon_seq_check(exon_seqs, pep_aligned, species, exon_name, verbose=True):

    [exon_pep_seq, trsl_from, trsl_to, exon_left_flank,
     exon_right_flank, exon_dna_seq] = exon_seqs   

    if exon_pep_seq == pep_aligned.replace ('-', ''):
        return True

    if verbose:
        print " foul "
        print species
        print exon_name
        print "in db:   ", exon_pep_seq
        print "in almt: ", pep_aligned.replace ('-', '')
        print

    return False

#########################################
def check_exon_order (cursor, ensembl_db_name, output_pep, names_of_exons):

    order_ok = True

    for name, seq in output_pep.iteritems():
        
        pepseqs = map (lambda pep: pep.replace('-', ''), seq.split('Z'))

        exon_name_ct = -1
        for pepseq in pepseqs:
            if pepseq:
                exon_name_ct += 1
                exon_name = names_of_exons[name][exon_name_ct]
                [species, exon_id, exon_known] = parse_aln_name(exon_name)
                exon_seqs   = get_exon_seqs(cursor, exon_id, exon_known, ensembl_db_name[species])[1:]
                if not exon_seq_check(exon_seqs, pepseq, species, exon_name):
                    order_ok = False
    return order_ok


#########################################
def expand_protein_to_dna_alnmt (cursor, ensembl_db_name, cfg, acg, sorted_seq_names, 
                                 names_of_exons, alnmt_pep, output_pep, flank_length):

    output_dna         = {}
    # which exons correspond to which name?
    [name_ct2exon_ct, exon_ct2name_ct] = name2count (output_pep, names_of_exons)
    if not name_ct2exon_ct: 
        return output_dna

    # in the peptide alignment find exon positions that are  global exon boundaries, 
    # and the ones that are local
    [global_bdry_position, local_bdry_position] = find_exon_boundaries(output_pep)

    # for each sequence
    # 1) check whether there are global positions (inserts) that need to be accomodated
    # 2) cut according to local boundary positions
    # 3) expand the exon
    # 4) recalculate the position of the boundary in this exon's frame
    # 5) add extra space for the flanks
    output_dna = {}
    for name, full_aligned_pepseq in output_pep.iteritems():

        output_dna[name] = ""

        prev_local_pos = -3
        exon_ct        = -1
        for local_pos in local_bdry_position[name]+[len(full_aligned_pepseq)]:

            # find inserts
            inserts          = []
            inserts_relative = []
            for global_pos in global_bdry_position:
                if global_pos >= local_pos: break
                if global_pos > prev_local_pos and global_pos< local_pos:
                    inserts.append(global_pos)

            # get the dna seqence (aligned)
            pep_seq     = full_aligned_pepseq[prev_local_pos+3:local_pos]
            exon_ct    += 1
            name_ct     = exon_ct2name_ct[name][exon_ct]
            if (name_ct < 0):
                dna_aligned = expand_pepseq (pep_seq, [], flank_length)
            else:
                exon_seq_name = names_of_exons[name][name_ct]
                [species, exon_id, exon_known] = parse_aln_name(exon_seq_name)
                exon_seqs   = get_exon_seqs(cursor, exon_id, exon_known, ensembl_db_name[species])[1:]
                if  exon_seq_check (exon_seqs, pep_seq, species, exon_seq_name):
                    dna_aligned = expand_pepseq (pep_seq, exon_seqs, flank_length)
                else:
                    dna_aligned = expand_pepseq (pep_seq, [], flank_length)
                    
            # recalculate the insert position in the dna version
            if inserts:
                inserts_relative = map(lambda pos: flank_length + 3*(pos - prev_local_pos-3), inserts)

            # put in the gaps if/where needed
            prev_ins = 0
            for ins in inserts_relative:
                output_dna[name] += dna_aligned[prev_ins:ins]
                output_dna[name] += '-'*flank_length + '-Z-' + '-'*flank_length
                prev_ins          = ins+9
            output_dna[name] += dna_aligned[prev_ins:]

            if local_pos < len(full_aligned_pepseq):
                output_dna[name] += '-Z-'
 
            prev_local_pos = local_pos


    return output_dna
   
#########################################
def check_Z_right(pepseq, pound=False):
    
    if not pepseq:
        return pepseq

    if pepseq[-1] == 'Z': 
        return pepseq[:-1]

    if pound:
        pattern = re.compile('Z[\-\#]*$')
    else:
        pattern = re.compile('Z[\-]*$')
        
    match   = re.search(pattern, pepseq)
    if (match):
        new_pepseq  = pepseq[:match.start()]
        if pound:
            new_pepseq += '#'
        else:
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
def decorate_and_concatenate (pepseqs):
    decorated_seq = ""
    count = 1
    for  pepseq in pepseqs:
        padded_count = "{0:03d}".format(count)
        decorated_seq += 'B'+padded_count+pepseq+'Z'
        count += 1

    return decorated_seq

########################################
def realign_slice_clustal_omega (pep_slice, seq_to_fix, pep_seq_pieces):

    new_pep_slice  = pep_slice

    #afa_fnm  = 'slice.afa'
    #output_fasta (afa_fnm, pep_slice.keys(), pep_slice)
    #print afa_fnm
    # exit(1)

    return new_pep_slice
    

########################################
def realign_slice (pep_slice, seq_to_fix, pep_seq_pieces):

    new_pep_slice  = pep_slice

    human_seq      =  decorate_and_concatenate(pep_slice['human'].split('Z'))
    pep_seq_to_fix =  decorate_and_concatenate(pep_seq_pieces)


    [aligned_human, aligned_ortho] \
        = smith_waterman_context (human_seq, pep_seq_to_fix)


    aligned_ortho = check_Z_right(re.sub('[\dB\#]','-', aligned_ortho))
    aligned_human = check_Z_right(re.sub('[\dB]','#', aligned_human), pound=True)

    if not aligned_ortho[-1] == '-': # so the boundary mark ends up being  '-Z-'
        aligned_ortho += '-'
        aligned_human += '#'


    # put the gaps into the remaining sequences
    new_pep_slice = {}
    new_pep_slice[seq_to_fix] = aligned_ortho
    for name in pep_slice.keys():            
        if (name== seq_to_fix or name == 'human'):
            pass
        else:
            new_pep_slice[name] = pep_slice[name]

    for pos in range(len(aligned_human)):
        if not aligned_human[pos] == '#': continue
        for name in pep_slice.keys():            
            if (name== seq_to_fix or name == 'human'):
                pass
            else:
                if pos == len(new_pep_slice[name]):
                    new_pep_slice[name] +=  '-'
                else:
                    temp = new_pep_slice[name] 
                    new_pep_slice[name] = temp[:pos]+'-'+temp[pos:]

    new_pep_slice['human'] = aligned_human.replace('#','-')

    return new_pep_slice

########################################
def fix_one2many (cursor, ensembl_db_name, cfg, acg, sorted_seq_names, canonical_human_exons, human_exon_map,
                  names_of_exons, seq_to_fix, overlapping_maps, alnmt_pep, output_pep):


    count = 0

    if not overlapping_maps: return [output_pep, names_of_exons]

    if not output_pep.has_key('human'): return [output_pep, names_of_exons]


    # 
    new_alignment_pep  = {}
    new_names_of_exons = []
    for human_exon in canonical_human_exons:
        has_map = False
        for ortho_exon in names_of_exons[seq_to_fix]:
            if alnmt_pep[human_exon].has_key(ortho_exon):
                if not ortho_exon in new_names_of_exons:
                    new_names_of_exons.append(ortho_exon)
                has_map = True

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
        for human_exon in human_exons:
            [template_name, template_seq]  = find_human_template(alnmt_pep[human_exon])
            sequence_pieces = []
            seq_piece_names = []
            for exon_seq_name in ortho_exons:
                if not alnmt_pep[human_exon].has_key(exon_seq_name): continue
                sequence_pieces.append(alnmt_pep[human_exon][exon_seq_name])
                seq_piece_names.append(exon_seq_name)
                new_names_of_exons = check_seq_overlap(template_seq, sequence_pieces, seq_piece_names, new_names_of_exons)

        # join sequences that are deemed to be ok
        pep_seq_pieces = [] 
        for ortho_exon in new_names_of_exons:
            for human_exon in human_exons:
                if alnmt_pep[human_exon].has_key(ortho_exon):
                    pep_seq_pieces.append( alnmt_pep[human_exon][ortho_exon].replace("-", "") )
                    break
     
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
            return [output_pep, names_of_exons]

        ####################################
        # find the slice position
        pep_slice_start = exon_aln_start_pep[smallest_id]
        pep_slice_end   = exon_aln_end_pep  [largest_id]

        ####################################
        # if the slice size becomes comparable to the alignment size, drop the sequence
        if  number_of_human_exons > 3 and \
                float(pep_slice_end-pep_slice_start)/len(output_pep['human']) > 0.3 or 'elephant' in seq_to_fix:
            #print "deleting", seq_to_fix
            del output_pep[seq_to_fix]
            return [output_pep, names_of_exons]

        ####################################
        # cut out the slice
        pep_slice = {}
        for name in current_pep.keys():            
            pep_slice[name] = current_pep[name][pep_slice_start:pep_slice_end]

        # slice realign
        new_pep_slice = realign_slice (pep_slice, seq_to_fix, pep_seq_pieces)

        # strip gaps and output
        boundary_cleanup(new_pep_slice, sorted_seq_names)

        if not check_seq_length(new_pep_slice, "new_pep_slice"):
            del output_pep[seq_to_fix]
            return [output_pep, names_of_exons]

        #################################### 
        # replace the slice with the re-aligned one
        new_alignment_pep = {}
        for name, pepseq in current_pep.iteritems():
            new_alignment_pep[name]  = pepseq[:pep_slice_start]
            new_alignment_pep[name] += new_pep_slice[name]
            new_alignment_pep[name] += pepseq[pep_slice_end:] 

        if not check_seq_length(new_alignment_pep, "new_alignment_pep"):
            del output_pep[seq_to_fix]
            return [output_pep, names_of_exons]

        boundary_cleanup(new_alignment_pep, new_alignment_pep.keys())
        new_alignment_pep = strip_gaps(new_alignment_pep)
        
        current_pep = new_alignment_pep


    pep_exons = new_alignment_pep[seq_to_fix].split ('Z') 
    exon_ct   = 0
    for pe in pep_exons:
        pe = pe.replace('-','')
        if pe: 
            exon_ct += 1

    if not exon_ct == len(new_names_of_exons):
        print seq_to_fix, 'length mismatch'
        print "lengths:", exon_ct, len(new_names_of_exons)
        print "\n".join( map (lambda seq: seq.replace('-','') + " *** ", pep_exons) )
        print new_names_of_exons
        return [output_pep, names_of_exons]

    output_pep = new_alignment_pep
    names_of_exons[seq_to_fix] = new_names_of_exons

    return [output_pep, names_of_exons]


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
def remove_ghosts (output_pep, names_of_exons):

    delimiter      = re.compile("Z")

    #find positions that are all-gaps-but-one
    for name, seq in output_pep.iteritems():
        if name=='human': continue

        # find Z positions
        start    = 0
        prev_end = 0
        exon_ct  = 0
        removed_exons = []
        for match in delimiter.finditer(seq):
            start    = prev_end 
            end      = match.start()
            prev_end = match.end()

            pepseq = seq[start:end].replace('-','')
            if not len(pepseq): continue

            exon_ct +=1

            # are perhaps all other seqs gap in that range?
            is_ghost = True
            for pos in range(start,end):
                for name2, seq2 in output_pep.iteritems():
                    if name2==name: continue
                    if not seq2[pos] == '-':
                        is_ghost = False
                        break
                        
            # if yes replace with gaps 
            if not is_ghost: continue

            # check the Z itself
            if start:
                remove_start = True
                for name2, seq2 in output_pep.iteritems():
                    if name2==name: continue
                    if not seq2[start-1] == '-':
                        remove_start = False
                        break
                if remove_start: start -= 1

            if end < len(seq):
                remove_end = True
                for name2, seq2 in output_pep.iteritems():
                    if name2==name: continue
                    if not seq2[end] == '-':
                        remove_end = False
                        break
                if remove_end: end += 1

            # remove
            temp = output_pep[name]
            output_pep[name] = temp[:start]+'-'*(end-start)+temp[end:]
            # which exon is it
            removed_exons.append(exon_ct)
            # remove the name from the list
            new_names = []
            for exon_ct in range(len(names_of_exons[name])):
                if not exon_ct in removed_exons:
                    new_names.append(names_of_exons[name][exon_ct])
            names_of_exons[name] = new_names

    # strip gaps
    output_pep = strip_gaps(output_pep)

    return [output_pep, names_of_exons]



#########################################
#########################################
#########################################
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
    gene_list.reverse()
    for gene_id in gene_list:
    #for gene_id in [387298]:
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)

        gene_ct += 1

        afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        #afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        #if (os.path.exists(afa_fnm) and os.path.getsize(afa_fnm) > 0):
        #     continue
        notes_fnm = "{0}/notes/{1}.txt".format(cfg.dir_path['afs_dumps'], stable_id)
        if (os.path.exists(notes_fnm) and os.path.getsize(notes_fnm) > 0):
            continue

        if verbose: 
            print gene_id, stable_id, get_description (cursor, gene_id)
        elif (not gene_ct%100): 
            print gene_ct, "out of ", len(gene_list)

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
        bad_exons = []
        for human_exon in canonical_human_exons:
            [alnmt_pep[human_exon], alnmt_dna[human_exon]]  = \
                make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial, flank_length)   
            if not alnmt_pep[human_exon]: 
               bad_exons.append(human_exon)

        canonical_human_exons = filter (lambda x: not x in bad_exons, canonical_human_exons)
        
        # >>>>>>>>>>>>>>>>>>
        # bail out if there is a problem
        if not canonical_human_exons: 
            #print "\t botched"
            continue

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
                
        # >>>>>>>>>>>>>>>>>>
        names_of_exons = {}
        human_exon_map     = {}
        for human_exon in canonical_human_exons:
            for exon_seq_name in alnmt_pep[human_exon].keys():
                concat_seq_name = parent_seq_name[exon_seq_name]

                if not names_of_exons.has_key(concat_seq_name): names_of_exons[concat_seq_name] = []
                if not exon_seq_name in  names_of_exons[concat_seq_name]:
                    names_of_exons[concat_seq_name].append(exon_seq_name)
                    
                if not human_exon_map.has_key(concat_seq_name): human_exon_map[concat_seq_name] = {}
                if not human_exon_map[concat_seq_name].has_key(human_exon): 
                    human_exon_map[concat_seq_name][human_exon] = []
                human_exon_map[concat_seq_name][human_exon].append(exon_seq_name)


        # >>>>>>>>>>>>>>>>>>
        # flag the cases when one orthologue exon maps to many human (and vice versa) for later
        for concat_seq_name, concat_exons in names_of_exons.iteritems():
            overlapping_maps[concat_seq_name] = find_overlapping_maps (ortho_exon_to_human_exon, concat_exons, alnmt_pep)


        # >>>>>>>>>>>>>>>>>>
        # concatenate the aligned exons for each species, taking into account that the alignment
        # doesn't have to be one to one
        headers     = []
        output_pep  = {}
        output_dna  = {}

        for concat_seq_name in names_of_exons.keys():

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
        boundary_cleanup(output_pep, sorted_seq_names)
        output_pep = strip_gaps(output_pep)

        for seq_to_fix in overlapping_maps.keys():
            # fix_one2many changes bout output_pep and names_of_exons
            [output_pep, names_of_exons] = fix_one2many (cursor, ensembl_db_name, cfg, acg, sorted_seq_names, 
                                                         canonical_human_exons, human_exon_map, 
                                                         names_of_exons, seq_to_fix, 
                                                         overlapping_maps[seq_to_fix], 
                                                         alnmt_pep, output_pep)
 
            # we may have chosen to delete some sequences
            sorted_seq_names = sort_names (sorted_trivial_names['human'], output_pep)

        if not check_seq_length (output_pep, "ouput_pep"): 
            #print "length check failure"
            continue

        if (0):
            # >>>>>>>>>>>>>>>>>>
            boundary_cleanup(output_pep, sorted_seq_names)
            output_pep = strip_gaps(output_pep)

            # get rid of the ghost exons that do not correpond to anything in any other species
            #[output_pep, names_of_exons] = remove_ghosts(output_pep, names_of_exons)

            afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
            output_fasta (afa_fnm, sorted_seq_names, output_pep)
            print afa_fnm

            # >>>>>>>>>>>>>>>>>>
            output_dna = expand_protein_to_dna_alnmt (cursor, ensembl_db_name, cfg, acg, 
                                                      sorted_trivial_names, names_of_exons,  
                                                      alnmt_pep, output_pep, flank_length)
            if not output_dna:
                continue

            output_dna = strip_gaps(output_dna)

            afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
            output_fasta (afa_fnm, sorted_seq_names, output_dna)
            print afa_fnm

        #continue

        # notes to accompany the alignment:
        print_notes (cursor, cfg,  ensembl_db_name, output_pep, names_of_exons,  sorted_seq_names, stable_id)

    return 

#########################################
def get_theme_ids(cursor, theme_name):

    fnm = '../resources/'+theme_name+'.txt'
    if not os.path.exists(fnm):
        print fnm, "not found"
        exit(1)

    if not os.path.getsize(fnm) > 0:
        print fnm, "empty"
        exit(1)
        
    inf = erropen(fnm, "r")
    gene_ids = []
    for line in inf:
        line.rstrip()
        [stable_id, name] = line.split("\t")
        qry = "select gene_id, description from gene where stable_id='%s'"% stable_id
        rows = search_db (cursor, qry)
        if not rows: continue
        gene_ids.append(rows[0][0])
    inf.close()

    return gene_ids


#########################################
def main():
    
    no_threads = 15

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
    #gene_list                      =  get_theme_ids(cursor, 'wnt_pathway')
    cursor.close()
    db.close()

    parallelize (no_threads, make_alignments, gene_list, [local_db, ensembl_db_name])
    
    return True


#########################################
if __name__ == '__main__':
    main()





'''
    #for random_check in range(100):
    #    gene_id = choice(gene_list)
    #for gene_id in [412667]: #  wls   
    #for gene_id in [416066]:  #  BRCA1   
    #for gene_id in [378768]:  #  p53
    #for gene_id in [389337]: #inositol polyphosphate-4-phosphatase
    #for gene_id in [418590]: # titin
    #for gene_id in [412362]: # complement component 1,s subcomponent
'''
