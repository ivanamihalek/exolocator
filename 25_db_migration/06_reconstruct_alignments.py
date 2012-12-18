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
def merged_sequence (template_seq, sequence_pieces):
    
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
 

    # if not, go ahead and merge
    merged = "ZZZZZZZZZZZZ".join(sequence_pieces)

    return merged


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
    
    aln_length = None
    for name, seq in sequence.iteritems():
        aln_length = len(seq)
        break

    if aln_length is None or aln_length==0:
        return sequence

    for pos in range(aln_length):
        all_gaps[pos] = True
        for name, seq in sequence.iteritems():
            if not len(seq): continue
            if pos >=  len(seq): continue
            if (not seq[pos]=='-'):
                all_gaps[pos] = False
                break

    for name, seq in sequence.iteritems():
        if not len(seq): continue
        if pos >=  len(seq): continue
        seq_stripped[name] = ""
        for pos in range(aln_length):
            if all_gaps[pos]: continue
            seq_stripped[name] += seq[pos]

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
def flags_init( flags, seq_name, flag_name):
    if not flags.has_key(seq_name): 
        flags[seq_name] = {}
    if not flags[seq_name].has_key(flag_name):
        flags[seq_name][flag_name] = []
    return

#########################################
def cleanup_exon_boundaries(sequence):

    all_gaps = {}  
    
    aln_length = None
    for name, seq in sequence.iteritems():
        aln_length = len(seq)
        break
    if aln_length is None or aln_length==0:
        return sequence

    for pos in range(aln_length):
        all_gaps[pos] = True
        for name, seq in sequence.iteritems():
            if not len(seq): continue
            if pos >=  len(seq): continue
            if (not seq[pos]=='-' and not seq[pos]=='Z'):
                all_gaps[pos] = False
                break
   
    delimiter  = re.compile('ZZ+')

    
    for name, seq in sequence.iteritems():
        bdry_start = []
        bdry_end   = []
        for match in delimiter.finditer(seq):
            start = match.start()
            end   = match.end()
            bdry_start.append(start)
            bdry_end.append(end)

        for bdry_ct in range(len(bdry_start)):
            start = bdry_start[bdry_ct]
            end   = bdry_end[bdry_ct]
            # is there a single position in that range that is all gaps
            good_pos = None
            for pos in range(start,end):
                if all_gaps[pos]: 
                    good_pos = pos
                    break
            # keep Z there
            # otherwise keep any old Z
            if good_pos is None: good_pos = start
            new_seq = sequence[name][0:start]+'-'*(good_pos-start)+'Z'+'-'*(end-good_pos-1)+sequence[name][end:] 
            sequence[name] = new_seq

#########################################
def expand (aligned_peptide, unaligned_dna):   
    output_dna = {}

    left_flank_pattern  = re.compile('[atcgn]+[ACTGN]')
    right_flank_pattern = re.compile('[ACTGN][actgn]+')

    aligned_peptide = re.sub('Z+', 'Z',aligned_peptide)
    aligned_exons          = aligned_peptide.split('Z')
    aligned_exons_stripped = map (lambda aln: aln.replace("-",""), aligned_exons) 

    unaligned_dna          = re.sub('Z+', 'Z', unaligned_dna)
    unaligned_exons_dna    = unaligned_dna.split('Z')

    print aligned_exons_stripped
    print unaligned_exons_dna

    if 0:
        for human_exon, exon_names in exon_dict.iteritems():
            for exon_name in exon_names:
                peptide         = alnmt_pep[human_exon][exon_name].replace("-","")
                aligned_peptide = ""
                for ct in range(len(aligned_exons_stripped)):
                    if (aligned_exons_stripped[ct] == peptide):
                        aligned_peptide = aligned_exons[ct]
                        break
                if not aligned_peptide: continue
                dna             = alnmt_dna[human_exon][exon_name].replace("-","")
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


                aligned_dna   = align_nucseq_by_pepseq(aligned_peptide, dna[len_left:len_right])

                if not aligned_dna:
                    print concat_name
                    print peptide
                    print aligned_exons_stripped
                    print aligned_peptide
                    print "dna:", dna
                    print "dna no flanks:",dna[len_left:len_right]
                    print left_flank
                    print right_flank
                    exit(1)

                left_flank  =  left_flank.rjust (flank_length, '-')
                right_flank =  right_flank.ljust(flank_length, '-')
                if output_dna[concat_name]:  output_dna[concat_name] += "-Z-"
                output_dna[concat_name] += left_flank + aligned_dna + right_flank

                print ">> %35s   %5d  %5d " % (concat_name, len(aligned_peptide)*3+20, len( left_flank + aligned_dna + right_flank ))

        for concat_name  in concatenated_exons.keys():
            print "%35s   %5d  %5d " % (concat_name, len(output_pep[concat_name])*3+240, len(output_dna[concat_name]))
          

    exit(1)


    return output_dna


#########################################
def fix_one2many (cfg, acg, sorted_trivial_names, exon_seq_names, concat_name, 
                  list_of_human_exons, seqid, alnmt_pep, output_pep, alnmt_dna, output_dna):
    
    new_alignment_pep = {}
    new_alignment_dna = {}

    exon_numbers = map (lambda ex: seqid[ex], list_of_human_exons)
    smallest_id  = exon_numbers[0]
    largest_id   = exon_numbers[-1]

    if len(list_of_human_exons) == 0: # error
        return output_pep # the same thing that came in 
    elif len(list_of_human_exons) == 1: # many ortho to one human
        many2one = True
        human_exon = list_of_human_exons[0]
        # resolve the possibility of an overlap
        sequence_pieces = []
        sequence_pieces_dna = []
        for exon_seq_name in exon_seq_names:
            sequence_pieces.append( alnmt_pep[human_exon][exon_seq_name])
            sequence_pieces_dna.append(alnmt_dna[human_exon][exon_seq_name])
        [template_name, template_seq] = find_human_template(alnmt_pep[human_exon])
        [seq_to_fix, seq_to_fix_dna]  = merged_sequence (template_seq, sequence_pieces, sequence_pieces_dna)
    else: # one ortho to many human -- the more complicated cases I pretend I do not see
        # sanity
        many2one = False
        prev = ""
        exon_seq_name =  exon_seq_names[0]
        for human_exon in list_of_human_exons:
            current = alnmt_pep[human_exon][exon_seq_name].replace ("-", "")
            if prev and not prev == current: # should be all one and the same
                print "oink? "
                exit (1)
                prev = current
        seq_to_fix = current
        seq_to_fix_dna = alnmt_dna[human_exon][exon_seq_name].replace ("-", "")

    print seq_to_fix
    print seq_to_fix_dna

    # pull  the slice out of the alignment
    # use human as the reference - in other species the boundaries might
    # be at different positions
    # 
    slice_seq = {}
    delimiter = re.compile("Z")

    # slice position in the peptide alignment
    exon_aln_start_pep = []
    exon_aln_end_pep   = []
    start              = 0
    prev_end           = 0
    for match in delimiter.finditer(output_pep['human']):
        start = prev_end 
        end   = match.start()
        exon_aln_start_pep.append(start)
        exon_aln_end_pep.append(end)
        prev_end   =  match.end()

    start = prev_end + 1
    end   = len(output_pep['human'])
    exon_aln_start_pep.append(start)
    exon_aln_end_pep.append(end)

    if len(exon_aln_start_pep) <= smallest_id or  len(exon_aln_start_pep) <= largest_id:
        return [output_pep, output_dna]

    # slice position in the dna alignment
    exon_aln_start_dna = []
    exon_aln_end_dna   = []
    start              =  0
    prev_end           =  0
    for match in delimiter.finditer(output_dna['human']):
        start = prev_end 
        end   = match.start()
        exon_aln_start_dna.append(start)
        exon_aln_end_dna.append(end)
        prev_end   =  match.end()

    start = prev_end + 1
    end   = len(output_dna['human'])
    exon_aln_start_dna.append(start)
    exon_aln_end_dna.append(end)

    if len(exon_aln_start_dna) <= smallest_id or  len(exon_aln_start_dna) <= largest_id:
        return [output_pep, output_dna]

    ####################################
    slice_start = exon_aln_start_pep[smallest_id]
    slice_end   = exon_aln_end_pep  [largest_id]
    for name, seq in output_pep.iteritems():
        if (name== concat_name): continue
        slice_seq[name] = seq[slice_start:slice_end]
        
    exon_seq_name = exon_seq_names[0]
    tmp_name  = exon_seq_name+"_"+"_".join(map (lambda x: str(x), exon_numbers))
    if len(tmp_name) > 50: tmp_name= tmp_name[0:50]
    afa_fnm   = "{0}/{1}.afa".format(cfg.dir_path['scratch'], tmp_name)

    tmp_name = exon_seq_name+"_"+"_".join(map (lambda x: str(x), exon_numbers))
    if len(tmp_name) > 50: tmp_name= tmp_name[0:50]
    out_fnm  = "{0}/{1}.out.afa".format(cfg.dir_path['scratch'], tmp_name)

    if 0: # profile alignment - makes mistakes ...
        output_fasta (afa_fnm, slice_seq.keys(), slice_seq) 

        tmp_name  = exon_seq_name
        if len(tmp_name) > 50: tmp_name= tmp_name[0:50]
        fasta_fnm = "{0}/{1}.fa".format(cfg.dir_path['scratch'], tmp_name)
        output_fasta (fasta_fnm, [concat_name], {concat_name:seq_to_fix})

        mafftcmd = acg.generate_mafft_profile (afa_fnm, fasta_fnm, out_fnm)
        ret      = commands.getoutput(mafftcmd)
        ret      = commands.getoutput("cat "+out_fnm)

    else:
        slice_seq[concat_name] = seq_to_fix
        output_fasta (afa_fnm, slice_seq.keys()+[concat_name], slice_seq) 
        mafftcmd = acg.generate_mafft_command (afa_fnm)
        ret      = commands.getoutput(mafftcmd)
        
    realigned = {}
    pepseq = ""
    for line in ret.split('\n'):
        if '>' in line:
            if ( pepseq ):
                realigned[name] = pepseq
            name = line.replace (">", "")
            name = name.replace (" ", "")
            pepseq = ""
        else:
            pepseq += line
    if pepseq: realigned[name] = pepseq

    # replace the slice with the re-aligned one
    new_alignment_pep = {}
    new_alignment_dna = {}
    dna_slice_start = exon_aln_start_dna[smallest_id]
    dna_slice_end   = exon_aln_end_dna  [largest_id]
    for name, pepseq in output_pep.iteritems():
        if realigned.has_key(name):
            new_alignment_pep[name]  = pepseq[:slice_start] # this should presumably include "Z"
            new_alignment_pep[name] += realigned[name]
            new_alignment_pep[name] += pepseq[slice_end:]
            # cook up dna
            new_alignment_dna[name]  = output_dna[name][:dna_slice_start]
            new_alignment_dna[name] += expand (realigned[name], seq_to_fix_dna)
            new_alignment_dna[name] += output_dna[name][dna_slice_end:]
        else:
            new_alignment_pep[name]  = pepseq
            new_alignment_dna[name]  = output_dna[name]
    
    return [new_alignment_pep, new_alignment_dna]


#########################################
def  sort_concatenated_exons( cursor, ensebl_db_name, concatenated_exons):
    for concat_seq_name in concatenated_exons.keys():
        for human_exon in concatenated_exons[concat_seq_name].keys():
            if len(concatenated_exons[concat_seq_name][human_exon])<=1: continue
            start = {}
            for name in concatenated_exons[concat_seq_name][human_exon]:
                [species, exon_id, exon_known] = parse_aln_name (name)
                ortho_exon = get_exon (cursor, exon_id, exon_known, ensebl_db_name[species])
                #print name, exon_id, ortho_exon.exon_id, ortho_exon.start_in_gene
                start[name] = ortho_exon.start_in_gene
            names_sorted = concatenated_exons[concat_seq_name][human_exon] # this is just an alias now
            names_sorted.sort(key=lambda name: start[name])


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
    #for gene_id in gene_list:
    for gene_id in [412667]: #  wls   
    #for gene_id in [378768]: #  p53


        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)
        print gene_id, stable_id,  get_description (cursor, gene_id)

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
        has_a_map = True
        ct        = 0
        seqid     = {}
        for human_exon in canonical_exons:
            [alnmt_pep[human_exon], alnmt_dna[human_exon]]  = \
                make_exon_alignment(cursor, ensembl_db_name, human_exon, mitochondrial)   
            if not alnmt_pep[human_exon]: 
                has_a_map=False
                break
            seqid[human_exon] = ct
            ct += 1

        # >>>>>>>>>>>>>>>>>>
        if not has_a_map: continue

        # do we have a sequence mapping to multiple human exons?
        ortho_exon_to_human_exon = {}
        for human_exon in canonical_exons:
            for exon_seq_name in alnmt_pep[human_exon].keys():
                if not ortho_exon_to_human_exon.has_key(exon_seq_name):
                    ortho_exon_to_human_exon[exon_seq_name] = [human_exon]
                else:
                    ortho_exon_to_human_exon[exon_seq_name].append(human_exon)

        # >>>>>>>>>>>>>>>>>>
        # find which species we have, and for how many exons
        # we may have two orthologues for the same species
        concatenated_exons = {}
        seq_name = {}
        flags    = {}
        parent_seq_name = {}
        for human_exon in canonical_exons:
            for exon_seq_name, exon_seq in alnmt_pep[human_exon].iteritems():
                (species, exon_id, exon_known) = parse_aln_name(exon_seq_name)
                ortho_gene_id                  = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
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
            for exon_seq_name in alnmt_pep[human_exon].keys():
                concat_seq_name = parent_seq_name[exon_seq_name]
                if not concatenated_exons.has_key(concat_seq_name): 
                    concatenated_exons[concat_seq_name] = {}
                if not concatenated_exons[concat_seq_name].has_key(human_exon): 
                    concatenated_exons[concat_seq_name][human_exon] = []
                concatenated_exons[concat_seq_name][human_exon].append(exon_seq_name)

        # make sure the exons that concatenated_exons refers to are sorted according to the order 
        # in which they appear in the gene
        sort_concatenated_exons( cursor, ensembl_db_name, concatenated_exons)
        # >>>>>>>>>>>>>>>>>>
        # concatenate the aligned exons for each species, taking into account that the alignment
        # doesn't have to be one to one
        headers     = []
        output_pep  = {}
        output_dna  = {}
        for concat_seq_name in concatenated_exons.keys():
            output_pep[concat_seq_name] = ""
            output_dna[concat_seq_name] = ""

            flagged_exons = []
            if concat_seq_name in flags.keys() and flags[concat_seq_name].has_key('one_ortho2many_human'):
                for to_fix in flags[concat_seq_name]['one_ortho2many_human']:
                    flagged_exons = flagged_exons+to_fix[1]

            for human_exon in canonical_exons:

                aln_length = len(alnmt_pep[human_exon].itervalues().next())
                if human_exon in flagged_exons:
                    # one ortho seq maps to multiple human exons
                    pep = '-'*aln_length
                    dna = '-'*(3*aln_length+2*flank_length)
                elif concatenated_exons[concat_seq_name].has_key(human_exon):
                    if (len(concatenated_exons[concat_seq_name][human_exon]) == 1):
                        # we have a neat one-to-one mapping
                        exon_seq_name = concatenated_exons[concat_seq_name][human_exon][0]
                        pep = alnmt_pep[human_exon][exon_seq_name]
                        dna = alnmt_dna[human_exon][exon_seq_name]
                    else:
                        # if two sequences map to the same human exon, merge
                        pep = '-'*aln_length
                        dna = '-'*(3*aln_length+2*flank_length)
                        # flag for later
                        flags_init(flags, concat_seq_name, 'many_ortho2one_human')
                        to_fix = [human_exon, concatenated_exons[concat_seq_name][human_exon]]
                        flags[concat_seq_name]['many_ortho2one_human'].append(to_fix)
                else: 
                    # no exon in this species
                    pep =  '-'*aln_length
                    dna = '-'*(3*aln_length+2*flank_length)
                   
                if output_pep[concat_seq_name]: output_pep[concat_seq_name] += '-Z-'
                if output_dna[concat_seq_name]: output_dna[concat_seq_name] += '--Z--'
                output_pep[concat_seq_name] += pep
                output_dna[concat_seq_name] += dna
                   
            headers.append(concat_seq_name)

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
                        [output_pep, output_dna]  = fix_one2many (cfg, acg, sorted_trivial_names, [exon_seq_name], seq_to_fix, 
                                                   list_of_human_exons, seqid, alnmt_pep, output_pep, alnmt_dna, output_dna)
        if 1:
            for seq_to_fix in flags.keys():
                 if flags[seq_to_fix].has_key('many_ortho2one_human'): 
                    for to_fix in flags[seq_to_fix]['many_ortho2one_human']:
                        human_exon          = to_fix[0]
                        list_of_ortho_exon_seq_names = to_fix[1]
                        [output_pep, output_dna] = fix_one2many (cfg, acg, sorted_trivial_names, list_of_ortho_exon_seq_names, 
                                                   seq_to_fix, [human_exon], seqid, alnmt_pep, output_pep, alnmt_dna, output_dna)
        # cleanup
        cleanup_exon_boundaries(output_pep)
        output_pep = strip_gaps(output_pep)
        output_dna = strip_gaps(output_dna)

        # >>>>>>>>>>>>>>>>>>
        # find place for best_afa -- for the moment can put it to the scratch space:
        afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        # afa_fnm  = 'test.afa'
        sorted_seq_names = sort_names (sorted_trivial_names['human'], output_pep)

        output_fasta (afa_fnm, sorted_seq_names, output_pep)
        print afa_fnm
        afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        output_fasta (afa_fnm, sorted_seq_names, output_dna)
        print afa_fnm

        # continue
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
        # see if sorting the species will help any -- it won't
        #if concat_name[-1].isdigit():
        #    aux = concat_name.split("_")
        #    base_name = "_".join(aux[:-1])
        #else:
        #    base_name = concat_name
        #sorted_names   = sort_names (sorted_trivial_names[base_name], slice_seq)
  '''
