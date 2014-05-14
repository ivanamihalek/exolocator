#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import pdb
#pdb.set_trace()

import MySQLdb, commands, re, os, time
import random, string
import inspect
from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.utils   import  *
from el_utils.el_specific   import  *
from el_utils.map     import  Map, get_maps, map2exon
from el_utils.tree    import  species_sort
from el_utils.ncbi    import  *
from el_utils.almt_cmd_generator import AlignmentCommandGenerator
from el_utils.config_reader      import ConfigurationReader
from el_utils.translation        import phase2offset, translation_bounds, crop_dna, translate
from el_utils.special_gene_sets  import *
from el_utils.processes import parallelize
from el_utils.exon_boundary_hacks import *
from bitstring import Bits
from alignment import * # C implementation of smith waterman
from   random  import choice
# BioPython
from Bio          import  SeqIO
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

import pdb

verbose = True


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
def sort_exon_names (list_of_exon_names):

    start_in_gene = {}
    for exon_seq_name in list_of_exon_names:
        [species, exon_id, exon_known, exon_start] = parse_aln_name(exon_seq_name)
        start_in_gene[exon_seq_name]  = exon_start
    list_of_exon_names.sort(key=lambda en: start_in_gene[en])

    return 
    
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

            exon_seq_name_2   = exon_seq_names[j]  
            human_exons_2     = ortho_exon_to_human_exon[exon_seq_name_2]

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
            ortho_exons = list(set(join_group) ) # this should take care of duplicates
            overlapping_maps.append([human_exons, ortho_exons])



    return overlapping_maps

#########################################
def remove_ghosts (output_pep, sequence_name_to_exon_names):

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
            for exon_ct in range(len(sequence_name_to_exon_names[name])):
                if not exon_ct in removed_exons:
                    new_names.append(sequence_name_to_exon_names[name][exon_ct])
            sequence_name_to_exon_names[name] = new_names

    # strip gaps
    output_pep = strip_gaps(output_pep)
    if not output_pep:  
        c=inspect.currentframe()
        print " in %s:%d" % (c.f_code.co_filename, c.f_lineno)
        return [None, None]
            

    return [output_pep, sequence_name_to_exon_names]



#########################################
def check_notes_directory (cfg):
    
    directory = "{0}/notes".format(cfg.dir_path['afs_dumps'])
    if not os.path.exists(directory):
        try:
            os.makedirs(directory) 
        except:
            print "error making", directory
            exit(1) # exit after an  error making the 'notes' directory

    return directory

#########################################
def find_maps_to (cursor, ensembl_db_name,  human_exon_to_ortho_exon, concat_seq_name, exon_seq_name):
    
    stable_ids_str = ""
    stable_ids     = []


    switch_to_db ( cursor, ensembl_db_name['homo_sapiens'])

    for human_exon in human_exon_to_ortho_exon[concat_seq_name].keys():
        if not exon_seq_name in human_exon_to_ortho_exon[concat_seq_name][human_exon]: continue
        stable_id = exon2stable (cursor, human_exon.exon_id)
        stable_ids.append(stable_id)


    stable_ids_str = ";".join(stable_ids)

    return stable_ids_str


#########################################
#########################################
def print_notes (cursor, cfg,  ensembl_db_name, output_pep, sequence_to_exons, sorted_seq_names, 
                 human_stable_id, human_exon_to_ortho_exon, assorted_notes):

    gene_id        = {}
    stable_gene_id = {}
    sci_name       = {}
    for name in output_pep.keys():
        for exon_name in sequence_to_exons[name]:
            [species, exon_id, exon_known, exon_start] = parse_aln_name(exon_name)
            ortho_gene_id                  = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
            # sequence can contain pieces from "different" genes - if they are from different scaffold pieces
            # or from not-so-distant pieces on a chromosome
            if not gene_id.has_key(name):  
                gene_id[name]        = []
                stable_gene_id[name] = []

            if not ortho_gene_id in gene_id[name]:
                gene_id[name].append(ortho_gene_id)
                stable_gene_id[name].append(gene2stable(cursor, ortho_gene_id, ensembl_db_name[species]))
                sci_name[name]       = species

    descr = get_description (cursor, gene_id['human'][0], ensembl_db_name['homo_sapiens'])

    # write to string
    out_string  = "% Notes to accompany the alignment of (tentative) orthologues\n"
    out_string += "%% for the canonical transcript of the human gene %s, \n" % human_stable_id
    out_string += "%% %s\n"  % descr
    out_string += "%%\n"
    out_string += "% The alignment shows the exons from the following genes: \n" 
    out_string += "%% %-30s  %-30s  %-30s \n" % ('species', 'common_name', 'gene_id')

    sorted_seq_names = filter (lambda name: name in output_pep.keys(), sorted_seq_names)
    for name in sorted_seq_names:
        if not sci_name.has_key(name) or not stable_gene_id.has_key(name): continue
        out_string += " %-30s  %-30s  %s" % ( name, sci_name[name],  stable_gene_id[name][0])
        for stable_id in stable_gene_id[name][1:]:
            out_string += ", %s" % stable_id
        out_string += "\n"
    
    if assorted_notes:
        out_string += assorted_notes

    out_string += "\n" 
    out_string += "% The following exons appear in the alignment\n" 
    out_string += "% Note: the exons assigned to the same peptide sequence might belong to several \"genes\"\n"
    out_string += "% If the yare found split across several scaffolds\n"
    out_string += "% (in which case Ensembl assigned two different identifiers to the two exons sets).\n"

    novel       = []
    novel_annot = {}
 
    for name in sorted_seq_names:
        if not stable_gene_id.has_key(name): 
            #print name, "bleep"
            continue
        out_string += "\n" 
        out_string += "%% sequence name: %s   corresponding to the gene: %s" % (name, stable_gene_id[name][0])
        for stable_id in stable_gene_id[name][1:]:
            out_string += ", %s" % stable_id
        out_string += "\n"

        out_string += "%% %50s  %10s  %10s  %6s  %6s    %-s  %-s\n" % \
           ('exon_id', 'gene_from', 'gene_to', 'coding', 'canon', 'source', 'maps_to_human_exon')

        for exon_name in sequence_to_exons[name]:
            [species, exon_id, exon_known, exon_start] = parse_aln_name(exon_name)
            exon = get_exon (cursor, exon_id, exon_known, ensembl_db_name[species])
            if exon_known == 1:
                exon_stable_id = exon2stable(cursor, exon_id, ensembl_db_name[species])
            else:
                exon_stable_id = exon_name

            if (exon_known ==2 or exon_known ==3):

                source = "SW# " if (exon_known ==2) else "usearch"
                
                template_species     = exon.template_species
                template_exon_seq_id = exon.template_exon_seq_id
                [template_exon_id, template_exon_known] = exon_seq_id2exon_id (cursor, 
                                                        template_exon_seq_id, 
                                                        ensembl_db_name[template_species])
                if (template_exon_known==1):
                    template_stable = exon2stable(cursor, int(template_exon_id), ensembl_db_name[template_species])
                else:
                    template_stable = 'novel'

                novel.append(exon_stable_id)
                if exon.has_NNN is None:
                     unseq = '-'
                elif exon.has_NNN: 
                     unseq = 'Y'
                else:
                     unseq = 'N'

                if exon.has_stop is None:
                     has_stop = '-'
                elif exon.has_stop: 
                     has_stop = 'Y'
                else:
                     has_stop = 'N'

                if exon.has_5p_ss is None:
                    fivep_ss = '-'
                else:
                    fivep_ss = exon.has_5p_ss

                if exon.has_3p_ss is None:
                    threep_ss = '-'
                else:
                    threep_ss = exon.has_3p_ss

                novel_annot[exon_stable_id] =  [exon_stable_id, source, template_stable, template_species, unseq, 
                                                has_stop, threep_ss, fivep_ss]
            else:
                source = get_logic_name (cursor, exon.analysis_id,  ensembl_db_name[species])
                
            maps_to_human_stable = find_maps_to (cursor, ensembl_db_name, human_exon_to_ortho_exon, name, exon_name)
                
            out_string += "  %50s   %10s  %10s     %-6s  %-10s   %-s   %-s \n" % \
                (exon_stable_id,  exon.start_in_gene, exon.end_in_gene,
                 exon.is_coding, exon.is_canonical, source, maps_to_human_stable)
               
    if novel:
        out_string += "\n" 
        out_string += "%% exons found by the exolocator pipeline: \n" 
        out_string += "%%    unseq_regions: the  region contains NNNN stretches \n" 
        out_string += "%%    me_score:  MaxEntScan score for the intron splice signal\n" 
        out_string += "%%    3pss: 3' splice signal check\n" 
        out_string += "%%    5pss: 5' splice signal check\n" 
        out_string += "%%    (note: by convention 5' and 3' refer to the intervening introns, not the exons;\n" 
        out_string += "%%     thus 3' check refers to checking the region flanking the exon's 5' end)\n" 

        out_string += "%% %50s  %10s  %15s  %30s    %-10s  %-10s   %-20s  %-20s\n" % \
           ('name   ', 'source   ',' template_stable', 'template_species', 
            'unseq_regions', 'has_stop', '3ss', '5pss' )
        for exon_name in novel:
            out_string += " %50s  %10s  %15s  %30s   %-10s  %-10s   %-s  %-s\n" %  tuple(novel_annot[exon_name])
    else:
        print 'no novel exons ound by the exolocator pipeline'

    directory = check_notes_directory (cfg)
    notes_fnm = directory + '/'+human_stable_id+'.txt'
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
def name2count (output_pep, sequence_to_exons):

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
        if not nonempty_exon_ct == len( sequence_to_exons[name]):
            print " foul:    nonempty_exon_ct={0}".format(nonempty_exon_ct),
            print " len(sequence_to_exons[name])={0}".format(len(sequence_to_exons[name])),
            print name
            return [{},{}]

    return [name_ct2exon_ct, exon_ct2name_ct]

#########################################
def exon_seq_check(exon_seqs, pep_aligned, species, exon_name):

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
def expand_protein_to_dna_alnmt (cursor, ensembl_db_name, cfg, acg, sorted_seq_names, 
                                 sequence_to_exons, alnmt_pep, output_pep, flank_length):

    output_dna         = {}
    # which exons correspond to which name?
    [name_ct2exon_ct, exon_ct2name_ct] = name2count (output_pep, sequence_to_exons)
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

        prev_local_pos   = -3
        exon_ct          = -1
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
                exon_seq_name = sequence_to_exons[name][name_ct]
                [species, exon_id, exon_known, exon_start] = parse_aln_name(exon_seq_name)
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
def realign_slice (pep_slice, seq_to_fix, pep_seq_pieces):

    new_pep_slice  = pep_slice

    human_seq      =  decorate_and_concatenate(pep_slice['human'].split('Z'))
    pep_seq_to_fix =  decorate_and_concatenate(pep_seq_pieces)


    [aligned_human, aligned_ortho] \
        = smith_waterman_context (human_seq, pep_seq_to_fix, -3, 0)


    aligned_ortho = check_Z_right(re.sub('[\dB\#]', '-', aligned_ortho))
    aligned_human = check_Z_right(re.sub('[\dB]'  , '#', aligned_human), pound=True)

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

#########################################
def check_overlap (alnmt_length, seqs):

    overlap = []

    for pos in range(alnmt_length): 
        for i in range(len(seqs)):
            if (seqs[i][pos] == '-'): continue
            for j in range (i+1, len(seqs)):
                if (seqs[j][pos] == '-'): continue
                index = str(i) + " " + str(j)
                if not index in overlap:
                    overlap.append(index)
                break
    return overlap

#########################################
def check_seq_overlap (cursor, ensembl_db_name, cfg, acg, template_seq, pep_seq_pieces, pep_seq_names, sequence_to_exons):
    
    #I should resolve this at some earlier place - like when I am writing the exon maps
    seq_names_to_remove = []

    template_length = len(template_seq)

    for piece in pep_seq_pieces:
        if (not len(piece) == template_length):
            #print "length mismatch for aligned exons (?)"
            #print piece
            #print template_seq
            return []

    # check whether any two pieces overlap
    overlap = check_overlap (template_length, pep_seq_pieces)
    if overlap: 
        # if there is an overlap, check that it is not the artefact of the multiple sequence alignment
        # (sometimes it clears up when only two seqs are aligned)
        randstr = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(10))
        fasta_fnm = "{0}/{1}.fa".format( cfg.dir_path['scratch'], randstr)

        sequences = {'template':template_seq.replace('-','')}
        #concatenate pieces by force - we trust here they are  given in the order in which they appear in the gene
        sequences['other'] = 'Z'.join( [pepseq.replace('-','') for pepseq in  pep_seq_pieces] )    
        output_fasta (fasta_fnm, sequences.keys(), sequences)

        # align
        afa_fnm  = "{0}/{1}.afa".format( cfg.dir_path['scratch'], randstr)
        mafftcmd = acg.generate_mafft_command (fasta_fnm, afa_fnm)
        ret      = commands.getoutput(mafftcmd)

        new_pep_seq_pieces = []
        inf = erropen(afa_fnm, "r")
        for record in SeqIO.parse(inf, "fasta"):
            if record.id == 'other':
                other_seq = record.seq
            else:
                new_template_seq = record.seq
        inf.close()
        commands.getoutput("rm "+afa_fnm+" "+fasta_fnm)

        prev = 0
        new_pep_seq_pieces = []
        template_pieces    = []
        for i in range(len(other_seq)):
            if other_seq[i] == 'Z' and i>0 :
                new_pep_seq_pieces.append(other_seq[prev:i])
                template_pieces.append(new_template_seq[prev:i])
                prev = i+1

        if prev < len(other_seq):
            new_pep_seq_pieces.append(other_seq[prev:len(other_seq)])
            template_pieces.append(new_template_seq[prev:len(other_seq)])
            
        # check the similarity of the obtained pieces
        min_similarity = cfg.get_value('min_accptbl_exon_sim')
        for i in range(len(new_pep_seq_pieces)):
            if ( pairwise_tanimoto (template_pieces[i], new_pep_seq_pieces[i]) < min_similarity ):
                seq_names_to_remove.append(pep_seq_names[i])
        
        new_sequence_to_exons = filter (lambda exon: exon not in seq_names_to_remove, sequence_to_exons)


    else:
        new_pep_seq_pieces     = pep_seq_pieces
        new_sequence_to_exons = sequence_to_exons

    # make sure 
             
    return new_sequence_to_exons

########################################
def fix_one2many (cursor, ensembl_db_name, cfg, acg, sorted_seq_names, canonical_human_exons, human_exon_to_ortho_exon,
                  sequence_name_to_exon_names, seq_to_fix, overlapping_maps, alnmt_pep, output_pep):

    count = 0

    #pdb.set_trace()
    # no overlapping maps - nothing to resolve"
    if not overlapping_maps: return [output_pep, sequence_name_to_exon_names]
    # not sure how this could happen
    if not output_pep.has_key('human'): return [output_pep, sequence_name_to_exon_names]

    # get rid of things that are duplicates or have not been associated with any humam exon
    list_of_ok_exon_names = []
    for ortho_exon in sequence_name_to_exon_names[seq_to_fix]:
        for human_exon in canonical_human_exons:
            if not alnmt_pep[human_exon].has_key(ortho_exon): continue
            if ortho_exon in list_of_ok_exon_names: continue
            list_of_ok_exon_names.append(ortho_exon)
     
    # find sequential numbers of exons that we have in this story
    seqid = {}
    ct    = 0
    for human_exon in canonical_human_exons:
        seqid[human_exon] = ct
        ct += 1
    number_of_human_exons = ct

    current_pep = output_pep
    
    if seq_to_fix == 'xenopus':
        print output_pep['human']
        for [human_exons, ortho_exons] in overlapping_maps:
            print ">>",
            for he in human_exons:
                print he.exon_id,
            print ortho_exons

    # for each unresolved "map"  cut out the slice and re-align
    for  [human_exons, ortho_exons] in overlapping_maps:

        exon_numbers = sorted(map(lambda ex: seqid[ex], human_exons))
        smallest_id  = exon_numbers[0]
        largest_id   = exon_numbers[-1]


        # check sequence overlap, if several map to the same  human exon
        for human_exon in human_exons:
            [template_name, template_seq]  = find_human_template(alnmt_pep[human_exon])
            sequence_pieces = []
            sequence_piece_names = []
            for exon_seq_name in ortho_exons:
                if not alnmt_pep[human_exon].has_key(exon_seq_name):  continue
                sequence_pieces.append(alnmt_pep[human_exon][exon_seq_name])
                sequence_piece_names.append(exon_seq_name)
            list_of_ok_exon_names = check_seq_overlap(cursor, ensembl_db_name, cfg, acg, template_seq, 
                                                      sequence_pieces, sequence_piece_names, list_of_ok_exon_names)
            if seq_to_fix=="xenopus":
                print 'ooooooooooooooooooooooooooooooooooooooooooo'
                print list_of_ok_exon_names
                switch_to_db(cursor, ensembl_db_name['homo_sapiens'])
                exon_seqs = get_exon_seqs (cursor, human_exon.exon_id, 1)
                [exon_pep_seq, trsl_from, trsl_to, exon_left_flank,
                 exon_right_flank, exon_dna_seq] = exon_seqs [1:]
                print "exon:", human_exon.exon_id, "covering exon:", human_exon.covering_exon,  "pepseq:", exon_pep_seq
 
        # join sequences that are deemed to be ok
        pep_seq_pieces = [] 
        for ortho_exon in list_of_ok_exon_names:
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
            print "abort 1"
            return [output_pep, sequence_name_to_exon_names]

        ####################################
        # find the slice position
        pep_slice_start = exon_aln_start_pep[smallest_id]
        pep_slice_end   = exon_aln_end_pep  [largest_id]

        ####################################
        # if the slice size becomes comparable to the alignment size, drop the sequence
        if  False and number_of_human_exons > 5 and \
                ( float(pep_slice_end-pep_slice_start)/len(output_pep['human']) > 0.3 or 'elephant' in seq_to_fix):
            #print "deleting", seq_to_fix
            print "abort 2: map apparently spread over too many human exons:"
            print seq_to_fix, " pep slice", pep_slice_start, pep_slice_end, 
            print "  len: ", len(output_pep['human']), "no human exons:", number_of_human_exons
            del output_pep[seq_to_fix]
            return [output_pep, sequence_name_to_exon_names]


        ####################################
        # cut out the slice
        pep_slice = {}
        for name in current_pep.keys():            
            pep_slice[name] = current_pep[name][pep_slice_start:pep_slice_end]

        # slice realign
        new_pep_slice = realign_slice (pep_slice, seq_to_fix, pep_seq_pieces)

        output_fasta ("old_slice.afa", pep_slice.keys(), pep_slice);
        output_fasta ("new_slice.afa", new_pep_slice.keys(), new_pep_slice);

        # strip gaps and output
        # boundary_cleanup(new_pep_slice, new_pep_slice.keys())

        if not check_seq_length(new_pep_slice, "new_pep_slice"):
            del output_pep[seq_to_fix]
            del sequence_name_to_exon_names[seq_to_fix]
            del human_exon_to_ortho_exon[seq_to_fix]
            print "abort 3"
            c=inspect.currentframe()
            print " in %s:%d" % (c.f_code.co_filename, c.f_lineno)
            return [output_pep, sequence_name_to_exon_names] # ie return as it was

        #################################### 
        # replace the slice with the re-aligned one
        new_alignment_pep = {}
        for name, pepseq in current_pep.iteritems():
            new_alignment_pep[name]  = pepseq[:pep_slice_start]
            new_alignment_pep[name] += new_pep_slice[name]
            new_alignment_pep[name] += pepseq[pep_slice_end:] 

        if not check_seq_length(new_alignment_pep, "new_alignment_pep"):
            del output_pep[seq_to_fix]
            del sequence_name_to_exon_names[seq_to_fix]
            del human_exon_to_ortho_exon[seq_to_fix]
            print "abort 4"
            c=inspect.currentframe()
            print " in %s:%d" % (c.f_code.co_filename, c.f_lineno)
            return [output_pep, sequence_name_to_exon_names]

        boundary_cleanup(new_alignment_pep, new_alignment_pep.keys())
        new_alignment_pep = strip_gaps(new_alignment_pep)
        if not new_alignment_pep:  
            print "abort 5"
            c=inspect.currentframe()
            print " in %s:%d" % (c.f_code.co_filename, c.f_lineno)
            return [output_pep, sequence_name_to_exon_names]
        
        current_pep = new_alignment_pep


    pep_exons = new_alignment_pep[seq_to_fix].split ('Z') 
    exon_ct   = 0
    for pe in pep_exons:
        pe = pe.replace('-','')
        if pe:   exon_ct += 1

    if not exon_ct == len(list_of_ok_exon_names):
        print "///////////////////////////////////////////////////"
        c=inspect.currentframe()
        print "in %s:%d" % (c.f_code.co_filename, c.f_lineno), ", after religning slice"
        print seq_to_fix + ' :',
        print "exon_ct (%d) not equal to the length of list_of_ok_exon_names (%d)" % ( exon_ct, len(list_of_ok_exon_names))
        print "ok exon names: ", list_of_ok_exon_names
        print "non-zero petide sequences (for exons)"
        print "\n".join( map (lambda seq: seq.replace('-','') + " *** ", pep_exons) )
        print "==================================================="
        return [output_pep, sequence_name_to_exon_names] # theses are empty

    output_pep = new_alignment_pep
    sequence_name_to_exon_names[seq_to_fix] = list_of_ok_exon_names

    return [output_pep, sequence_name_to_exon_names]



#########################################
def fix_split_codons (cursor, ensembl_db_name, cfg, acg, sorted_seq_names, 
                      mitochondrial, sequence_name_to_exon_names, alnmt_pep, output_pep, flank_length):

    ##########
    output_pep_new  = {}

    ##########
    if not 'human' in output_pep.keys(): # this shouldn't really happen ...
        return output_pep_new

    # which exons correspond to which name?
    [name_ct2exon_ct, exon_ct2name_ct] = name2count (output_pep, sequence_name_to_exon_names)
    if not name_ct2exon_ct: 
        return output_pep_new

    # in the peptide alignment find exon positions that are  global exon boundaries, 
    # and the ones that are local
    [global_bdry_position, local_bdry_position] = find_exon_boundaries(output_pep)


    output_pep_new = {}

    # "insert" here will be the codon that we will always glue to the right of Z denoting the boundary
    insert_global  = []
    insert         = {}
    for name  in output_pep.keys():
       insert[name] = {}

    for pos in global_bdry_position:
        for name  in output_pep.keys():
            insert[name][pos] = False

    for name, full_aligned_pepseq in output_pep.iteritems():
 
        prev_local_pos   = -3
        exon_ct          = -1
        prev_right_flank = ""
        for local_pos in local_bdry_position[name]+[len(full_aligned_pepseq)]:

            
            human_seq      = output_pep['human'][prev_local_pos+3:local_pos]
            pep_seq        = full_aligned_pepseq[prev_local_pos+3:local_pos]

            exon_ct    += 1
            name_ct     = exon_ct2name_ct[name][exon_ct]
            if (name_ct < 0):
                prev_right_flank = None
                continue

            exon_seq_name = sequence_name_to_exon_names[name][name_ct]
            [species, exon_id, exon_known, exon_start] = parse_aln_name(exon_seq_name)
            exon_seqs     = get_exon_seqs(cursor, exon_id, exon_known, ensembl_db_name[species])[1:]

            if not exon_seqs:
                prev_right_flank = None
                continue

            [exon_pep_seq, trsl_from, trsl_to, exon_left_flank,
             exon_right_flank, exon_dna_seq] = exon_seqs  


            phase = get_exon_phase (cursor, exon_id, exon_known)


            if phase > 0 and prev_right_flank:

                offset    = (3-phase)%3
                cary      = prev_right_flank[:phase]
                tmp_patch = exon_left_flank.lower() + exon_dna_seq[:trsl_from]
                flanking_nucleotides   = tmp_patch [-offset:]
                codon                  = cary + flanking_nucleotides
                [phase_suggested, res] = translate (codon, 0, mitochondrial, strip_stop = False)

                if res: # we are going to have an insert here
                    insert [name][prev_local_pos+3] = res
                    if not (prev_local_pos+3) in insert_global:
                        insert_global.append(prev_local_pos+3)

            # so here are the ad-hoc rules we will use: we will say that we have a valid previous flank
            # 1) the previous exon maps to previous human exon
            #        we are mkaing sure that this is the case by setting prev_right_flank to None
            #        if there is no mapping
            # 2) the end of the previous exon agrees with the end of the previous human exon
            # now check that the right ends of the exon match:
            end_match = False
            for i in range(len(human_seq)):
                h = human_seq[-i]
                p = pep_seq[-i]
                if h == '-' and p == '-': continue
                if not h == '-' and not p == '-':
                    end_match = True
                break
            if end_match:
                prev_right_flank  = exon_dna_seq [trsl_to:]
                prev_right_flank += exon_right_flank.lower()
            else:
                pass

            prev_local_pos = local_pos


    #######################################################################
    #rerun one more time - put in the inserts    
    insert_global.append(len(full_aligned_pepseq))
    insert_global.sort()


    for name, full_aligned_pepseq in output_pep.iteritems():

        pos                  = insert_global[0]
        pep_seq              = full_aligned_pepseq[:pos]
        output_pep_new[name] = pep_seq
        prev_pos             = pos

        for  pos in insert_global[1:]:
            pep_seq        = full_aligned_pepseq[prev_pos:pos]

            if prev_pos in insert[name].keys() and insert[name][prev_pos]:
                new_pep = ''
                for i in range(len(pep_seq)):
                    if not pep_seq[i] == '-':  break
                    new_pep += '-'
                if (len(new_pep) == len(pep_seq)) : # ie if these are gaps only 
                    # not sure how this could happen
                    new_pep = new_pep[:-1]
                new_pep += insert[name][prev_pos]
                new_pep += pep_seq [i:]
                output_pep_new[name] += new_pep
            
            else:
                output_pep_new[name] += "-"
                output_pep_new[name] += pep_seq

            prev_pos = pos

    return output_pep_new
   


#########################################
def find_pairs (some_list):
    list_length = len(some_list)
    pairs = []
    for i in range(list_length):
        for j in range(i+1,list_length):
            pairs.append([some_list[i], some_list[j]])
    return pairs

#########################################
def find_lower_denom_name(para1, para2):
    [name_to_keep, name_to_drop] = [para1, para2]

    name_pieces = para1.split("_")
    if not isinteger(name_pieces[-1]): 
        name_to_keep = para1
        name_to_drop = para2
    else:
        id1 = int(name_pieces[-1])
        name_pieces = para2.split("_")
        if not isinteger(name_pieces[-1]): 
            name_to_keep = para2
            name_to_drop = para1
        else:
           id2 = int(name_pieces[-1]) 
           if (id1 < id2):
               name_to_keep = para1
               name_to_drop = para2
           else:
               name_to_keep = para2
               name_to_drop = para1

    return [name_to_keep, name_to_drop]

#########################################
#########################################
def delete (name_to_drop, output_pep, sequence_name_to_exon_names, human_exon_to_ortho_exon, deleted):
    deleted.append(name_to_drop)
    del output_pep[name_to_drop]
    del sequence_name_to_exon_names[name_to_drop]
    del human_exon_to_ortho_exon[name_to_drop]


#########################################
def fuse_seqs_split_on_scaffolds (cursor, acg,  ensembl_db_name, output_pep, sequence_name_to_exon_names, 
                                  ortho_exon_to_human_exon, canonical_human_exons, human_exon_to_ortho_exon):

    notes = ""

    # find species that have multiple orthologues
    
    mulitple_orthos = []
    if ( type(output_pep) is str ):
        print "output_pep:", output_pep
        return ""
    for seq_name in output_pep.keys():
        name_pieces = seq_name.split("_")
        if not isinteger(name_pieces[-1]): continue
        species_common = "_".join(name_pieces[:-1])
        if species_common not in mulitple_orthos: mulitple_orthos.append(species_common)
    
    for species_common in mulitple_orthos:
        paralogues = filter (lambda seq_name: species_common in seq_name,  output_pep.keys())
        # for all pairs of paralogues
        deleted = []
        for [para1, para2] in find_pairs (paralogues):
            if para1 in deleted or para2 in deleted: continue

            # do they map to disjoint set of human exons?
            human_exons_1 = set([])
            for exon_name in sequence_name_to_exon_names[para1]:
                human_counterparts = ortho_exon_to_human_exon[exon_name]
                if human_counterparts: 
                    human_exons_1.update(set(human_counterparts))

            if not human_exons_1: # again, no idea how this happens
                delete (para1, output_pep, sequence_name_to_exon_names, human_exon_to_ortho_exon, deleted)
                continue
           
            human_exons_2 = set([])
            for exon_name in sequence_name_to_exon_names[para2]:
                human_counterparts = ortho_exon_to_human_exon[exon_name]
                if human_counterparts: 
                    human_exons_2.update(set(human_counterparts))

            if not human_exons_2: # again, no idea how this happens
                delete (para2, output_pep, sequence_name_to_exon_names, human_exon_to_ortho_exon, deleted)
                continue

            if  human_exons_1 & human_exons_2:continue # there is intersection - we move on

            # do these seqs belong to different pieces of sequence?
            # if one is on the left from an exon on the other group, then they should be all

            he1 = iter(human_exons_1).next()
            he2 = iter(human_exons_2).next()
            interspersed = False
            if he1.end_in_gene < he2.start_in_gene:
                for other_he1 in human_exons_1:
                    if interspersed: break
                    for he2 in human_exons_2:
                        if not other_he1.end_in_gene < he2.start_in_gene: 
                            interspersed = True
                            break
            elif  he1.start_in_gene > he2.end_in_gene:
                 for other_he1 in human_exons_1:
                    if interspersed: break
                    for he2 in human_exons_2:
                        if not other_he1.start_in_gene > he2.end_in_gene: 
                            interspersed = True
                            break
            else:
                interspersed = True
               
            if interspersed: continue         
            
            # one last check, the most expensive if I'm not mistaken:
            # are these two 'genes'  on different scaffolds?
            # or perhaps we might still take them to actually represent a single gene
            # if they are on the same chromosome, not too far apart
            # we definitely do not want them to be on different chromosomes
            # use get_gene_seq function; note that file_names might be several names separated by space
            # [gene_seq, canonical_exon_pepseq, file_names] = get_gene_seq(acg, cursor, gene_id, species)
            
            # at this point, para1 and para2 should belong to the same 'gene'
            [exon_id, exon_known] = sequence_name_to_exon_names[para1][0].split ("_")[-3:-1]
            species    = "_".join (sequence_name_to_exon_names[para1][0].split ("_")[:-3]) 
            gene_id_1   = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
            # I have no idea how this could happen, but it does:
            if not gene_id_1: continue
            stable_id_1 = gene2stable(cursor, gene_id_1, ensembl_db_name[species])
            [gene_seq, canonical_exon_pepseq, file_name_1, seq_name_1, start_1, end_1] = \
                get_gene_seq(acg, cursor, gene_id_1, species)
 
            [exon_id, exon_known] = sequence_name_to_exon_names[para2][0].split ("_")[-3:-1]
            gene_id_2   = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
            # I have no idea how this could happen, but it does:
            if not gene_id_2: continue
            stable_id_2 = gene2stable(cursor, gene_id_2, ensembl_db_name[species])
            [gene_seq, canonical_exon_pepseq,  file_name_2, seq_name_2, start_2, end_2] = \
                get_gene_seq(acg, cursor, gene_id_2, species)
            
            if file_name_1 == file_name_2 and seq_name_1 == seq_name_2:
                    continue
                # I should also check here how far apart are these two seqs
                # if they are two far apart, the assumption that this is actually the same gene
                # might be valid; but then the qiestion is how far apart is too far apart
                # (and then how many cases like that we have)
            else:
                pass
  
            # if we got so far, join the two seqs under the lower denominator name
            [name_to_keep, name_to_drop] = find_lower_denom_name(para1, para2)

            print "keep:", name_to_keep, "  drop:", name_to_drop
            new_seq = ""
            for pos in range(len(output_pep[para1])):
                if ( output_pep[para1][pos] == '-'):
                    new_seq += output_pep[para2][pos] 
                else:
                    new_seq += output_pep[para1][pos] 
            output_pep[name_to_keep] = new_seq
 
 
            # make sure we have all exons listed correctly under the new name - we'll need them to expand dna
            new_exon_set = []
            new_map      = {}
            for human_exon in canonical_human_exons:
                if human_exon in human_exon_to_ortho_exon[para1].keys():
                    for ex1 in  sequence_name_to_exon_names[para1]:
                        if ex1 in human_exon_to_ortho_exon[para1][human_exon]:
                            new_exon_set.append(ex1)
                            if not new_map.has_key(human_exon):  new_map[human_exon] = []
                            new_map[human_exon].append(ex1)
                
                if human_exon in human_exon_to_ortho_exon[para2].keys():
                    for ex2  in  sequence_name_to_exon_names[para2]:
                        if ex2 in human_exon_to_ortho_exon[para2][human_exon]:
                            new_exon_set.append(ex2)
                            if not new_map.has_key(human_exon):  new_map[human_exon] = []
                            new_map[human_exon].append(ex2)

            sequence_name_to_exon_names[name_to_keep] = new_exon_set
            human_exon_to_ortho_exon[name_to_keep] = new_map
            delete (name_to_drop, output_pep, sequence_name_to_exon_names, human_exon_to_ortho_exon, deleted)
            
            notes += "{0}  {1}:{2},{3}-{4}   {5}:{6},{7}-{8}\n".format(name_to_keep,  
                                                                       stable_id_1, seq_name_1, start_1, end_1, 
                                                                       stable_id_2, seq_name_2, start_2, end_2)
    if notes:
        header  = "% The following pieces of sequence were found split across different scaffolds/contigs\n" 
        header += "% and assumed to actually belong to the same gene:\n" 
        notes   = "\n"+header+notes

    return notes

#########################################
def remove_dubious_paralogues (cursor, ensembl_db_name, output_pep, sequence_name_to_exon_names, human_exon_to_ortho_exon):

    notes = ""
    mulitple_orthos = []

    if ( isinstance(output_pep,str)) : return notes

    for seq_name in output_pep.keys():
        name_pieces = seq_name.split("_")
        if not isinteger(name_pieces[-1]): continue
        trivial_name = "_".join(name_pieces[:-1])
        if trivial_name not in mulitple_orthos:  mulitple_orthos.append(trivial_name)
    
    for trivial_name in mulitple_orthos:
        # find_mammals can handle a list of names, but we have only one here
        if not find_mammals(cursor, [trivial_name]): continue
        paralogues = filter (lambda seq_name: trivial_name in seq_name,  output_pep.keys())
        tanimoto   = {}
        max_tanimoto = -1
        for para in  paralogues:
            # what is the is the  highest tanimoto with human?
            tanimoto[para] = pairwise_tanimoto (output_pep['human'],
                                                output_pep[para], use_heuristics=False)
            # keep everything within say 90% of the max
            if ( max_tanimoto < tanimoto[para]): 
                max_tanimoto = tanimoto[para]
        # sort paralogues by their tanimoto score
        sorted_para = sorted(paralogues, key=lambda para: -tanimoto[para])
        ct = 0
        tmp_names     = []
        dropped_paras = []
        for para in sorted_para:

            if tanimoto[para] < 0.9*max_tanimoto:
                # drop
                [exon_id, exon_known] = sequence_name_to_exon_names[para][0].split ("_")[-3:-1] # the last number is the start in the gene
                species   = "_".join (sequence_name_to_exon_names[para][0].split ("_")[:-3])   
                gene_id   = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
                stable_id = gene2stable(cursor, gene_id, ensembl_db_name[species])
                dropped_paras.append(stable_id)
                # 
                del output_pep[para]
                del sequence_name_to_exon_names[para]
                del human_exon_to_ortho_exon[para]

            else:
                ct += 1
                tmp_name = "tmp"
                if ct > 1: tmp_name += "_"+str(ct)
                output_pep[tmp_name] = output_pep.pop(para)
                sequence_name_to_exon_names[tmp_name] = sequence_name_to_exon_names.pop(para)
                human_exon_to_ortho_exon[tmp_name] = human_exon_to_ortho_exon.pop(para)
                tmp_names.append(tmp_name)

        ct = 0
        for tmp_name in tmp_names:
            ct += 1
            new_name = trivial_name
            if ct > 1: new_name += "_"+str(ct)
            output_pep    [new_name] = output_pep.pop(tmp_name)
            sequence_name_to_exon_names[new_name] = sequence_name_to_exon_names.pop(tmp_name)
            human_exon_to_ortho_exon[new_name] = human_exon_to_ortho_exon.pop(tmp_name)
            
        if dropped_paras:
            notes += "for {0}: ".format(trivial_name)
            first = True
            for stable_id in dropped_paras:
                # find stable id
                if not first:
                    notes += ";"
                notes += stable_id
                first = False
            notes += "\n"
    if notes:
        header  = "% The following sequences, labeled in Ensembl as one2many orthologues,  were dropped\n"
        header += "% because they were deemed problematic: too short, or too different compared to human sequence\n"
        notes  = "\n" + header + notes
    return notes

#########################################
def sort_trivial_names (cursor, all_species):
    # in the afa headers use 'trivial' names for the species: cow, dog, pig, ...
    trivial_name   = translate_to_trivial(cursor, all_species)

    # walk the taxonomical tree, and sort the species according to
    # the (distance of) the last common ancestor
    sorted_species       = {}
    sorted_trivial_names = {}
    for qry_species in ['homo_sapiens']:
        sorted_species[qry_species] = species_sort(cursor, all_species, qry_species)
        trivial = trivial_name[qry_species]
        sorted_trivial_names[trivial] = map(lambda species: trivial_name[species], sorted_species[qry_species])

    return [trivial_name, sorted_trivial_names]

#########################################
def check_afa_age (cfg, stable_id):

    max_days = 30

    afa_age = "old"
    afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
    if (os.path.exists(afa_fnm) and os.path.getsize(afa_fnm) > 0 ):
        time_modified = os.path.getmtime(afa_fnm)
        number_of_days_since_modified = (time.time() - time_modified)/(60*60*24)
        if number_of_days_since_modified < max_days:
            print "\t %s last modified %s. Moving on." % (stable_id, time.ctime(os.path.getmtime(afa_fnm) ))
            afa_age  = "new"
    return afa_age

#########################################
def make_atlas(cursor, ensembl_db_name, canonical_human_exons, alnmt_pep, trivial_name):
    # >>>>>>>>>>>>>>>>>>
    # find which species we have, and for how many exons
    # we may have two orthologues from the same species
    seq_name = {}
    parent_seq_name  = {}
    for human_exon in canonical_human_exons:
        if not alnmt_pep[human_exon]: continue # -- this should not happen; we should have at least human exons
        for exon_seq_name in alnmt_pep[human_exon].keys():
            [species, exon_id, exon_known, exon_start] = parse_aln_name(exon_seq_name)
            ortho_gene_id                  = exon_id2gene_id(cursor, ensembl_db_name[species], exon_id, exon_known)
            # gene --> name -- retrieve old name, or construct a new one
            parent_seq_name[exon_seq_name] = get_name (seq_name, trivial_name[species], ortho_gene_id) 
                
    # >>>>>>>>>>>>>>>>>>
    sequence_name_to_exon_names = {}
    human_exon_to_ortho_exon = {}
    for human_exon in canonical_human_exons:
        if not alnmt_pep[human_exon]: continue
        for exon_seq_name in alnmt_pep[human_exon].keys():
            concat_seq_name = parent_seq_name[exon_seq_name]

            if not sequence_name_to_exon_names.has_key(concat_seq_name): sequence_name_to_exon_names[concat_seq_name] = []
            if not exon_seq_name in  sequence_name_to_exon_names[concat_seq_name]:
                sequence_name_to_exon_names[concat_seq_name].append(exon_seq_name)
                    
            if not human_exon_to_ortho_exon.has_key(concat_seq_name): human_exon_to_ortho_exon[concat_seq_name] = {}
            if not human_exon_to_ortho_exon[concat_seq_name].has_key(human_exon): 
                human_exon_to_ortho_exon[concat_seq_name][human_exon] = []

            human_exon_to_ortho_exon[concat_seq_name][human_exon].append(exon_seq_name)

    # make sure we are sorted, will come handy down the line
    for seq_name in sequence_name_to_exon_names.keys():
        sort_exon_names (sequence_name_to_exon_names[seq_name]  )

    # >>>>>>>>>>>>>>>>>>
    # flag the cases when one orthologue exon maps to many human (and vice versa) for later
    # do we have a sequence mapping to multiple human exons?
    ortho_exon_to_human_exon = {}
    for human_exon in canonical_human_exons:
        if not alnmt_pep[human_exon]: continue
        for exon_seq_name in alnmt_pep[human_exon].keys():
            if not ortho_exon_to_human_exon.has_key(exon_seq_name):
                ortho_exon_to_human_exon[exon_seq_name] = [human_exon]
            else:
                ortho_exon_to_human_exon[exon_seq_name].append(human_exon)
    # 
    overlapping_maps = {}
    for concat_seq_name, concat_exons in sequence_name_to_exon_names.iteritems():
        overlapping_maps[concat_seq_name] = find_overlapping_maps (ortho_exon_to_human_exon, concat_exons, alnmt_pep)
        # could the overlapping maps end up not being sorted?
        for [human_exons, ortho_exons] in overlapping_maps[concat_seq_name]:
            sort_exon_names (ortho_exons)
            human_exons.sort(key=lambda exon: exon.start_in_gene)

    return [human_exon_to_ortho_exon, sequence_name_to_exon_names, ortho_exon_to_human_exon, overlapping_maps]





#########################################
def make_exon_alignments(cursor, ensembl_db_name, canonical_human_exons,
                         mitochondrial, min_similarity, flank_length):
    alnmt_pep = {}
    alnmt_dna = {}
    first_human_exon = True
    for human_exon in canonical_human_exons:
        # make_exon_alignment defined in el_utils/el_specific.py
        # we need the info about the first human exon
        # to get a bit more lenient whent the first exon consists of M only
        # (I'll have to take a look at one point what's up with that)
        [alnmt_pep[human_exon], alnmt_dna[human_exon]]  =   make_exon_alignment(cursor, ensembl_db_name,  human_exon.exon_id, 
                                                                                human_exon.is_known,  mitochondrial, min_similarity, 
                                                                                flank_length, first_human_exon)   
        first_human_exon = False

    return [alnmt_pep, alnmt_dna] 


#########################################
#########################################
#########################################
def make_alignments ( gene_list, db_info):

    [local_db, ensembl_db_name] = db_info

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
    # find triviaal names (cow, pig etc) and sort them according to tax distance from human
    [trivial_name, sorted_trivial_names] = sort_trivial_names(cursor, all_species)
    # find minimum acceptable similarity between mapped exons (used in this pipeline)
    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    min_similarity = cfg.get_value('min_accptbl_exon_sim') 

    ##########################################################################
    # for each  gene in the provided list
    for gene_id in gene_list:

        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)
        # if we are running this pipe repetedly we want to skip if
        # the last time we worked on this gene was recently enough
        #if  check_afa_age (cfg, stable_id) == "new": continue                               
 
        if verbose: 
            print gene_id, stable_id, get_description (cursor, gene_id)
        elif (not gene_list.index(gene_id)%100): 
            print gene_list.index(gene_id), "out of ", len(gene_list)

        # mitochondrial?
        mitochondrial = is_mitochondrial(cursor, gene_id)
        # find all canonical coding  human exons 
        canonical_human_exons = filter (lambda x:  x.is_canonical and x.is_coding, gene2exon_list(cursor, gene_id))
        # bail out if there is a problem
        if not canonical_human_exons: continue
        # the exons are not guaranteed to be in order
        canonical_human_exons.sort(key=lambda exon: exon.start_in_gene)
        # reconstruct  per-exon alignments with orthologues
        [alnmt_pep, alnmt_dna] = make_exon_alignments(cursor, ensembl_db_name, canonical_human_exons,
                                                      mitochondrial, min_similarity, flank_length)

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # we want to be able to retrieve the info starting from whichever end, so we construct the following maps:
        # to find all exons from an ortohologue, that map to a given human exon:
        #      human_exon_to_ortho_exon:  human_exon_to_ortho_exon[concat_seq_name][human_exon].append(exon_seq_name)
        # given a sequence name, retrieve all exons that belong to it
        #     sequence_to_exon: sequence_name_to_exon_names[concat_seq_name].append(exon_seq_name)
        # in the other direction - to find all human exons that a given exon  from orthologous sequence maps to
        #      ortho_exon_to_human_exon:  ortho_exon_to_human_exon[exon_seq_name].append(human_exon)
        # finally, collect in one place info about maps between human and orthologue that are not one-to-one
        #      overlapping_maps: overlapping_maps[concat_seq_name].append([human_exons, ortho_exons])
        [human_exon_to_ortho_exon, sequence_name_to_exon_names, 
         ortho_exon_to_human_exon, overlapping_maps] = make_atlas(cursor, ensembl_db_name, canonical_human_exons, 
                                                                  alnmt_pep, trivial_name)
        # the alignment always has human sequence, but if it is the only one
        # (see for example RPL41, ENSG00000229117, a 25 residue peptide,  for which NCBI REfseq
        # reports a single  confirmed  homologue in mouse, but Ensembl reports no orthologues at all)
        if ( len(sequence_name_to_exon_names) <= 1):
            print "\t no orthologues for",  gene_id, stable_id, " (?)"
            continue
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # concatenate the aligned exons for each species, taking into account that the mapping
        # doesn't have to be one to one
        headers     = []
        output_pep  = {}
        output_dna  = {}
        output_pep_ok = True
        for concat_seq_name in sequence_name_to_exon_names.keys():
            
            if not human_exon_to_ortho_exon.has_key(concat_seq_name): continue # this shouldn't happen but oh well

            output_pep[concat_seq_name] = ""
            output_pep_ok = True

            # compile human exons which are in one-to-many relationship in this species (this sequence, actually)
            flagged_human_exons = set([])
            if overlapping_maps.has_key(concat_seq_name):
                for  [human_exons, ortho_exons] in overlapping_maps[concat_seq_name]:
                    if len(human_exons) == 1 and len(ortho_exons) == 1: continue 
                    flagged_human_exons |= set(human_exons)

            for human_exon in canonical_human_exons:

                if ( not alnmt_pep.has_key(human_exon)): 
                    continue
                if ( isinstance(alnmt_pep[human_exon], str)): continue
                
                aln_length = len(alnmt_pep[human_exon].itervalues().next())
                pep = '-'*aln_length
                if human_exon in flagged_human_exons:
                    # one ortho seq maps to multiple human exons, or vice versa
                    pep = '-'*aln_length

                elif human_exon_to_ortho_exon[concat_seq_name].has_key(human_exon):
                    # we have a neat one-to-one mapping
                    exon_seq_name = human_exon_to_ortho_exon[concat_seq_name][human_exon][0]
                    pep = alnmt_pep[human_exon][exon_seq_name]
                    if not pep:  # what's this?
                        pep = '-'*aln_length
                else: 
                    # no exon in this species
                    pep =  '-'*aln_length
                   
                if output_pep[concat_seq_name]: output_pep[concat_seq_name] += '-Z-'
                output_pep[concat_seq_name] += pep
                if not output_pep:  
                    c=inspect.currentframe()
                    print " in %s:%d output peptide empty" % (c. f_code.co_filename, c.f_lineno)
                    output_pep_ok = False
                    # Oct 13: I am not sure of the full implication of this, so I'll just abort

            if output_pep_ok:  headers.append(concat_seq_name)


        #########################################################
        if not output_pep_ok: continue

        #########################################################
        # >>>>>>>>>>>>>>>>>>
        sorted_seq_names = sort_names (sorted_trivial_names['human'], output_pep)
        boundary_cleanup(output_pep, sorted_seq_names)
        output_pep = strip_gaps(output_pep)

        assorted_notes = ""
        for seq_to_fix in overlapping_maps.keys():
            #if not seq_to_fix=='chimpanzee': continue
            # fix_one2many changes both output_pep and sequence_name_to_exon_names
            if not overlapping_maps[seq_to_fix]: continue
            [output_pep, sequence_name_to_exon_names] = fix_one2many (cursor, ensembl_db_name, cfg, acg, sorted_seq_names, 
                                                         canonical_human_exons, human_exon_to_ortho_exon, 
                                                         sequence_name_to_exon_names, seq_to_fix, 
                                                         overlapping_maps[seq_to_fix], alnmt_pep, output_pep)

        # check if any two pieces of seqeunce ended up on different scaffolds/contigs
        fusion_notes = fuse_seqs_split_on_scaffolds(cursor, acg, ensembl_db_name,  output_pep, sequence_name_to_exon_names, 
                                     ortho_exon_to_human_exon, canonical_human_exons, human_exon_to_ortho_exon)

        assorted_notes += fusion_notes + "\n"
        # get rid of dubious paralogues (multiple seqs from the same species)
        para_notes = remove_dubious_paralogues (cursor, ensembl_db_name, output_pep, sequence_name_to_exon_names, human_exon_to_ortho_exon)
        assorted_notes += para_notes + "\n"
        # we may have chosen to delete some sequences
        sorted_seq_names = sort_names (sorted_trivial_names['human'], output_pep)

        if not check_seq_length (output_pep, "ouput_pep"): 
            print "length check failure"
            #continue

        if (1):
            # >>>>>>>>>>>>>>>>>>
            boundary_cleanup(output_pep, sorted_seq_names)
            output_pep = strip_gaps(output_pep)
            if not output_pep:
                print "no output pep 3: moving on"
                continue

        if (1):
            # >>>>>>>>>>>>>>>>>>
            output_dna = expand_protein_to_dna_alnmt (cursor, ensembl_db_name, cfg, acg, 
                                                      sorted_trivial_names, sequence_name_to_exon_names,  
                                                      alnmt_pep, output_pep, flank_length)
            if not output_dna:
                print "no output dna: moving on"
                continue

            output_dna = strip_gaps(output_dna)
            afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
            ret = output_fasta (afa_fnm, sorted_seq_names, output_dna)
            #if not ret: continue
            if verbose: print afa_fnm

        # >>>>>>>>>>>>>>>>>>
        output_pep = fix_split_codons (cursor, ensembl_db_name, cfg, acg, 
                                           sorted_trivial_names, mitochondrial, sequence_name_to_exon_names,  
                                           alnmt_pep, output_pep, flank_length)

        if not output_pep:
            print "no output pep 4: moving on"
            continue

        afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        ret = output_fasta (afa_fnm, sorted_seq_names, output_pep)
 
        if verbose: print afa_fnm
            
        # notes to accompany the alignment:
        print_notes (cursor, cfg,  ensembl_db_name, output_pep, sequence_name_to_exon_names,  
                     sorted_seq_names, stable_id, human_exon_to_ortho_exon, assorted_notes)
       
    return 


#########################################
def main():
    
    no_threads = 1
    special    = 'test'

    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <set name> <number of threads> <method>" % sys.argv[0]
        exit(1) # exit after the usage statement
    elif len(sys.argv)==3:

        special = sys.argv[1]
        special = special.lower()
        if special == 'none': special = None

        no_threads = int(sys.argv[2])

    local_db   = False
    
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    print "running ", sys.argv[0]
    if special:

        print "using", special, "set"
        if special == 'complement':
            gene_list = get_complement_ids(cursor, ensembl_db_name, cfg)
        else:
            gene_list = get_theme_ids (cursor, ensembl_db_name, cfg, special )
    else:
        print "using all protein coding genes"
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        
    cursor.close()
    db.close()

    parallelize (no_threads, make_alignments, gene_list, [local_db, ensembl_db_name])
    
    return True


#########################################
if __name__ == '__main__':
    main()


