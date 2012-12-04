#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re

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
def merged_sequence (template_seq, sequence_pieces, nucseq_pieces, flank_length):
    
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
            print "length mismatch: ", len(aligned_pepseq.replace('-',''))*3, len(nucseq)
            exit (1)
	codons = iter(map(''.join, zip(*[iter(nucseq)]*3)))
	aligned_nucseq = ''.join(('---' if c=='-' else next(codons) for c in aligned_pepseq))
	return aligned_nucseq

#########################################
def expand_aligned_pepseq(cursor, aligned_pepseq, exon_seqs,  exon_id, exon_is_known, flank_length):

    # All sequences refer to the same exon:
    # aligned_pepseq is known, and so is unaligned_pepseq and the original dna it came from
    # from this the positions of the gaps in the dna version should be reconstructed.
    # WHen stripped of the gaps, the aligned_pepseq  should be the same as the unaligned_pepseq.
    # There are two problems, however  with the correspondence between the unaligned_pepseq and dnda:
    #   (i) the translation is sometimes not possible in the region advertised by Ensembl
    #  (ii) the nucleotides needed to complete the codon have been added from the previous exon.
    # To be on the safe side, in this step at least, we reconstruct the translation step/

    [unaligned_pepseq, left_flank, right_flank, dna_seq] = exon_seqs
    empty = ["", "", ""]


    if len(dna_seq) < 10:
        return empty

    aligned_dna_seq = '-'*(10+3*len(aligned_pepseq)+10)

    orig_unaligned_pepseq = unaligned_pepseq

    exon                 = get_exon           (cursor, exon_id, exon_is_known)
    mitochondrial        = is_mitochondrial   (cursor, exon.gene_id)
    [seq_start, seq_end] = translation_bounds (cursor, exon_id)
    dna_cropped          = crop_dna           (seq_start, seq_end, dna_seq)
    if ( exon.phase > 0) :
        offset = phase2offset (exon.phase)
        unaligned_pepseq = unaligned_pepseq[1:]
        aligned_pepseq   = re.sub ("[A-Z]", '-', aligned_pepseq, 1)
    else:
        offset = 0

    if seq_start is None: 
        seq_start = 0
    elif seq_start > 0:
        seq_start -= 1
    if seq_end   is None: 
        seq_end   = 0
    elif seq_end > 0:
        seq_end -= 1

    # this is kinda crass, but the most straightforward
    # starting from the suggested phase, translate the
    # the dna sequence in all three frames, until the unaligned_pepseq is found
    done = False
    ct   = 0


    match_found = 0
    while not done:

        offset = offset%3
        dnaseq = Seq (dna_cropped[offset:], generic_dna)
        if mitochondrial:
            pepseq = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
        else:
            pepseq = dnaseq.translate().tostring()
            
        if ( pepseq[-1] == '*'):
            pepseq = pepseq[:-1]

        if ( 'U' in unaligned_pepseq):
            pepseq  = pepseq.replace('*', 'U')
 
        if ( not '*' in pepseq):

            if (unaligned_pepseq in pepseq):
                # where exactly does the match start?
                pattern = re.compile(unaligned_pepseq)
                for match in pattern.finditer(pepseq):
                    if ( match.end()-match.start()== len(unaligned_pepseq)):
                        match_found += 1
                        match_start = seq_start+offset+3*match.start()
                        #print "<<<", match.start(), match.end(), match_start, offset
                        #print "ul ", unaligned_pepseq
                        #print "pep", pepseq
                        #print "----------"

            elif pepseq in unaligned_pepseq:
                pattern = re.compile(pepseq)
                for match in pattern.finditer(unaligned_pepseq):
                    if ( match.end()-match.start()== len(unaligned_pepseq)):
                        match_found += 1
                        match_start = seq_start+offset+3*match.start()
                        #print ">>>"
                    
        else:
            nonstop =  re.compile("[A-Z]+")
            for subseq in re.findall (nonstop, pepseq):
                if len(subseq) < 4: continue
                pattern = re.compile(subseq)

                for match in pattern.finditer(unaligned_pepseq):
                    if (match.end()== len(unaligned_pepseq)):
                        match_found += 1
                        match_start = seq_start+offset+3*match.start()+3
                        #print "ooo"
        ct     += 1
        offset += 1


        done = (match_found or ct==3)

    if (not match_found == 1):
        print "------------------------------"
        print unaligned_pepseq
        print "number of  matches ", match_found
        print
        return empty
    else:
        # some more sanity checking
        match_end  = match_start + 3*len(unaligned_pepseq)
        dnaseq     = Seq (dna_seq[match_start:match_end], generic_dna)
        if mitochondrial:
            pepseq = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
        else:
            pepseq = dnaseq.translate().tostring()

        if ( 'U' in unaligned_pepseq):
            pepseq  = pepseq.replace('*', 'U')

        if (not pepseq == unaligned_pepseq):
            print exon.phase, offset
            print orig_unaligned_pepseq
            print unaligned_pepseq
            print pepseq
            exit(1)
        # expand
        aligned_dna_seq = align_nucseq_by_pepseq(aligned_pepseq, dnaseq.tostring())

        # add the left and the right flank

        if (match_start > flank_length):
            effective_left_flank = dna_seq[match_start-flank_length:match_start]
        elif (match_start == 0):
            effective_left_flank = left_flank[-flank_length:].lower()
        else:
            effective_left_flank = left_flank[-(flank_length-match_start):].lower()+\
                dna_seq[:match_start]
            
        delta = len(dna_seq)-match_end 
        if ( delta  > flank_length):
            effective_right_flank = dna_seq[match_end:match_end+flank_length]
        elif (match_end == len(dna_seq)):
            effective_right_flank = right_flank[:flank_length].lower()
        else:
            effective_right_flank = dna_seq[match_end:]+right_flank[:flank_length-delta].lower()
            
        # pad the flanking seqs to the needed length
        effective_left_flank  =  effective_left_flank.rjust(flank_length, '-')
        effective_right_flank = effective_right_flank.ljust(flank_length, '-')


    # whatever we return must be 3*the input peptide long + the flanks
    return [effective_left_flank, aligned_dna_seq, effective_right_flank]


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
    for gene_id in gene_list:
       
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)

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
        sequence     = {}
        dna_sequence = {}
        seq_name     = {}
        has_a_map    = False
        for human_exon in canonical_exons:

            # find all entries in the exon_map with this exon id
            # (= find all orthologous exons that the human exon maps to)
            maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
            if not maps:
                continue
            else:
                has_a_map = True

            for map in maps:
                if not map.bitmap: continue
                species = map.species_2
                #if not species=='petromyzon_marinus': continue
                switch_to_db (cursor,  ensembl_db_name[species])

                # get the raw (unaligned) sequence for the exon that maps onto human
                exon_seqs = get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
                exon_seqs = exon_seqs[1:] # the first entry is database id
                [pep_seq, left_flank, right_flank, dna_seq] = exon_seqs

                # inflate the compressed sequence
                unaligned_sequence = pep_seq
                if not unaligned_sequence:continue

                bs = Bits(bytes=map.bitmap)
                # check bitmap has correct number of 1s
                if ( not bs.count(1) == len(unaligned_sequence)):
                    print "bitmap check fails (?)"
                    continue
                # rebuild aligned sequence
                usi = iter(unaligned_sequence)
                aligned_pep_seq = "".join(('-' if c=='0' else next(usi) for c in bs.bin))
                # rebuild aligned dna sequence, while we are at that
                [left_flank, aligned_dna_seq, right_flank] =  expand_aligned_pepseq (cursor, aligned_pep_seq,  
                                                                                     exon_seqs, map.exon_id_2, 
                                                                                     map.exon_known_2, flank_length)
                #
                reconstructed_sequence = aligned_pep_seq
                reconstructed_nucseq   = left_flank +aligned_dna_seq+right_flank

                # come up with a unique name for this sequence
                species = map.species_2                
                # exon --> gene 
                ortho_gene_id = exon_id2gene_id(cursor, ensembl_db_name[species], 
                                               map.exon_id_2, map.exon_known_2)
                # gene --> name
                if ( not seq_name.has_key(species) ): 
                    sequence_name     = species
                    seq_name[species] = {ortho_gene_id:sequence_name}
                else:
                    if seq_name[species].has_key(ortho_gene_id):
                        sequence_name = seq_name[species][ortho_gene_id]
                    else:
                        # if sequence for the species exists, 
                        # but not for this gene, mangle the name
                        number_of_genes = len(seq_name[species].keys()) + 1
                        sequence_name   = species + "_" + str(number_of_genes)
                        seq_name[species][ortho_gene_id] = sequence_name

                # use human exon as a label, no matter which species the sequence is
                # actually coming from 
                if not sequence.has_key(sequence_name):
                    sequence[sequence_name]     = {}
                    dna_sequence[sequence_name] = {}
 
                if not sequence[sequence_name].has_key(human_exon):
                    sequence[sequence_name][human_exon]     = reconstructed_sequence
                    dna_sequence[sequence_name][human_exon] = reconstructed_nucseq
                else:
                    # there is a number of sequences that map onto this human sequence
                    if type(sequence[sequence_name][human_exon]) is list:
                        sequence[sequence_name][human_exon].append(reconstructed_sequence)
                        dna_sequence[sequence_name][human_exon].append(reconstructed_nucseq)
                    else: # we need to start storing as a list
                        tmp = sequence[sequence_name][human_exon]
                        sequence[sequence_name][human_exon]     = [tmp, reconstructed_sequence]
                        
                        tmp = dna_sequence[sequence_name][human_exon]
                        dna_sequence[sequence_name][human_exon] = [tmp, reconstructed_nucseq]

        # >>>>>>>>>>>>>>>>>>
        if not has_a_map: continue

        # how the hell this happened:
        if not sequence.has_key('homo_sapiens'):
            print "no human seq found in the alignment!!!!!!!"
            continue

        to_remove = []
        for  i in range(len(canonical_exons)): 
            human_exon = canonical_exons[i]
            # not sure how could this happen, but ...
            if (not sequence['homo_sapiens'].has_key(human_exon) ):
                to_remove.append(i)

        for i in range( len(to_remove)-1, -1, -1):
            del canonical_exons[to_remove[i]]

        ###############################################################
        # stitch  the exons together
        headers        = []
        output_pep     = {}
        output_dna     = {}
        ortholog_count = {}
        for species in sorted_species:
            if not seq_name.has_key(species):
                continue
            trivial         = trivial_name[species]
            # how many genes homologous to the humna query do we have in this particular species?
            number_of_genes = len(seq_name[species].keys())
            
            for ct in range (1,number_of_genes+1):

                # construct the name for this gene (to be used in afa file), and put it in the header
                sequence_name = species
                if (ct > 1):
                    sequence_name = species + "_" + str(ct)
                new_name      = trivial
                if (number_of_genes > 1):
                    new_name  = trivial + "_" + str(ct)                
                headers.append(new_name)

                # reconstruct the full length alignment from exon alignments
                output_pep[new_name] = ""
                output_dna[new_name] = ""
                first_exon = True
                for  human_exon in canonical_exons:
                    if not first_exon:
                        output_pep[new_name] += 'Z'
                        output_dna[new_name] += 'Z'
                    # this species has a counterpart human exon
                    if sequence[sequence_name].has_key(human_exon):

                        # multiple candidates for the exon that the human exon maps to
                        if ( type(sequence[sequence_name][human_exon]) is list):
                            template = sequence['homo_sapiens'][human_exon]
                            # merged will delete the pieces that are covered and less similar to the template
                            # it will also return the dna version of the merged sequence
                            [pep, dna] = merged_sequence (template, sequence[sequence_name][human_exon], 
                                                          dna_sequence[sequence_name][human_exon], flank_length)
                            output_pep[new_name] += pep
                            output_dna[new_name] += dna
                        # a single homologue
                        else:
                            output_pep[new_name] +=     sequence[sequence_name][human_exon]
                            output_dna[new_name] += dna_sequence[sequence_name][human_exon]

                    # the exon in this species is an insert wrt human sequence
                    else:
                        if (not sequence['homo_sapiens'].has_key(human_exon) ):
                            print ct, human_exon
                            exit(1)
                        output_pep[new_name] += '-'*len(    sequence['homo_sapiens'][human_exon])
                        output_dna[new_name] += '-'*len(dna_sequence['homo_sapiens'][human_exon])
                    first_exon = False
                    

        # find place for best_afa -- for the moment can put it to the scratch space:
        afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        output_fasta (afa_fnm, headers, output_pep)
        print afa_fnm
        afa_fnm  = "{0}/dna/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        output_fasta (afa_fnm, headers, output_dna)
        print afa_fnm

        continue

        # notes to accompany the alignment:
        notes_fnm  = "{0}/notes/{1}.txt".format(cfg.dir_path['afs_dumps'], stable_id)
        print notes_fnm
        print_notes (notes_fnm, orthologues, exons, sorted_species, specid2name, human_stable_id, source)
        
        

#########################################
def main():
    
    no_threads = 10

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
