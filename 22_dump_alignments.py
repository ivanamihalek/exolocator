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

from bitstring import Bits


#########################################
def translate_to_trivial(cursor, all_species):

    trivial_name = {}
    for species in all_species:
        taxid                 = species2taxid (cursor, species)
        trivial_name[species] = taxid2trivial(cursor, taxid)

    return trivial_name

#########################################
def  fract_identity (seq1, seq2):
    
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
            print i, deletable
            del sequence_pieces[deletable]

    # if not, go ahead and merge
    merged = ""
    for pos in range(template_length):
        new_char = '-'
        for piece in sequence_pieces:
            if piece[pos] == '-': 
                continue
            else:
                new_char =  piece[pos]
                break
        merged += new_char
    
    return merged   


#########################################
def main():

    verbose  = False
    local_db = False

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
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    # for each human gene
    for gene_id in gene_ids:
    #for gene_id in [412667]: #  wls
    #for gene_id in [378768]: #  p53
    #for gene_id in [378766]: #  dynein
       
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)
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

        # one2one   orthologues
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        known_orthologues      = get_orthos (cursor, gene_id, 'orthologue')
        # not-clear orthologues
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        unresolved_orthologues = get_orthos (cursor, gene_id, 'unresolved_ortho')

        # reconstruct the alignment with orthologues
        sequence  = {}
        seq_name  = {}
        has_a_map = False
        for human_exon in canonical_exons:

            # find all entries in the exon_map with this exon id
            # (= find all orthologous exons that the human exon maps to)
            maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
            if not maps:
                continue
            else:
                has_a_map = True
            for map in maps:
                species = map.species_2
                # get the unaligned sequence
                unaligned_sequence = get_exon_pepseq(cursor, map.exon_id_2, map.exon_known_2, 
                                                     ensembl_db_name[map.species_2])
                #print "###################################"
                #print map
                #print get_exon_pepseq(cursor, map.exon_id_1, map.exon_known_1, ensembl_db_name['homo_sapiens'])
                #print get_exon_pepseq(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[species])
                #print 

                if map.bitmap and unaligned_sequence:
                    bs = Bits(bytes=map.bitmap)
                    # check bitmap has correct number of 1s
                    if ( not bs.count(1) == len(unaligned_sequence)):
                        print "bitmap check fails (?)"
                        continue
                    # rebuild aligned sequence
                    usi = iter(unaligned_sequence)
                    reconstructed_sequence = "".join(('-' if c=='0' else next(usi) for c in bs.bin))
                else:
                    #print " No bitmap! Error in the database! Skip past this entry"
                    #print map.species_2
                    #print map.exon_id_2
                    continue

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
                        sequence_name   = seq_name[species][ortho_gene_id]
                    else:
                        # if sequence for the species exists, 
                        # but not for this gene, mangle the name
                        number_of_genes = len(seq_name[species].keys()) + 1
                        sequence_name   = species + "_" + str(number_of_genes)
                        seq_name[species][ortho_gene_id] = sequence_name

                # use human exon as a label, no matter which species the sequence is
                # actually coming from 
                if not sequence.has_key(sequence_name):
                    sequence[sequence_name] = {}
 
                if not sequence[sequence_name].has_key(human_exon):
                    sequence[sequence_name][human_exon] = reconstructed_sequence
                    
                else:
                    # there is a number of sequences that map onto this human sequence
                    if type(sequence[sequence_name][human_exon]) is list:
                        sequence[sequence_name][human_exon].append(reconstructed_sequence)
                    else: # we need to start storing as a list
                        tmp = sequence[sequence_name][human_exon]
                        sequence[sequence_name][human_exon] = [tmp, reconstructed_sequence]

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

          
        # stitch together different exons
        headers        = []
        output_seq     = {}
        ortholog_count = {}
        for species in sorted_species:
            if not seq_name.has_key(species):
                continue
            trivial         = trivial_name[species]
            number_of_genes = len(seq_name[species].keys())
            
            for ct in range (1,number_of_genes+1):

                sequence_name = species
                if ( ct > 1):
                    sequence_name = species + "_" + str(ct)
                new_name      = trivial
                if ( number_of_genes > 1):
                    new_name = trivial + "_" + str(ct)
                
                headers.append(new_name)
                output_seq[new_name] = ""
                first_exon = True
                for  human_exon in canonical_exons:
                    if not first_exon:
                        output_seq[new_name] += 'Z'
                    if sequence[sequence_name].has_key(human_exon):
                        if (  type(sequence[sequence_name][human_exon]) is list):
                            template = sequence['homo_sapiens'][human_exon]
                            # merged will delete the pieces that are covered and less similar to the template
                            output_seq[new_name] +=  merged_sequence (template, sequence[sequence_name][human_exon])

                        else:
                            output_seq[new_name] += sequence[sequence_name][human_exon]
                    else:
                        if (not sequence['homo_sapiens'].has_key(human_exon) ):
                            print ct, human_exon
                            exit
                        output_seq[new_name] += '-'*len(sequence['homo_sapiens'][human_exon])
                    first_exon = False
                    

        # find place for best_afa -- for the moment can put it to the scratch space:
        afa_fnm  = "{0}/pep/{1}.afa".format(cfg.dir_path['afs_dumps'], stable_id)
        output_fasta (afa_fnm, headers, output_seq)
        print afa_fnm

        # reconstruct the dna alignment - output (? how are we going to do that?)
        # there is exon_seq table in each genome database, containing both the dna and the peptide sequence ...
        # as well as the left and the right flanks of the exon


#########################################
if __name__ == '__main__':
    main()
