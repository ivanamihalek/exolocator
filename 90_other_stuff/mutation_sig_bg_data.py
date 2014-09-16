#!/usr/bin/python -u
# extract all info one needs to estimate the background mutation rate when analyzing tumors

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
def align_nucseq_by_pepseq(aligned_pepseq, nucseq):
    if (not len(aligned_pepseq.replace('-',''))*3 == len(nucseq)):
        print "length mismatch:", len(aligned_pepseq.replace('-',''))*3,  len(nucseq)
        return " -xxxx- "
    codons = iter(map(''.join, zip(*[iter(nucseq)]*3)))
    aligned_nucseq = ''.join(('---' if c=='-' else next(codons) for c in aligned_pepseq))
    return aligned_nucseq

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

    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1, ref_only=True)

    # for each human gene
    gene_ct = 0
    for gene_id in gene_ids[:3]:
       
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)

        gene_ct += 1
        if (not gene_ct%100): print gene_ct, "out of ", len(gene_ids)
        if verbose: print gene_id, stable_id, get_description (cursor, gene_id)

        # find all canonical coding  human exons 
        # get_canonical_coding_exons also sorts exons by the start in the gene
        canonical_human_exons = get_canonical_coding_exons (cursor, gene_id, ensembl_db_name['homo_sapiens'])
        # bail out if there is a problem
        if not canonical_human_exons: continue

        # reconstruct the alignment with orthologues
        sequence  = {}
        seq_name  = {}
        has_a_map = False
        for human_exon in canonical_human_exons:

            maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
            if not maps:
                continue
            else:
                has_a_map = True

            for map in maps:
                species = map.species_2
                # get the raw (unaligned) sequence for the exon that maps onto human
                [exon_seq_id, unaligned_sequence, pepseq_transl_start, pepseq_transl_end, left_flank, right_flank, nucseq] = \
                    get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
                # inflate the compressed sequence
                if map.bitmap and unaligned_sequence:
                    bs = Bits(bytes=map.bitmap)
                    # check bitmap has correct number of 1s
                    if ( not bs.count(1) == len(unaligned_sequence)):
                        print "bitmap check fails (?)"
                        continue
                    # rebuild aligned sequence
                    usi = iter(unaligned_sequence)
                    reconstructed_sequence = "".join(('-' if c=='0' else next(usi) for c in bs.bin))
                    # rebuild aligned dna sequence, while we are at that
                    reconstructed_nucseq   = align_nucseq_by_pepseq(reconstructed_sequence, nucseq)

                    print unaligned_sequence
                    print nucseq
                    print reconstructed_nucseq
                    print
                else:
                    print fail
                    continue
 

#########################################
if __name__ == '__main__':
    main()

