#!/usr/bin/python -u
# extract all info one needs to estimate the background mutation rate when analyzing tumors

import MySQLdb, commands, re

from   el_utils.mysql       import  *
from   el_utils.ensembl     import  *
from   el_utils.el_specific import  *
from   el_utils.utils   import  erropen
from   el_utils.map     import  Map, get_maps
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader

from Bio.Seq      import Seq
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

    verbose  = True
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
    
    logf = erropen("error.log", "w") 
    if not logf: exit(1)

    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    #gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1, ref_only=True)

    # for each human gene
    gene_ids = [10093176]
    gene_ct = 0
    for gene_id in gene_ids[:100]:
       
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)

        gene_ct += 1
        if (not gene_ct%100): print gene_ct, "out of ", len(gene_ids)

        # find all canonical coding  human exons 
        # get_canonical_coding_exons also sorts exons by the start in the gene
        canonical_human_exons = get_canonical_coding_exons (cursor, gene_id, ensembl_db_name['homo_sapiens'])

        # bail out if there is a problem
        if not canonical_human_exons: continue

        full_reconstituted_cDNA = ""
        prev_right_flank = ""
        for human_exon in canonical_human_exons:
            [exon_seq_id, pepseq, pepseq_transl_start, pepseq_transl_end, left_flank, right_flank, nucseq] = \
                    get_exon_seqs(cursor, human_exon.exon_id, human_exon.is_known)
            print " %10d  %s " % (human_exon.exon_id, pepseq)
            print "lengths:  %4d  %4d " % (len(pepseq)*3, len(nucseq[pepseq_transl_start:pepseq_transl_end]))
            print "pepseq_transl_start pepseq_transl_end ",  pepseq_transl_start, pepseq_transl_end, " exon len", len(nucseq)
            # add the split codon
            phase = get_exon_phase (cursor, human_exon.exon_id, human_exon.is_known)
            split_codon = ""
            if phase > 0 and prev_right_flank and left_flank:
                offset      = (3-phase)%3
                split_codon = prev_right_flank[:phase] + left_flank[-offset:]

            full_reconstituted_cDNA += split_codon + nucseq[pepseq_transl_start:pepseq_transl_end]
            prev_right_flank = right_flank

        
        if ( is_mitochondrial(cursor, gene_id)):
            full_reconstituted_seq = Seq(full_reconstituted_cDNA).translate(table="Vertebrate Mitochondrial").tostring()
        else:
            full_reconstituted_seq = Seq(full_reconstituted_cDNA).translate().tostring()
            
        canonical = get_canonical_transl (acg, cursor, gene_id, 'homo_sapiens', strip_X = False)
        if ( len(full_reconstituted_seq) != len(canonical) ):
            
            print >> logf, gene_id, stable_id, get_description (cursor, gene_id)
            print >> logf, "error reassembling,  len(full_reconstituted_seq) != len(canonical) ",  len(full_reconstituted_seq) , len(canonical) 
            print >> logf, "canonical:"
            print >> logf, canonical
            print >> logf, "reconstituted:"
            print >> logf, full_reconstituted_seq
            continue

        print gene_id, stable_id, get_description (cursor, gene_id)
        codons = map(''.join, zip(*[iter(full_reconstituted_cDNA)]*3))
        for i in range(len(codons)):
            print " %5d  %s  %s " % (i, full_reconstituted_seq[i], codons[i])

    logf.close()

#########################################
if __name__ == '__main__':
    main()

