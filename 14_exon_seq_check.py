#!/usr/bin/python

import MySQLdb
import sys, commands
from   random                 import  choice
from   el_utils.mysql         import  connect_to_mysql, search_db, switch_to_db, check_null
from   el_utils.ensembl       import  get_species, get_gene_ids, gene2stable
from   el_utils.ensembl       import  gene2exon_list, get_exon_seqs, is_mitochondrial
from   el_utils.exon          import  Exon
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.threads       import  parallelize
from   el_utils.utils         import  erropen
from   el_utils.translation   import  crop_dna, translation_bounds, translate, phase2offset

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


###########################################
def main():

    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)    
    all_species = ['homo_sapiens']

    for species in all_species:

        print species

        switch_to_db (cursor,  ensembl_db_name[species])

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        tot         = 0
        ct          = 0
        no_pepseq   = 0
        exon_seq_ok = 0
        mismatch    = 0
       
        for tot in range(500):
 
            gene_id = choice(gene_ids)
            #gene_id = 380431
            # get _all_ exons
            exons = gene2exon_list(cursor, gene_id)
            if (not exons):
                print 'no exons for gene', gene_id
                sys.exit(1)

            for exon in exons:

                # exons seqs are its aa translation, left_flank, right_flank, and dna_seq
                exon_seqs = get_exon_seqs(cursor, exon.exon_id, exon.is_known)
                if (not exon_seqs):
                    ct += 1
                    #print  'no exon seq for exon',  exon.exon_id
                    #sys.exit(1)
                    continue

                elif (exon.is_coding and exon.covering_exon < 0):

                    [exon_seq_id, protein_seq, left_flank, right_flank, dna_seq] = exon_seqs

                    if (not protein_seq or not len(protein_seq) ):
                        continue

                    exon_seq_ok += 1
                    if ( len(protein_seq)*3 >= len(dna_seq)-5):
                        continue

                    mismatch += 1
                    mitochondrial        = is_mitochondrial(cursor, gene_id)
                    [seq_start, seq_end] = translation_bounds (cursor, exon.exon_id)
                    dna_cropped          = crop_dna (seq_start, seq_end, dna_seq)
                    if ( len(protein_seq)*3 >= len(dna_cropped)-3):
                        continue

                    mismatch += 1

                    if (1):
                        [phase, pepseq]      = translate (dna_cropped, exon.phase, mitochondrial)
                     
                        print "gene id ", gene_id, gene2stable(cursor, gene_id), " exon_id ", exon.exon_id
                        print "is known: ", exon.is_known
                        print "phase:  ", exon.phase
                        print "mitochondrial: ", mitochondrial
                        print len(protein_seq)*3, len(dna_cropped), len(dna_seq)

                        print "phase suggested: ", phase
                        print " ** ", pepseq
                        offset = phase2offset(exon.phase)
                        dnaseq = Seq (dna_cropped[offset:], generic_dna)
                        print " ** ", dnaseq.translate()
 
                        print
                        print
                        exit (1)
 
        print
        print species
        print "tot number of genes checked: ", tot
        print "            without dna seq: ", ct
        print "   exons with petide seq ok: ", exon_seq_ok
        print "    dna seq length mismatch: ", mismatch
        print "   exons without petide seq: ", no_pepseq

    cursor.close()
    db    .close()



#########################################
if __name__ == '__main__':
    main()
