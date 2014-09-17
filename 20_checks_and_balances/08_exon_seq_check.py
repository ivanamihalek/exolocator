#!/usr/bin/python -u

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
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

###########################################
def test1 ():
    mismatch += 1
    mitochondrial        = is_mitochondrial(cursor, gene_id)
    [seq_start, seq_end] = translation_bounds (cursor, exon.exon_id)
    dna_cropped          = crop_dna (seq_start, seq_end, dna_seq)
    if ( len(protein_seq)*3 >= len(dna_cropped)-3):
        return ""

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
 
###########################################
def main():

    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)    

    for species in all_species:

        if not species=='homo_sapiens': continue

        print
        print species

        switch_to_db (cursor,  ensembl_db_name[species])

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        tot_exons   = 0
        no_exon_seq = 0
        short_dna   = 0
        pepseq_ok   = 0
        mismatch    = 0
        stored_incorrect = 0
        translation_fail = 0
        #####################################
        for gene_id in [10092907]:
        #for tot in range(1000):
            #gene_id = choice(gene_ids)

            # get _all_ exons
            exons = gene2exon_list(cursor, gene_id)
            if (not exons):
                print 'no exons for gene', gene_id
                sys.exit(1)

            for exon in exons:

                #####################################                
                if not exon.is_coding:
                    print exon.exon_id, " not coding "
                    continue
                if exon.covering_exon >0:
                    print exon.exon_id, " is covered by ", exon.covering_exon 
                    continue
                    

                tot_exons += 1
                # exons seqs are its aa translation, left_flank, right_flank, and dna_seq
                exon_seqs = get_exon_seqs(cursor, exon.exon_id, exon.is_known)
                if (not exon_seqs):
                    no_exon_seq += 1
                    print "no exon seqs for  ", gene_id, exon.exon_id
                    #exit(1)
                    continue                   

                [exon_seq_id, pepseq, pepseq_transl_start, 
                 pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs

                if len(dna_seq)<3:
                    short_dna += 1
                    print "short_dna:", dna_seq
                    continue

                if (pepseq_transl_start == -10): # ??? what is this shit? adn what happens downstream if the pepseq_transl_start is None?
                    translation_fail += 1
                    print "pepseq_transl_start:", pepseq_transl_start
                    continue

                mitochondrial        = is_mitochondrial(cursor, gene_id)
                dnaseq  = Seq (dna_seq[pepseq_transl_start:pepseq_transl_end], generic_dna)
                if (mitochondrial):
                    pepseq2 = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
                else:
                    pepseq2 = dnaseq.translate().tostring()

                if True:
                    print exon.exon_id
                    print "pep stored:", pepseq
                    print "dna transl:", pepseq2
                    print "dna begin:", dna_seq[:12]
                    print "start:" , pepseq_transl_start, 
                    print "end:",  pepseq_transl_end
                    print

                if (not pepseq == pepseq2):
                    stored_incorrect += 1
                else:
                    pepseq_ok += 1

        print "total coding exons ", tot_exons
        print "no exon seq info   ", no_exon_seq
        print "short dna          ", short_dna
        print "transl failure     ", translation_fail
        print "stored pepseq does not correspond to the translation of stored dna:   ", stored_incorrect
        print "pepseq ok          ", pepseq_ok

    cursor.close()
    db    .close()



#########################################
if __name__ == '__main__':
    main()
