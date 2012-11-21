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
# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


#########################################
def translation_bounds(cursor, exon_id):

    seq_start = None
    seq_end   = None

    qry  = "select seq_start, start_exon_id"
    qry += " from translation "
    qry += " where start_exon_id = %d"  % exon_id
    rows = search_db(cursor, qry, verbose=False)
    if rows:
        [seq_start, start_exon_id] = rows[0]

    qry  = "select seq_end, end_exon_id"
    qry += " from translation "
    qry += " where end_exon_id = %d"  % exon_id
    rows = search_db(cursor, qry, verbose=False)
    if rows:
        [seq_end, end_exon_id] = rows[0]

    return [seq_start, seq_end]

#########################################
def phase2offset(phase):
    if phase > 2:
        phase = phase%3
    if phase==0:
        offset = 0
    else:
        offset = 3-phase
    return offset
#########################################
def  translate (dna_seq, phase, mitochondrial=False):
    pepseq = ""
    if phase < 0: phase = 0
    offset = phase2offset(phase)
    dnaseq = Seq (dna_seq[offset:], generic_dna)
    if mitochondrial:
        pepseq = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
    else:
        pepseq = dnaseq.translate().tostring()

    if  pepseq and pepseq[-1]=='*':
        pepseq = pepseq[:-1]
    print "orig phase: ", phase, offset, pepseq

    if not '*' in pepseq:
        return pepseq

    phase = (phase+1)%3
    offset = phase2offset(phase)
    dnaseq = Seq (dna_seq[offset:], generic_dna)
    if mitochondrial:
        pepseq = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
    else:
        pepseq = dnaseq.translate().tostring()
    if  pepseq and pepseq[-1]=='*':
        pepseq = pepseq[:-1]
    
    print "phase+1: ",  phase, offset, pepseq

    if not '*' in pepseq:
        return pepseq

    phase = (phase+1)%3
    offset = phase2offset(phase)
    dnaseq = Seq (dna_seq[offset:], generic_dna)
    if mitochondrial:
        pepseq = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
    else:
        pepseq = dnaseq.translate().tostring()
    if  pepseq and pepseq[-1]=='*':
        pepseq = pepseq[:-1]
    
    print "phase+2: ", phase,  offset, pepseq

    if not '*' in pepseq:
        return pepseq

    return pepseq
#########################################
def main():

    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)    
    all_species = ['tupaia_belangeri', 'tetraodon_nigroviridis']

    for species in all_species:

        switch_to_db (cursor,  ensembl_db_name[species])

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        tot         = 0
        ct          = 0
        no_pepseq   = 0
        exon_seq_ok = 0
       
        for tot in range(500):
 
            gene_id = choice(gene_ids)

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
                    
                    if ( exon_seqs[0]):
                        exon_seq_ok += 1
                    else:

                        [protein_seq, left_flank, right_flank, dna_seq] = exon_seqs
                        mitochondrial =  is_mitochondrial(cursor, gene_id)
                        print "no translation "
                        print "gene id ", gene_id, gene2stable(cursor, gene_id), " exon_id ", exon.exon_id
                        print "mitochondrial: ", mitochondrial

                        # check if there is annotation about translation starting
                        # or ending in this exon
                        [seq_start, seq_end] = translation_bounds (cursor, exon.exon_id)
                        seq_start = check_null(seq_start)
                        seq_end   = check_null(seq_end)

                        if seq_start is None and  seq_end is None:
                            pass
                        else:
                            if seq_start is None:
                                dna_seq = dna_seq[:seq_end]
                            elif seq_end is None:
                                if seq_start>0:
                                    dna_seq = dna_seq[seq_start-1:]
                            else:
                                if seq_start>0:
                                    dna_seq = dna_seq[seq_start-1:seq_end]
                                else:
                                    dna_seq = dna_seq[:seq_end]

                        pepseq = translate (dna_seq, exon.phase, mitochondrial)

                        print " ** ", pepseq
                        print
                        print
                        no_pepseq += 1

 
        print
        print species
        print "tot number of genes checked: ", tot
        print "            without dna seq: ", ct
        print "   exons with petide seq ok: ", exon_seq_ok
        print "   exons without petide seq: ", no_pepseq

    cursor.close()
    db    .close()



#########################################
if __name__ == '__main__':
    main()
