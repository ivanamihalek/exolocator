#!/usr/bin/python

import MySQLdb
import commands, sys
from   el_utils.mysql   import connect_to_mysql, search_db, switch_to_db, check_null
from   el_utils.mysql   import store_or_update
from   el_utils.ensembl import get_gene_ids, get_species, is_mitochondrial
from   el_utils.ensembl import gene2exon_list, get_exon_seqs, gene2stable
from   el_utils.threads import parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#########################################
def get_phase(cursor, exon_id):

    qry  = "select is_coding, phase, gene_id from gene2exon where exon_id = %d" % exon_id
    rows = search_db(cursor, qry)
    if (rows):
        [is_coding, phase, gene_id] = rows[0]
    else:
        [is_coding, phase, gene_id] = [0,0,0]

    return [is_coding, phase, gene_id]

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

    if pepseq and pepseq[-1]=='*':
        pepseq = pepseq[:-1]
    if not '*' in pepseq:
        return [phase, pepseq]

    phase = (phase+1)%3
    offset = phase2offset(phase)
    dnaseq = Seq (dna_seq[offset:], generic_dna)
    if mitochondrial:
        pepseq = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
    else:
        pepseq = dnaseq.translate().tostring()
    
    if pepseq and  pepseq[-1]=='*':
        pepseq = pepseq[:-1]
    if not '*' in pepseq:
        return [phase, pepseq]

    phase = (phase+1)%3
    offset = phase2offset(phase)
    dnaseq = Seq (dna_seq[offset:], generic_dna)
    if mitochondrial:
        pepseq = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
    else:
        pepseq = dnaseq.translate().tostring()
    
    if pepseq and  pepseq[-1]=='*':
        pepseq = pepseq[:-1]
    if not '*' in pepseq:
        return [phase, pepseq]

    return [phase, pepseq]

#########################################
def pep_exon_seqs(species_list, db_info):

    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    for species in species_list:
        
        print
        print "############################"
        print  species

        if not switch_to_db(cursor, ensembl_db_name[species]):
            return False

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')



        tot         = 0
        gene_ct     = 0
        no_pepseq   = 0
        pepseq_ok   = 0
        short_dna   = 0
        no_exon_seq = 0
        translation_fail = 0

        # for all protein coding genes in a species
        for gene_id in gene_ids:
            gene_ct += 1
            if (not  gene_ct%1000):
                print species, "tot:", gene_ct

            # for all exons in the gene
            exons = gene2exon_list(cursor, gene_id)
            if (not exons):
                print 'no exons for gene', gene_id
                sys.exit(1)

            for exon in exons:

                tot += 1
                #####################################                
                if (not exon.is_coding or  exon.covering_exon > 0):
                    continue 
                exon_seqs = get_exon_seqs(cursor, exon.exon_id, exon.is_known)
                if (not exon_seqs):
                    no_exon_seq += 1
                    continue                   
                [exon_seq_id, protein_seq, left_flank, right_flank, dna_seq] = exon_seqs
                if len(dna_seq)<3:
                    short_dna += 1
                    continue
                protein_seq = check_null (protein_seq)
                if ( protein_seq and len(protein_seq)>0):
                    pepseq_ok += 1
                    continue

                #####################################                
                mitochondrial        = is_mitochondrial(cursor, gene_id)
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
                [phase, pepseq] = translate (dna_seq, exon.phase, mitochondrial)

                if ( not pepseq): # usually some short pieces (end in pos 4 and such)
                    translation_fail += 1
                    continue
 
                if pepseq[-1] == '*':
                    pepseq = pepseq[:-1]

                if  '*' in pepseq: #store
                    translation_fail += 1
                    continue

                qry  = "update exon_seq "
                qry += "set protein_seq   = '%s' " % pepseq
                qry += "where exon_seq_id =  %d  " % exon_seq_id
                rows = search_db (cursor, qry)
                if (rows):
                    rows = search_db (cursor, qry, verbose = True)
                    exit(1)
             
                pepseq_ok += 1
           

        print species
        print "total coding exons ", tot
        print "no exon seq info   ", no_exon_seq
        print "short dna          ", short_dna
        print "transl failure     ", translation_fail
        print "pepseq ok          ", pepseq_ok
            

#########################################
def main():
    no_threads = 5

    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    all_species = ['equus_caballus', 'echinops_telfairi', 'oryctolagus_cuniculus', 
                   'ornithorhynchus_anatinus', 'oryzias_latipes']
    parallelize (no_threads, pep_exon_seqs, all_species, [local_db, ensembl_db_name] )



#########################################
if __name__ == '__main__':
    main()


'''
               if seq_start: print "start ", seq_start, 
                    if seq_end:   print "end ", seq_end 
                    print "exon_id ", exon_id
                    print "phase ", phase 
                    print "orig dna ", dna_seq_orig
                    print "orig", pepseq_orig
                    print "    ", pepseq
                    exit(1)
                    continue
'''
