#!/usr/bin/python

import MySQLdb
import commands
from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.mysql   import  store_or_update
from   el_utils.ensembl import get_species
from   el_utils.threads import parallelize
# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#########################################
def get_phase(cursor, exon_id):

    qry  = "select is_coding, phase from gene2exon where exon_id = %d" % exon_id
    rows = search_db(cursor, qry)
    if (rows):
        [is_coding, phase] = rows[0]
    else:
        [is_coding, phase] = [0,0]

    return [is_coding, phase]
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
    if phase==0:
        offset = 0
    else:
        offset = phase+1
    return offset

#########################################
def  translate (dna_seq, phase):
    pepseq = ""
    if phase < 0: phase = 0
    offset = phase2offset(phase)
    dnaseq = Seq (dna_seq[offset:], generic_dna)
    pepseq = dnaseq.translate().tostring()

    if not '*' in pepseq:
        return pepseq

    phase = (phase+1)%3
    offset = phase2offset(phase)
    dnaseq = Seq (dna_seq[offset:], generic_dna)
    pepseq = dnaseq.translate().tostring()
    
    if not '*' in pepseq:
        return pepseq

    phase = (phase+1)%3
    offset = phase2offset(phase)
    dnaseq = Seq (dna_seq[offset:], generic_dna)
    pepseq = dnaseq.translate().tostring()
    
    if not '*' in pepseq:
        return pepseq

    return pepseq

#########################################
def pep_exon_seqs(species_list, ensembl_db_name):

    db     = connect_to_mysql()
    cursor = db.cursor()

    for species in species_list:
        
        if (species == 'homo_sapiens'):
            continue
        print
        print "############################"
        print  species

        if not switch_to_db(cursor, ensembl_db_name[species]):
            return False

        range_end = 0

        while True:
            range_start = range_end   +   1
            range_end   = range_start + 999

            qry  = "select exon_seq_id, exon_id, dna_seq, protein_seq "
            qry += " from exon_seq  where exon_seq_id >= %d " % range_start
            qry += " and exon_seq_id <= %d "  % range_end 
            rows = search_db(cursor, qry, verbose=False)
            if not rows:
                break
            ct   = 0
            fail = 0
            for row in rows:
                [exon_seq_id, exon_id, dna_seq, protein_seq] = row
                if protein_seq:
                    continue
                if len(dna_seq)<3:
                    continue
                [is_coding, phase]   = get_phase (cursor, exon_id)
                if not is_coding:
                    continue
                ct   += 1

                [seq_start, seq_end] = translation_bounds (cursor, exon_id)

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
   
                            
                pepseq = translate (dna_seq, phase)
 
                if ( not pepseq): # ususally some short pieces (end in pos 4 and such)
                    fail +=1
                    continue
 
                if pepseq[-1] == '*':
                    pepseq = pepseq[:-1]

                if  not '*' in pepseq: #store
                    qry  = "update exon_seq "
                    qry += "set protein_seq='%s' "   % protein_seq
                    qry += "where exon_seq_id = %d " % exon_seq_id
                    rows = search_db (cursor, qry)
                    if (rows):
                        rows = search_db (cursor, qry, verbose = True)
                        exit(1)
                else: # check if there is annotation about translation starting
                    # or eding in this exon - if not store with the stop codon (?)
                    fail +=1
                    #if seq_start: print "start ", seq_start, 
                    #if seq_end:   print "end ", seq_end, 
                    #print "exon_id ", exon_id,
                    #print "phase ", phase 
                    #print "new pepseq ", pepseq
                    pass

            print species, range_start, range_end, "\n\t", ct, fail

#########################################
def main():

    no_threads = 5

    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    parallelize (no_threads, pep_exon_seqs, all_species, ensembl_db_name)



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
