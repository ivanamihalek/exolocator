
from   el_utils.mysql   import *
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

#########################################
def phase2offset(phase):
    if phase > 2:
        phase = phase%3
    if phase==0:
        offset = 0
    else:
        offset = 3-phase
    return offset

##############################################
def crop_dna (seq_start, seq_end, dna_seq):

    dna_cropped = dna_seq

    seq_start   = check_null(seq_start)
    seq_end     = check_null(seq_end)

    if seq_start is None and  seq_end is None:
        pass
    else:
        if seq_start is None:
            dna_cropped = dna_seq[:seq_end]
        elif seq_end is None:
            if seq_start>0:
                dna_cropped = dna_seq[seq_start-1:]
        else:
            if seq_start>0:
                dna_cropped = dna_seq[seq_start-1:seq_end]
            else:
                dna_cropped = dna_seq[:seq_end]
 

    return dna_cropped

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

