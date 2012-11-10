#!/usr/bin/python

import MySQLdb
import commands
from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.mysql   import  store_or_update
from   el_utils.ensembl import  get_species, get_gene_ids, gene2exon_list
from   el_utils.ensembl import  gene2stable, gene2stable_canon_transl, is_mitochondrial
from   el_utils.ensembl import  gene2canon_transl, get_canonical_exons, get_selenocysteines
from   el_utils.exon    import  Exon
from   el_utils.threads import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#def is_mitochondrial (cursor, seq_region_id):

    


#########################################
def  get_primary_seq_info (cursor, gene_id, species):

    # seq identifier from gene table
    qry  = "select seq_region_id, seq_region_start, seq_region_end, seq_region_strand"
    qry += " from gene where gene_id = %d" % gene_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         exit(1)
    [seq_region_id, seq_region_start, seq_region_end, seq_region_strand] = rows[0]
    

    qry  = "select name, file_name from seq_region "
    qry += " where seq_region_id= %d" %  seq_region_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         return []
    [seq_name, file_names] = rows[0]
    # Ivana:
    # What is the proper way to find out whether the seq_region is mitochondrial?
    # EMily:
    # It's in the seq_region table, under name. Name will either be a chromosome
    # number, MT for mitochrondria or the name of a contig.
    mitochondrial = is_mitochondrial(cursor, gene_id)

    return [seq_name, file_names, seq_region_start, 
            seq_region_end, seq_region_strand, mitochondrial]

#########################################
def  get_alt_seq_info (cursor, gene_id, species):

    # seq identifier from gene table
    qry  = "select seq_region_id, seq_region_strand from gene where gene_id = %d" % gene_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         exit(1)
    [seq_region_id, seq_region_strand] = rows[0]
    
    # check whether we have "assembly exception"
    # we do not want 'PAR' regions, though:
    '''
    The pseudo-autosomal regions are homologous DNA sequences on the (human) X and Y chromosomes. 
    They allow the pairing and crossing-over of these sex chromosomes the same way the autosomal 
    chromosomes do during meiosis. 
    As these genomic regions are identical between X and Y, they are oftentimes only stored once.
    '''

    qry  = "select seq_region.name,  assembly_exception.exc_seq_region_start, assembly_exception.exc_seq_region_end "
    qry += "from seq_region, assembly_exception "
    qry += "where seq_region.seq_region_id = assembly_exception.exc_seq_region_id "
    qry += "and assembly_exception.seq_region_id = %d" % seq_region_id
    qry += " and not assembly_exception.exc_type = 'PAR'"
    rows = search_db (cursor, qry)
    if (rows):
        [seq_name, seq_region_start, seq_region_end] = rows[0]
        qry = " select distinct file_name from seq_region where seq_region.name = '%s' " % seq_name
        rows = search_db (cursor, qry)
        file_names = ""
        for row in rows:
            if file_names:
                file_names += " "
            file_names += row[0]

        mitochondrial = is_mitochondrial (cursor, gene_id)
        return [seq_name, file_names, seq_region_start, 
                seq_region_end, seq_region_strand, mitochondrial]
    else:
        return []

#########################################
def extract_gene_seq (acg, species, seq_name, file_names, seq_region_strand,  
                      seq_region_start, seq_region_end):

    # now the question is, which file do I use if there are several options?
    first_choice  = ""
    second_choice = ""
    for file_name in file_names.split(" "):
        if  '.chromosome.' in  file_name:
            first_choice = file_name
        elif '.toplevel.' in  file_name:
            second_choice = file_name
            
    fasta_db_file = ""
    if first_choice:
        fasta_db_file = first_choice
    elif second_choice:
        fasta_db_file = second_choice

    if not fasta_db_file:
        print "failed to decide on fasta_db_file:"
        print file_names
        exit(1)


    # extract gene sequence  from fasta db
    fastacmd = acg.generate_fastacmd_gene_command(species, seq_name, fasta_db_file,
                                                  seq_region_strand,  seq_region_start,    
                                                  seq_region_end)
    ret = commands.getoutput(fastacmd)
    if not ret:
        print "no refturn for fastacmd for", species, gene_id
        print "fastacmd: ", fastacmd
        exit (1)
    if ('ERROR' in ret):
        print "Error running fastacmd: ", fastacmd
        print ret
        if 'Ignoring sequence location' in ret:
            print 'will ignore'
            print
        else:
            exit (1)
    gene_seq = ""
    reading = 0
    for line in ret.split("\n"):
        if ('>' in line):
            reading = 1
            continue
        if (not reading):
            continue
        line.rstrip()
        gene_seq += line

    return gene_seq
#########################################
def strip_stop(pepseq):
    if (not pepseq or len(pepseq)==0):
        return pepseq
    if ( pepseq[-1] == '*'):
        pepseq = pepseq[:-1]
    return pepseq

#########################################
def canonical_transl_info (cursor, gene_id):
    
    qry  = "select canonical_transcript_id from gene  where gene_id = %d " % gene_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         return []
    
    canonical_transcript_id = int(rows[0][0])
    qry = "select start_exon_id, seq_start, end_exon_id,  seq_end "
    qry += " from translation where transcript_id = %d " % canonical_transcript_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         return []
 
    return rows[0]
#########################################
def get_canonical_transl (acg, cursor, gene_id, species):

    canonical_translation = ""

    canonical_transl_id = gene2stable_canon_transl(cursor, gene_id)
    if ( not canonical_transl_id):
        print "no canonical transl id found for ", gene_id
        exit(1)

    cmd = acg.generate_fastacmd_protein_command (canonical_transl_id, species, 
                                                 "all", None)
    fasta = commands.getoutput(cmd)
    if (not fasta):
        print gene2stable (cursor, gene_id = gene_id), 
        print "fasta not found for ", canonical_transl_id
        exit(1)

    canonical_translation = ""
    for line in fasta.split("\n"):
        if ('>' in line):
            continue
        line.rstrip()
        canonical_translation += line

    while (len(canonical_translation) and canonical_translation[0] == 'X'):
        canonical_translation = canonical_translation[1:]

    return canonical_translation

########################################
def  get_canonical_exons (cursor, gene_id):

    exons = gene2exon_list (cursor, gene_id)
    if (not exons):
        print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
        exit(1) # shouldn't happen at this point

    # sorting exons in place by their start in gene:
    exons.sort(key=lambda exon: exon.start_in_gene)

    canonical_coding_exons = []
    reading = False
    for exon in exons:
        if (not exon.is_canonical): 
             continue
        if (not exon.canon_transl_start is None):
            reading = True
        if (reading):
            canonical_coding_exons.append(exon)
        if (not exon.canon_transl_end is None):
            break  

    return canonical_coding_exons


#########################################
def  transl_reconstruct (cursor,  gene_id, gene_seq, canonical_coding_exons, 
                         is_mitochondrial, verbose = False):

    canonical_exon_pepseq = {}
    translated_seq = "" 

    [can_transl_start_exon, can_transl_start_position,
     can_transl_end_exon, can_transl_end_position] = canonical_transl_info (cursor, gene_id)


    # do we have any selenocysteines by any chance
    selenoC_pos = get_selenocysteines (cursor,  gene_id)

    carry = ""
    ok_so_far = True
    # sanity checking
    for exon in canonical_coding_exons:

        #print
        #print "exon", exon.exon_id
        #find exon sequence within the gene
        start = exon.start_in_gene
        if (exon is canonical_coding_exons[0]):
            if ( not exon.exon_id == can_transl_start_exon ):
                print " error start cantransl:  gene_id ",  gene_id,
                print  " exon_id ", exon.exon_id, " canon: ", can_transl_start_exon
                exit (1)
            start +=  exon.canon_transl_start

        if ( exon is canonical_coding_exons[-1]):
            if ( not exon.exon_id == can_transl_end_exon ):
                print " error end cantransl:  gene_id ",  gene_id,
                print  " exon_id ", exon.exon_id, " canon: ", can_transl_end_exon
                exit (1)
            end = exon.start_in_gene + exon.canon_transl_end
        else:
            end = exon.end_in_gene

        if (not exon.phase == -1 and not exon.phase == len(carry)):
            #print "Houston we have a problem: exon phase =", exon.phase,
            #print " the length of carry =", len(carry), 
            #print " (gene_id %d, exon_id %d) " % (gene_id, exon.exon_id)
            if (  exon.phase ):
                start += 3-exon.phase
            carry = ""


        exon_seq     =  gene_seq[ start: end+1]
        exon_seq_for_transl_purposes = carry + exon_seq

        remainder    = len(exon_seq_for_transl_purposes)%3
        if ( remainder == 0 ):
            carry = ""
        elif (remainder == 1 ):
            carry    = exon_seq_for_transl_purposes[-1:]
            exon_seq_for_transl_purposes  = exon_seq_for_transl_purposes[:-1]
        else:
            carry    = exon_seq_for_transl_purposes[-2:]
            exon_seq_for_transl_purposes = exon_seq_for_transl_purposes[:-2]

        dnaseq = Seq (exon_seq_for_transl_purposes, generic_dna)
        pepseq = dnaseq.translate()

        # turn to the corresponding BioPython object
        if ( is_mitochondrial ):
            pepseq = dnaseq.translate(table="Vertebrate Mitochondrial")
        else:
            pepseq = dnaseq.translate()
            
        # replace stop codons from selenoC positions, if there are such
        if (selenoC_pos):
            for pos_in_full_length_translation in selenoC_pos:
                pos_in_exon_translation = pos_in_full_length_translation-len(translated_seq)
                if pos_in_exon_translation<0 or pos_in_exon_translation>len(pepseq):
                    continue
                tempseq = pepseq.tomutable()
                tempseq[pos_in_exon_translation] = 'U'
                pepseq  = tempseq.toseq()
 
        # strip the last stop codon only
        if ( exon is canonical_coding_exons[-1]):
            pepseq = strip_stop(pepseq)  

        pepseq0 =  pepseq
        if verbose:
            print "phase 0", pepseq0

        # if there are stil stop codons we'll give another shot 
        # to the possibility that it is mitochondrial, (and we ddin't know it)
        # after that we cry foul
        if ( '*' in pepseq):
            pepseq = dnaseq.translate(table="Vertebrate Mitochondrial")

        # strip the last stop codon only
        if ( exon is canonical_coding_exons[-1]):
            pepseq = strip_stop(pepseq)  

        # some further  desperate measures 
        ok_so_far = True
        if ( '*' in pepseq):
            ok_so_far = False
            dnaseq = Seq (exon_seq_for_transl_purposes[1:], generic_dna)
            pepseq = dnaseq.translate()
            pepseq = strip_stop(pepseq)
            pepseq1 =  pepseq
            if verbose:
                print "phase 1", pepseq1
        else:
            ok_so_far = True

        if (not ok_so_far and '*' in pepseq):
            dnaseq = Seq (exon_seq_for_transl_purposes[2:], generic_dna)
            pepseq = dnaseq.translate()
            pepseq = strip_stop(pepseq)
            pepseq2 =  pepseq
            if verbose:
                print "phase 2", pepseq2
        else:
            ok_so_far = True

        if (not ok_so_far and  '*' in pepseq):
            if verbose:
                print "Error: stop codon "
        else:
            ok_so_far = True

        if ( not ok_so_far):
            return [{}, ""]
        translated_seq += pepseq # I need the seq in selenoC - to decide
                                 # where the position of  U should be
        canonical_exon_pepseq[exon.exon_id] = pepseq.tostring()

    return [canonical_exon_pepseq,translated_seq] 

#########################################
def compare_seqs (canonical_translation, translated_seq, verbose=False):

    comparison_ok = True

    while (len(translated_seq) and translated_seq[0] == 'X'):
        translated_seq = translated_seq[1:]

    difference = len(translated_seq) - len(canonical_translation)
    if ( abs(difference) > 3):
        comparison_ok = False
        if verbose:
            print
            print ">canon"
            print canonical_translation
            print ">exons"
            print translated_seq
            print
    else:
        diff  =  0
        start = -1
        for i in range(len(translated_seq)):
            if ( i >= len(canonical_translation)):
                break
            if (not translated_seq[i] ==  canonical_translation[i]):
                diff += 1
                if start < 0:
                    start = i
        if (diff > 2):
            comparison_ok = False
            if verbose:
                print
                print ">canon"
                print canonical_translation
                print ">exons"
                print translated_seq
                print translated_seq[start], canonical_translation[start]
                print "nuber of  diff sites: ", diff, " starting from ", start
                print

    return comparison_ok

#########################################
def  get_gene_seq (acg, cursor, gene_id, species):

    null = ["",{}]

    #########################################
    # which file should we be looking in, which sequence, from where to where
    ret = get_primary_seq_info (cursor, gene_id, species)
    if (not ret):
        return null
    [seq_name, file_names, seq_region_start, seq_region_end, 
     seq_region_strand, is_mitochondrial] = ret
    # i'm not quite clear why Ensembl is doing this, but sometimes we need the alternative
    # region - ("PATCH" deposited as tte right sequence, but its missing most of the gene)
    # so first establish that is it the case: find canonical translation
    canonical_translation  = get_canonical_transl (acg, cursor, gene_id, species)
    # find all canonical exons associated with the gene id
    canonical_coding_exons = get_canonical_exons (cursor, gene_id)
    # extract raw gene  region
    gene_seq = extract_gene_seq(acg, species, seq_name, file_names, seq_region_strand,  
                                seq_region_start, seq_region_end)
    # reconstruct the translation from the raw gene_seq and exon boundaries
    [canonical_exon_pepseq,translated_seq] = transl_reconstruct (cursor, gene_id, gene_seq, canonical_coding_exons, 
                                                 is_mitochondrial)
    if (translated_seq):
        # compare the two sequences and cry foul if they are not the same:
        comparison_ok = compare_seqs (canonical_translation, translated_seq)
    else:
        comparison_ok = False
    # if we succefully translated the exons, and came up with the same answer 
    # as the canonical translation, we are done here
    if (comparison_ok):
        return [gene_seq, canonical_exon_pepseq]
 
    #########################################
    # otherwise repeat the procedure with the alternative seq info:
    ret = get_alt_seq_info (cursor, gene_id, species)
    if (not ret):
        return null
    [seq_name, file_names, seq_region_start, seq_region_end, 
     seq_region_strand, is_mitochondrial] = ret
      # i'm not quite clear why Ensembl is doing this, but sometimes we need the alternative
    # region - ("PATCH" deposited as tte right sequence, but its missing most of the gene)
    # so first establish that is it the case: find canonical translation
    canonical_translation  = get_canonical_transl (acg, cursor, gene_id, species)
    # find all canonical exons associated with the gene id
    canonical_coding_exons = get_canonical_exons (cursor, gene_id)
    # extract raw gene  region
    gene_seq = extract_gene_seq(acg, species, seq_name, file_names, seq_region_strand,  
                                seq_region_start, seq_region_end)
    # reconstruct the translation from the raw gene_seq and exon boundaries
    [canonical_exon_pepseq,translated_seq] = transl_reconstruct (cursor, gene_id, gene_seq, canonical_coding_exons, 
                                                                  is_mitochondrial)
    if (translated_seq):
        # compare the two sequences and cry foul if they are not the same:
        comparison_ok = compare_seqs (canonical_translation, translated_seq)
    else:
        comparison_ok = False
    # if we succefully translated the exons, and came up with the same answer 
    # as the canonical translation, we are done here
    if (comparison_ok):
        return [gene_seq, canonical_exon_pepseq]
 
    return null 

#########################################
def get_exon_seqs (gene_seq, exons):

    exon_seq    = {} 
    left_flank  = {}
    right_flank = {}

    for exon in exons:
        start   = exon.start_in_gene
        end     = exon.end_in_gene
        exon_id = exon.exon_id
        left_flank[exon_id]  = gene_seq[start-15:start]
        exon_seq[exon_id]    = gene_seq[start:end+1]
        right_flank[exon_id] = gene_seq[end+1:end+16]
        
    return [exon_seq, left_flank, right_flank]
            
#########################################
def store (cursor, exons, exon_seq, left_flank, right_flank, canonical_exon_pepseq):
    
    for exon in exons:
        exon_id = exon.exon_id
        #####
        fixed_fields  = {}
        update_fields = {}
        fixed_fields['exon_id']          = exon_id
        fixed_fields['is_known']         = exon.is_known
        fixed_fields['is_sw']            = 0
        if (exon_seq[exon_id]):
            update_fields['dna_seq']     = exon_seq[exon_id]
        if (left_flank[exon_id] ):
            update_fields['left_flank']  = left_flank[exon_id]
        if (right_flank[exon_id] ):
            update_fields['right_flank'] = right_flank[exon_id]
        if ( canonical_exon_pepseq.has_key(exon_id) and canonical_exon_pepseq[exon_id]):
            update_fields['protein_seq'] = canonical_exon_pepseq[exon_id]
        #####
        store_or_update (cursor, 'exon_seq', fixed_fields, update_fields)



#########################################
def store_exon_seqs(species_list, ensembl_db_name):

    db     = connect_to_mysql()
    cursor = db.cursor()
    acg     = AlignmentCommandGenerator()

    for species in species_list:
        if (species == 'homo_sapiens'):
            continue
        print
        print "############################"
        print  species

        if not switch_to_db(cursor, ensembl_db_name[species]):
            return False      
 
        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        ###########################
        seqs_not_found = []
        ct  = 0
        tot = 0
        for gene_id in gene_ids:
            tot += 1
            if (not  tot%1000):
                print species, ct, tot
                
            # extract raw gene  region - bonus return from checking whether the 
            # sequence is correct: translation of canonical exons
            [gene_seq, canonical_exon_pepseq] = get_gene_seq(acg, cursor, gene_id, species)
            if (not gene_seq):
                ct += 1
                print 'no sequence found for ', gene_id, "   ",   ct, "out of ", tot
                seqs_not_found.append(gene_id)
                continue
            # get _all_ exons
            exons = gene2exon_list(cursor, gene_id, ensembl_db_name[species])
            if (not exons):
                print 'no exons for ', gene_id
                exit(1)

            # get the sequence for each of the exons, as well as for the flanks
            # (the return are three dictionaries, with exon_ids as keys)
            [exon_seq, left_flank, right_flank] = get_exon_seqs (gene_seq, exons)
            # store (exons, dna, protein)
            store (cursor, exons, exon_seq, left_flank, right_flank, canonical_exon_pepseq)

        print species, ct, tot
        if (seqs_not_found):
            outf = open(species+".seqs_not_found", "w")
            for not_found in seqs_not_found:
                print >> outf, str(not_found)+", ",
            print
            print
            outf.close
    cursor.close()
    db    .close()



#########################################
def main():

    no_threads = 5

    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    parallelize (no_threads, store_exon_seqs, all_species, ensembl_db_name)





#########################################
if __name__ == '__main__':
    main()
'''
    # Ivana:
    # What is the proper way to find out whether the seq_region is mitochondrial?
    # EMily:
    # It's in the seq_region table, under name. Name will either be a chromosome
    # number, MT for mitochrondria or the name of a contig.
Hi Ivana

To find which contigs are mitochondrial, use the assembly table.

select s1.name from seq_region s1, assembly asm, seq_region s2 where
s1.seq_region_id = asm.cmp_seq_region_id and s2.seq_region_id =
asm.asm_seq_region_id and s2.name = "MT" ;
in human, returns
+----------------+
| name |
+----------------+
| NC_012920 |
| J01415.1.16569 |
+----------------+

Emily

Hi Ivana

The reason for this is that patches are not always the same length as the
genomic region they patch over. In most cases, a patch corrects sequencing
errors but the number of bp stays the same, but in some cases, a patch adds a
chunk of sequence. We anchor the 5' (wrt the chromosome orientation) end of the
patch to the identical reference coordinates, and allow the 3' end to differ
slightly from the reference coordinates. We don't change the complete genomic
coordinates every time a patch is added (this happens when we bring out a new
assembly) as this would be more hassle than it's worth and most people don't
notice it anyway. However the one thing it does affect is the coordinates of
genes that overlap the 3' end of the patch.

ENSG00000261899, and presumably the other genes you've had a similar issue
with, overlaps the 3' end of a patch, so its coordinates on the patch are
shifted compared to its coordinates on the reference genome.

Here's the data from the assembly_exception table for the particular patch that
affects ENSG00000261899:
+-----------------------+---------------+------------------+----------------+-------------+------------------+-----+
| assembly_exception_id | seq_region_id | seq_region_start | seq_region_end |exc_type | exc_seq_region_id | exc_seq_region_start | exc_seq_region_end | ori
|
+-----------------------+---------------+------------------+----------------+-------------+-------------------+----------------------+--------------------+-----+
| 67                    | 1000057054    | 36453102 | 36596491 | PATCH_NOVEL | 27508 | 36453102 |36590458 | 1 |
+-----------------------+---------------+------------------+----------------+-------------+-------------------+----------------------+--------------------+-----+

To get over this you need to link over the assembly_exception table. This will
give you the exact relationship between the patch you are looking at and the
chromosome it is linked to. If you then use the
seq_region_start/exc_seq_region_start relationship, this will cover all cases,
whether there is a shift or not.

Hope this helps,

Emily
'''
