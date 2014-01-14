#!/usr/bin/python

import MySQLdb
import commands
from random import choice
from   el_utils.mysql       import  connect_to_mysql, search_db
from   el_utils.ensembl     import  *
from   el_utils.el_specific import  *
from   el_utils.exon        import  Exon
from   el_utils.threads     import  parallelize
from   el_utils.almt_cmd_generator  import AlignmentCommandGenerator
from   el_utils.special_gene_sets   import  get_theme_ids
from   el_utils.config_reader       import ConfigurationReader

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna



#########################################
def  get_seq_info (cursor, gene_id, species):

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

    return [seq_name, file_names, seq_region_start, seq_region_end, seq_region_strand]

#########################################
def  get_alt_seq_info (cursor, gene_id, species):

    # seq identifier from gene table
    qry  = "select seq_region_id from gene where gene_id = %d" % gene_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         exit(1)
    seq_region_id = rows[0][0]
    
    # Check whether we have "assembly exception."
    # The region in question should not be 'PAR' - that's something else:
    '''
    The pseudo-autosomal regions are homologous DNA sequences on the (human) X and Y chromosomes. 
    They allow the pairing and crossing-over of these sex chromosomes the same way the autosomal 
    chromosomes do during meiosis. 
    As these genomic regions are identical between X and Y, they are oftentimes only stored once.
    '''

    qry  = "select seq_region.name,  assembly_exception.exc_seq_region_start, assembly_exception.exc_seq_region_end "
    qry += "from seq_region, assembly_exception "
    qry += "where seq_region.seq_region_id = assembly_exception.exc_seq_region_id "
    qry += "and assembly_exception.seq_region_id = %d " % seq_region_id
    qry += "and not assembly_exception.exc_type = 'PAR'"

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
        return [seq_name, file_names, seq_region_start, seq_region_end]
    else:
        return []

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
def get_gene_seq (acg, species, seq_name, file_names, seq_region_strand,  seq_region_start, seq_region_end):

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
def  transl_reconstruct (cursor,  gene_id, gene_seq, canonical_coding_exons, verbose = False):

    translated_seq = "" 

    [can_transl_start_exon, can_transl_start_position,
     can_transl_end_exon, can_transl_end_position] = canonical_transl_info (cursor, gene_id)

    is_MT = is_mitochondrial (cursor, gene_id)

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
        if ( is_MT ):
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
            return ""
        translated_seq += pepseq

    return translated_seq

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
def check_canonical_sequence(local_db, species_list, ensembl_db_name):

    verbose = False

    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
        cfg    = ConfigurationReader()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    species_list = ['homo_sapiens']
    #for species in species_list:
        print
        print "############################"
        print  species

        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        
        #gene_ids = get_theme_ids(cursor, ensembl_db_name, cfg, 'missing_seq')
        
        ct  = 0
        tot = 0
        
        seen = {}
        #for gene_id in gene_ids:
        #for gene_id in [412667]:
        for tot in range(1000):
            gene_id = choice(gene_ids)
            #tot +=1 
 
            if seen.has_key(gene_id):
                continue
            seen[gene_id] = True
            #print 
            #print gene_id, gene2stable (cursor, gene_id), get_description (cursor, gene_id)
            # find canonical translation
            canonical_translation  = get_canonical_transl (acg, cursor, gene_id, species)

            # find all canonical exons associated with the gene id
            canonical_coding_exons = get_canonical_exons (cursor, gene_id)

            # find seq_name and region
            ret = get_seq_info (cursor, gene_id, species)
            if ret:
                [seq_name, file_names,  seq_region_start, 
                 seq_region_end, seq_region_strand] = ret
            else:
                print "no seq info "
                ct += 1 
                continue


            # extract raw gene  region
            gene_seq = get_gene_seq( acg, species, seq_name, file_names, seq_region_strand,  
                                     seq_region_start, seq_region_end)

            # reconstruct the translation from the raw gene_seq and exon boundaries
            translated_seq = transl_reconstruct (cursor, gene_id, gene_seq, 
                                                 canonical_coding_exons, verbose = verbose)
            if (translated_seq):
                # compare the two sequences and cry foul if they are not the same:
                comparison_ok = compare_seqs (canonical_translation, translated_seq, verbose = verbose)
                if (comparison_ok):
                    #print "translation ok"
                    continue
            ###################### if ok, we are done here ######################

            # if translation does not match the reported canonical sequence,
            # find the alternative seq_region
            [orig_name, orig_file,  orig_start, orig_end] = [seq_name, file_names, 
                                                             seq_region_start, seq_region_end]
            ret = get_alt_seq_info (cursor, gene_id, species)
            if ret:
                [seq_name, file_names, seq_region_start, seq_region_end] = ret
                print "alt seq info ", seq_name, file_names, seq_region_start, seq_region_end
            else:
                ct +=1 
                continue
            # extract raw gene  region
            gene_seq = get_gene_seq( acg, species, seq_name, file_names, seq_region_strand,  
                                     seq_region_start, seq_region_end)

            # reconstruct the translation from the raw gene_seq and exon boundaries
            translated_seq = transl_reconstruct (cursor, gene_id, gene_seq, 
                                                 canonical_coding_exons, verbose = verbose)
            if (translated_seq):
                # compare the two sequences and cry foul if they are not the same:
                comparison_ok = compare_seqs (canonical_translation, translated_seq, verbose = verbose)
                if (comparison_ok):
                    continue

            if (not translated_seq or not comparison_ok):
                print "error translating", gene_id, gene2stable(cursor, gene_id)
                print "original name: ",  orig_name, orig_file,  orig_start, orig_end
                print "attempted fix: ",  seq_name, file_names,  seq_region_start, seq_region_end
                ct += 1
                print ct, "out of",  tot
                print "==========================================="
                print canonical_translation
                print "========"
                print translated_seq

 
        print "\t translation fail: ", ct, "out of ", tot
       
    cursor.close()
    db    .close()

#########################################
def main():

    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
  

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    check_canonical_sequence (local_db, all_species, ensembl_db_name)





#########################################
if __name__ == '__main__':
    main()






#########################################
'''
        #for gene_id in gene_ids[6000:8000]:
        #for gene_id in [314408,  314604,  314656,  314728,  314736,  314756,  314794,  314805,  314845,  314954,  314978,  314990,  315225,  315324,  315616,  315722,  315802,  315982,  316001,  316194,  319848,  320075,  320236,  320285,  320404,  320891,  368524,  368526,  368549,  368639,  368646,  368651,  368669,  368684,  368687,  368698,  368707,  368743,  368762,  368766,  368767,   368985,  369163,  369184,  369185,  369189,  369191,  369194,  369197,  369266,  369306,  369333,  369359,  369385,  369413,  369474,  369524]:

            print 
            print "==========================================="
            print "can transl id: ", gene2stable_canon_transl (cursor, gene_id)
            print canonical_translation
            print "========"
            print translated_seq

'''
