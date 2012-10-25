#!/usr/bin/python

import MySQLdb
import commands
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  get_species, get_gene_ids, gene2exon_list
from   el_utils.ensembl import  gene2stable, gene2stable_canon_transl
from   el_utils.objects import  Exon
from   el_utils.threads import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator

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
def strip_stop(pepseq):
    if (not pepseq or len(pepseq)==0):
        return pepseq
    if ( pepseq[-1] == '*'):
        pepseq = pepseq[:-1]
    return pepseq

#########################################
def  canonical_transl_info (cursor, gene_id):
    
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
def check_canonical_sequence(species_list, ensembl_db_name):

    db     = connect_to_mysql()
    cursor = db.cursor()
    acg     = AlignmentCommandGenerator()

    for species in species_list:
        if (not species == 'homo_sapiens'):
            continue
        #if (not species == 'danio_rerio'):
        #    continue
        print
        print "############################"
        print  species


        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        
        
 
        #gene_ids.reverse()
        ct  = 0
        tot = 0
        #for gene_id in gene_ids:
        for gene_id in [314403]:
 
            [can_transl_start_exon, can_transl_start_position,
             can_transl_end_exon, can_transl_end_position] = canonical_transl_info (cursor, gene_id)

            canonical_transl_id = gene2stable_canon_transl(cursor, gene_id)
            # find canonical translation
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

            # find gene sequence:
            # (1) find seq_name and region
            tot +=1 
            ret = get_seq_info (cursor, gene_id, species)
            if ret:
                [seq_name, file_names,  seq_region_start, 
                 seq_region_end, seq_region_strand] = ret
            else:
                ct +=1 
                continue

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


            print ">>> ", fasta_db_file
            print ">>> ", file_names

            is_mitochondrial = ('.MT.' in fasta_db_file)

            # extract gene sequence  from fasta db
            fastacmd = acg.generate_fastacmd_gene_command(species, seq_name, fasta_db_file,
                                                          seq_region_strand,  seq_region_start,    
                                                          seq_region_end)

            print fastacmd
            ret = commands.getoutput(fastacmd)
            if not ret:
                print "no refturn for fastacmd for", species, gene_id
                print "fastacmd: ", fastacmd
                exit (1)
            if ('ERROR' in ret):
                print "Error running fastacmd: ", fastacmd
                print ret
                exit (1)
            gene_seq = ""
            for line in ret.split("\n"):
                if ('>' in line):
                    continue
                line.rstrip()
                gene_seq += line

                
            print seq_region_end-seq_region_start+1

            #print gene_seq

                
            ################################################################
            # find all exons associated with the gene id
            exons = gene2exon_list (cursor, gene_id)
            if (not exons):
                print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
                
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
                
            translated_seq = "" 
            carry = ""
            ok_so_far = True
            # sanity checking
            for exon in canonical_coding_exons:
                
                #print
                #print "exon", exon.exon_id
                # find exon sequence within the gene
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
                # strip the last stop codon only
                if ( exon is canonical_coding_exons[-1]):
                    pepseq = strip_stop(pepseq)  

                pepseq0 =  pepseq
                print "phase 0", pepseq0
                print "\t ", exon_seq_for_transl_purposes

                # if there are stil stop codons we'll give another shot to
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
                    print "phase 1", pepseq1
                else:
                    ok_so_far = True

                if (not ok_so_far and '*' in pepseq):
                    dnaseq = Seq (exon_seq_for_transl_purposes[2:], generic_dna)
                    pepseq = dnaseq.translate()
                    pepseq = strip_stop(pepseq)
                    pepseq2 =  pepseq
                    print "phase 2", pepseq2
                else:
                    ok_so_far = True

                if (not ok_so_far and  '*' in pepseq):
                    ct += 1
                    print 
                    print ct, tot
                    print "Error translating ", species, gene_id, canonical_transl_id,
                    print " (stop codon) "
                else:
                    ok_so_far = True
           
 
                if ( not ok_so_far):
                    break

                translated_seq += pepseq

            if ( not ok_so_far):
                continue

            while (len(translated_seq) and translated_seq[0] == 'X'):
                translated_seq = translated_seq[1:]
                
            difference = len(translated_seq) - len(canonical_translation)
            if ( abs(difference) > 3):
                ct += 1
                print
                print ct, tot
                print gene_id, canonical_transl_id
                print ">canon"
                print canonical_translation
                print ">exons"
                print translated_seq
                print
            else:
                diff = 0
                start = -1
                for i in range(len(translated_seq)):
                    if ( i >= len(canonical_translation)):
                        break
                    if (not translated_seq[i] ==  canonical_translation[i]):
                        diff += 1
                        if start < 0:
                            start = i

                if (diff > 2):
                    ct += 1 
                    print
                    print ct, tot
                    print gene_id, canonical_transl_id
                    print "canon"
                    print canonical_translation
                    print "exons"
                    print translated_seq
                    print canonical_transl_id,  translated_seq[start], canonical_translation[start]
                    print "nuber of  diff sites: ", diff, " starting from ", start
                    print
                   

            if (not  tot%2000):
                #print ct, tot
                break

        print species, ct, tot
       
    cursor.close()
    db    .close()

#########################################
def main():

    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    check_canonical_sequence (all_species, ensembl_db_name)





#########################################
if __name__ == '__main__':
    main()


'''
            tot += 1
                pass
                qry  = "select seq_region_start, seq_region_end from exon"
                qry += " where exon_id = %d " % exon.exon_id
                rows = search_db (cursor, qry)
                [exon_start, exon_end] = rows[0]
                fastacmd = acg.generate_fastacmd_gene_command(species, seq_name, fasta_db_file,
                                                          seq_region_strand, exon_start,    
                                                          exon_end)
                
                ret = commands.getoutput(fastacmd)
                if not ret:
                    print "no refturn for fastacmd for", species, gene_id
                    print "fastacmd: ", fastacmd
                    exit (1)

                exon_seq = ""
                for line in ret.split("\n"):
                    if ('>' in line):
                        continue
                    line.rstrip()
                    exon_seq += line

                exon_dnaseq = Seq (exon_seq, generic_dna)
                pepseq      = exon_dnaseq.translate()
                print pepseq
                
'''
