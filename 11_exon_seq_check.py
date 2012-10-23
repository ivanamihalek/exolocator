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
    
    # now for seq_name and seq_file
    if ( species == 'danio_rerio'):
        seq_region_id +=949 # beats me

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
    if ( pepseq[-1] == '*'):
        pepseq = pepseq[:-1]
    return pepseq

#########################################
def check_canonical_sequence(species_list, ensembl_db_name):

    db     = connect_to_mysql()
    cursor = db.cursor()
    acg     = AlignmentCommandGenerator()

    for species in species_list:
        if (not species == 'ailuropoda_melanoleuca'):
            continue
        print
        print "############################"
        print  species


        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        
        gene_ids.reverse()
        ct  = 0
        tot = 0
        for gene_id in gene_ids:
 
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

            canonical_translation = ""
            for line in fasta.split("\n"):
                if ('>' in line):
                    continue
                line.rstrip()
                canonical_translation += line

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
                print fasta_names
                exit(1)

            is_mitochondrial = ('.MT.' in fasta_db_file)

            # (2) extract seq from fasta db
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
                exit (1)
            gene_seq = ""
            for line in ret.split("\n"):
                if ('>' in line):
                    continue
                line.rstrip()
                gene_seq += line
                
            # find all exons associated with the gene id
            exons = gene2exon_list (cursor, gene_id)
            if (not exons):
                #ct +=1
                print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
                
            # sorting exons in place by their start in gene:
            exons.sort(key=lambda exon: exon.start_in_gene)

            spliced_seq = "" 
            carry = ""
            for exon in exons:
                if (not exon.is_canonical or not exon.is_coding): 
                     continue
                # find exon sequence within the gene
                exon_seq     = carry + gene_seq[ exon.start_in_gene:exon.end_in_gene+1]
                remainder    = len(exon_seq)%3
                if ( remainder == 0 ):
                    carry = ""
                elif (remainder == 1 ):
                    carry    = exon_seq[-1:]
                    exon_seq = exon_seq[:-1]
                else:
                    carry    = exon_seq[-2:]
                    exon_seq = exon_seq[:-2]

                dnaseq = Seq (exon_seq, generic_dna)
                pepseq = dnaseq.translate()
                spliced_seq += exon_seq
                
            # turn to the corresponding BioPython object
            dnaseq = Seq (spliced_seq, generic_dna)
            pepseq = dnaseq.translate()
            if ( is_mitochondrial ):
                pepseq = dnaseq.translate(table="Vertebrate Mitochondrial")
            else:
                pepseq = dnaseq.translate()
                pepseq = strip_stop(pepseq)  # strip the last stop codon only
                
                # if there are stil stop codons we'll give another shot to
                # to the possibility that tit si mitochondrial,
                # after that we cry foul
                if ( '*' in pepseq):
                    pepseq = dnaseq.translate(table="Vertebrate Mitochondrial")
                    pepseq = strip_stop(pepseq)
                         
                if ( '*' in pepseq): # some more desperate measures
                    dnaseq = Seq (spliced_seq[1:], generic_dna)
                    pepseq = dnaseq.translate()
                    pepseq = strip_stop(pepseq)

                if ( '*' in pepseq): # some more desperate measures
                    dnaseq = Seq (spliced_seq[2:], generic_dna)
                    pepseq = dnaseq.translate()
                    pepseq = strip_stop(pepseq)

                if ( '*' in pepseq):
                    ct += 1
                    print "Error translating ", species, gene_id, canonical_transl_id
                    print "canon: ", canonical_translation
                    print "exons: ", pepseq
                    print "dna:   ", spliced_seq

            difference = len(pepseq) - len(canonical_translation)
            if ( abs(difference) > 3):
                ct += 1
                print canonical_transl_id
                print "canon", canonical_translation
                print "exons", pepseq
                print
            else:
                for i in range(len(pepseq)):
                    if ( i >= len(canonical_translation)):
                        break
                    if (not pepseq[i] ==  canonical_translation[i]):
                        ct += 1
                        print canonical_transl_id
                        print "canon", canonical_translation
                        print "exons", pepseq
                        print canonical_transl_id, i, pepseq[i], canonical_translation[i]
                        print
                        break

            if (not  tot%500):
                print ct, tot

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
