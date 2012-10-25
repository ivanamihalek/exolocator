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
def store_exon_seqs(species_list, ensembl_db_name):

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
                print fasta_names
                exit(1)

            is_mitochondrial = ('.MT.' in fasta_db_file)

            # extract sequencce from fasta db
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
            for exon in exons:  # find exon sequence within the gene
 
                exon_seq                     =         gene_seq[ exon.start_in_gene:exon.end_in_gene+1]
                if (exon.is_canonical and exon.is_coding): 
               
                    exon_seq_for_transl_purposes = carry + gene_seq[ exon.start_in_gene:exon.end_in_gene+1]
                    remainder    = len(exon_seq)%3
                    if ( remainder == 0 ):
                        carry = ""
                    elif (remainder == 1 ):
                        carry    = exon_seq_for_transl_purposes[-1:]
                        exon_seq = exon_seq_for_transl_purposes[:-1]
                    else:
                        carry    = exon_seq_for_transl_purposes[-2:]
                        exon_seq = exon_seq_for_transl_purposes[:-2]
                
                else:
                    if (exon.phase == 0):
                        exon_seq_for_transl_purposes = exon_seq
                    elif (exon.phase == 1):
                        exon_seq_for_transl_purposes = exon_seq[1:]
                    else:
                        exon_seq_for_transl_purposes = exon_seq[2:]

                # BioPython object
                dnaseq = Seq (exon_seq_for_transl_purposes, generic_dna)

                if ( is_mitochondrial ):
                    pepseq = dnaseq.translate("Vertebrate Mitochondrial")
                else:
                    pepseq = dnaseq.translate()
                
                

                # store to exon_seq table
                is_smith_waterman = 0;
                store_seq (cursor, exon.exon_id, exon.is_known, is_smith_waterman,
                            exon_seq, left_seq, right_seq, pepseq.string())
                # store the resiprocal info: exon_seq_id to gene2exon table
                sore_exon_seq_id (cursor, exon, exon_seq_id)

            if (not  tot%500):
                print ct, tot

        print species, ct, tot
       
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
