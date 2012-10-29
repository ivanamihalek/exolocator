#!/usr/bin/python

import MySQLdb
import commands
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  get_species, get_gene_ids, gene2exon_list
from   el_utils.ensembl import  gene2stable, gene2stable_canon_transl
from   el_utils.exon    import  Exon
from   el_utils.threads import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator




            


#########################################
def check_exons(species_list, ensembl_db_name):

    db     = connect_to_mysql()
    cursor = db.cursor()

    acg = AlignmentCommandGenerator()
    species_list = ['danio_rerio']
    #species_list = ['callithrix_jacchus']
    #species_list = ['ailuropoda_melanoleuca']
    for species in species_list:
        print
        print "############################"
        print  species
        #if (species in ['ailuropoda_melanoleuca', 'anolis_carolinensis', 
        #                'bos_taurus','danio_rerio']):
        #    continue
        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        
        ct = 0
        tot = 0
        #for gene_id in [1]:
        #for gene_id in [21459]:
        for gene_id in gene_ids:

            tot += 1
            # find all exons associated with the gene id
            exons = gene2exon_list (cursor, gene_id)
            if (not exons):
                #ct +=1
                print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
                
   
            length = 0
            for exon in exons:
                if (not exon.is_canonical): 
                    continue
                
                if (not exon.is_coding): 
                    continue
                

                #print
                #print exon
                #print exon.start_in_gene, exon.end_in_gene

                if (exon.translation_ends is None):
                    length += exon.end_in_gene - exon.start_in_gene + 1
                else:
                    length += exon.translation_ends + 1

                if (not exon.translation_starts is None):
                    length -= exon.translation_starts
            
            if (not length):
                print gene2stable (cursor, gene_id = gene_id), " no exons marked as canonical"
                continue

                
            # what is the length of the canonical transcript according to Ensembl
            canonical_transl_id = gene2stable_canon_transl(cursor, gene_id)
            if ( not canonical_transl_id):
                print "no canonical transl id found for ", gene_id
                continue


            # get canonical transcript from ensembl fasta database
            cmd = acg.generate_fastacmd_protein_command (canonical_transl_id, species, 
                                                         "all", None)
            fasta = commands.getoutput(cmd)
            
            if (not fasta):
                print gene2stable (cursor, gene_id = gene_id), "fasta not found for ", canonical_transl_id
                continue
            
            translation_length = 0
            for line in fasta.split('\n'):
                if '>' in line:
                    continue
                translation_length += len(line)


            if ( abs(length/3-translation_length) > 3):
                ct +=1
                print gene2stable (cursor, gene_id = gene_id),
                print "exon length ", length/3, " transl len ", translation_length,
                print "(", ct, " out of ", tot, ")"
            

            if (tot==1000):
                break
            #print fasta
   
 
    cursor.close()
    db.close()

    return True



#########################################
def main():

    no_threads = 1

    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    parallelize (no_threads, check_exons, all_species, ensembl_db_name)



#########################################
if __name__ == '__main__':
    main()
