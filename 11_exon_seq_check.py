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
    for species in species_list:
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
        
        
        for gene_id in gene_ids[:10]:
        #for gene_id in [412667]:
        #for tot in range(1000):
            #gene_id = choice(gene_ids)
            tot +=1 
 
            # get _all_ exons
            exons = gene2exon_list(cursor, gene_id, ensembl_db_name[species])
            if (not exons):
                print 'no exons for ', gene_id
                ct += 1
                continue

            # extract raw gene  region - bonus return from checking whether the 
            # sequence is correct: translation of canonical exons
            [gene_seq, canonical_exon_pepseq, file_names] = get_gene_seq(acg, cursor, gene_id, species)
            if (not gene_seq or not canonical_exon_pepseq):
                ct += 1
                print 'no sequence found for ', gene_id, "   ",   ct, "out of ", tot
                continue


 
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
