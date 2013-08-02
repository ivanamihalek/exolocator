#!/usr/bin/python

import MySQLdb, subprocess, re, commands
import copy, pdb

from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.utils   import  *
from el_utils.map     import  Map, get_maps
from el_utils.tree    import  species_sort

from el_utils.special_gene_sets  import get_theme_ids, get_complement_ids
from el_utils.almt_cmd_generator import AlignmentCommandGenerator
from el_utils.config_reader      import ConfigurationReader
from el_utils.threads import  parallelize
from subprocess import Popen, PIPE, STDOUT

#########################################
verbose = True

#########################################
def annotate(gene_list, db_info):
    # 
    [local_db, ensembl_db_name, species] = db_info
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
        acg = AlignmentCommandGenerator()
    else:
        db  = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    print "thread %s annotating %s "% (get_thread_name(), species)

    # find annotation
    # find splices
    # output

    cursor.close()
    db.close()

#########################################
#########################################
def main():
    
    species    = 'oryctolagus_cuniculus'
    no_threads = 1
 

    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <species> <number of threads>" % sys.argv[0]
        exit(1)
    elif len(sys.argv)==3:

        species = sys.argv[1].lower()

        no_threads = int(sys.argv[2])
        

    local_db   = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)


    print '======================================='
    print "running %s for %s " % (sys.argv[0], species)

    switch_to_db (cursor,  ensembl_db_name['species'])
    gene_list = get_gene_ids (cursor, biotype='protein_coding')

    cursor.close()
    db.close()

    parallelize (no_threads, annotate, gene_list, [local_db, ensembl_db_name, species])
    
    return True

#########################################
if __name__ == '__main__':
	main()

