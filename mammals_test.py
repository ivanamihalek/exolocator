#!/usr/bin/python

import MySQLdb
from random import choice
from el_utils.mysql   import  *
from el_utils.ensembl import  *
from el_utils.ncbi    import  *
from el_utils.config_reader      import ConfigurationReader


#########################################
def main():
    

    local_db   = False

    if local_db:
        db  = connect_to_mysql()
        #cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        #cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    mammals = find_mammals(cursor, all_species)
    for mammal in mammals:
        print mammal
       

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

'''
       'telomere_maintenance',
                     'nonhom_end_joining', 'egfr_signaling', 'cell_cycle_checkpoints',
                     'genecards_top500', 'wnt_pathway',  'enzymes'

       prev_end = 0
        for human_exon in human_exons:

            first_exon = (human_exons.index(human_exon) == 0)
            if  first_exon: origin = human_exon.start_in_gene
            print " &&&   %10d    %10d   %10d" % (human_exon.start_in_gene-origin, 
                                                     human_exon.end_in_gene-human_exon.start_in_gene+1,
                                                     human_exon.start_in_gene - prev_end)
            prev_end = human_exon.end_in_gene

'''
