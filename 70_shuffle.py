#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re, sys
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.utils import  *
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
# BioPython
from Bio.Seq      import Seq


#########################################
def main():

    local_db = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
        acg = AlignmentCommandGenerator()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    print "running ", sys.argv[0]

    [species, gene_id] = sys.argv[1:3]
    gene_id = int(gene_id)
    print species, gene_id

    switch_to_db (cursor, ensembl_db_name[species])

    gene_coords =  get_gene_coordinates (cursor, gene_id)  
    [gene_seq_region_id, gene_start, gene_end, gene_strand] = gene_coords
    qry  = "select name, file_name "
    qry += " from seq_region where seq_region_id = %d " % int(gene_seq_region_id)
    rows = search_db(cursor, qry)
    [name, file_names] = rows[0]
    filename = get_best_filename(file_names)
   
    exon_list = gene2exon_list(cursor, gene_id)

    for e in exon_list:
        print
        print e.exon_id, e.is_known, e.start_in_gene, e.phase, e.strand

        if ( gene_strand >  0 ):
            region_start = gene_start + e.start_in_gene
            region_end   = gene_start + e.end_in_gene

        else:
            region_end   = gene_end - e.start_in_gene
            region_start = gene_end - e.end_in_gene
        print "\t", region_start, region_end
        fasta    = get_fasta (acg, species, name, filename, gene_strand, region_start, region_end)
        qry_seq  = "".join(fasta.splitlines()[1:])
        pepseq   = Seq(qry_seq).translate().tostring()
        print pepseq

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

'''
    #for gene_id in [412667]: #  wls
    #for gene_id in [378768]: #  p53
    #for gene_id in [378766]: #  dynein
'''

