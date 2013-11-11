#!/usr/bin/python

import MySQLdb
import commands
from random import choice
from   el_utils.mysql       import  *
from   el_utils.ensembl     import  *
from   el_utils.exon        import  Exon
from   el_utils.threads     import  parallelize
from   el_utils.almt_cmd_generator  import AlignmentCommandGenerator
from   el_utils.special_gene_sets   import  get_theme_ids
from   el_utils.config_reader       import ConfigurationReader

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
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

    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    qry  = "select gene_id, stable_id, seq_region_id, seq_region_start, seq_region_end, description "
    qry += "from gene where stable_id"
    qry += " = 'ENSG00000182809' or stable_id = 'ENSG00000270931'"
    rows  = search_db (cursor, qry)
    for row in rows:
        print row

    qry  = "select seq_region.name,  assembly_exception.exc_seq_region_start, assembly_exception.exc_seq_region_end "
    qry += "from seq_region, assembly_exception  where seq_region.seq_region_id = assembly_exception.exc_seq_region_id "
    qry += "and assembly_exception.seq_region_id = 1001061114 and not assembly_exception.exc_type = 'PAR'"
    rows  = search_db (cursor, qry)
    for row in rows:
        print row

    qry  = "select gene_id, alt_allele_group_id from alt_allele "
    qry += "where gene_id = %d " % 672092
    qry += "or  gene_id = %d "   % 695872
    rows  = search_db (cursor, qry)
    for row in rows:
        print row

    # how many human genes have multiple "alleles"
    group = {}
    qry = "select  alt_allele_group_id, group_concat(alt_allele_id)  from alt_allele group by alt_allele_group_id"
    rows  = search_db (cursor, qry)
    for row in rows:
        group[row[0]] = row[1]

    print "there are ", len(group), " alt_allele_groups"
    for alt_allele_group_id, group in group.iteritems():
        al_ids = map (lambda ai:  int(ai), group.split(","))
        print alt_allele_group_id, len(al_ids)
        for al_id in al_ids:
            qry = "select attrib from alt_allele_attrib "
            qry += "where alt_allele_id = %d " % al_id
            rows  = search_db (cursor, qry)
            if rows and rows[0]:
                attrib = rows[0][0]
            else:
                attrib = ""
            print "\t", al_id,  attrib

    cursor.close()
    db    .close()

 



#########################################
if __name__ == '__main__':
    main()


