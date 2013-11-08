#!/usr/bin/python

import MySQLdb
import commands
from random import choice
from   el_utils.mysql       import  *
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
def main():

    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
  

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    swtich_to_db (cursor, ensemb_db['homo_sapiens'])
    qry  = "select gene_id, stable_id, seq_region_id, seq_region_start, seq_region_end from gene where stable_id"
    qry += " = 'ENSG00000182809' or stable_id = 'ENSG00000270931'"
    rows  = search_db (cursor, qry)
    for row in rows:
        print row


    cursor.close()
    db    .close()

 



#########################################
if __name__ == '__main__':
    main()






#########################################
        #for gene_id in gene_ids[6000:8000]:
        #for gene_id in [314408,  314604,  314656,  314728,  314736,  314756,  314794,  314805,  314845,  314954,  314978,  314990,  315225,  315324,  315616,  315722,  315802,  315982,  316001,  316194,  319848,  320075,  320236,  320285,  320404,  320891,  368524,  368526,  368549,  368639,  368646,  368651,  368669,  368684,  368687,  368698,  368707,  368743,  368762,  368766,  368767,   368985,  369163,  369184,  369185,  369189,  369191,  369194,  369197,  369266,  369306,  369333,  369359,  369385,  369413,  369474,  369524]:
