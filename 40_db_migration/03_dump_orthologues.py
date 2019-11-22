#!/usr/bin/python -u

import MySQLdb
import sys, commands
from   el_utils.mysql         import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl       import  *
from   el_utils.exon          import  Exon
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.threads       import  parallelize
from   el_utils.utils         import  erropen
from el_utils.ncbi    import  taxid2trivial

# jupiter

#+----------------------+-------------+------+-----+---------+----------------+
#| Field                | Type        | Null | Key | Default | Extra          |
#+----------------------+-------------+------+-----+---------+----------------+
#| orth_pair_id         | int(10)     | NO   | PRI | NULL    | auto_increment |
#| gene_id              | int(10)     | YES  | MUL | NULL    |                |
#| cognate_gene_id      | int(10)     | YES  |     | NULL    |                |
#| cognate_genome_db_id | int(11)     | YES  |     | NULL    |                |
#| source               | varchar(20) | YES  |     | NULL    |                |
#+----------------------+-------------+------+-----+---------+----------------+
#5 rows in set (0.00 sec)

#mysql> select * from orthologue limit 3;
#+--------------+---------+-----------------+----------------------+---------+
#| orth_pair_id | gene_id | cognate_gene_id | cognate_genome_db_id | source  |
#+--------------+---------+-----------------+----------------------+---------+
#|            1 |  370480 |           91863 |                  125 | ensembl |
#|            2 |  370480 |           29281 |                  123 | ensembl |
#|            3 |  370480 |           18633 |                  115 | ensembl |
#+--------------+---------+-----------------+----------------------+---------+

# reindeer

#+-----------------+-------------+------+-----+---------+----------------+
#| Field           | Type        | Null | Key | Default | Extra          |
#+-----------------+-------------+------+-----+---------+----------------+
#| id              | int(10)     | NO   | PRI | NULL    | auto_increment | 
#| ensembl_gene_id | varchar(30) | YES  | MUL | NULL    |                | 
#| cognate_gene_id | varchar(30) | YES  |     | NULL    |                | 
#| species         | varchar(50) | YES  |     | NULL    |                | 
#| common_name     | varchar(50) | YES  |     | NULL    |                | 
#+-----------------+-------------+------+-----+---------+----------------+
#+----+-----------------+--------------------+------------------------+-------------+
#| id | ensembl_gene_id | cognate_gene_id    | species                | common_name |
#+----+-----------------+--------------------+------------------------+-------------+
#|  1 | ENSG00000215906 | ENSCPOG00000022624 | cavia_porcellus        | guinea_pig  | 
#|  2 | ENSG00000215906 | ENSFCAG00000018472 | felis_catus            | cat         | 
#|  3 | ENSG00000215906 | ENSGACG00000004509 | gasterosteus_aculeatus | stickleback | 
#+----+-----------------+--------------------+------------------------+-------------+

#########################################
def translate_to_trivial(cursor, all_species):
    trivial_name = {}
    for species in all_species:
        taxid                 = species2taxid (cursor, species)
        trivial_name[species] = taxid2trivial(cursor, taxid)

    return trivial_name

#########################################
def orthos_tabstring(orthos_info):
    

    ret = "\t".join( map((lambda token:  type(token) is str and token or str(token)), orthos_info) )

    return ret

  
#########################################
def dump_orthos (species_list, db_info):

    
    [local_db, ensembl_db_name] = db_info
    db     = connect_to_mysql()
    cfg    = ConfigurationReader()
    cursor = db.cursor()

     # find db ids adn common names for each species db
    [all_species, ensembl_db_name] = get_species (cursor)

    # in the afa headers use 'trivial' names for the species: cow, dog, pig, ...
    trivial_name   = translate_to_trivial(cursor, all_species)

    out_path = cfg.get_path('afs_dumps')
    outfile  = "{0}/orthologue_dump.txt".format(out_path)
    print outfile
    of       = erropen (outfile,"w")

    species  = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])


    qry = "select * from orthologue"
    rows = search_db (cursor, qry)
    for row in rows:
        [pair_id, human_gene_id, cognate_gene_id, genome_db_id, source] =  row
        species = genome_db_id2species (cursor, genome_db_id)
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        human_stable_id = gene2stable(cursor, human_gene_id)
        switch_to_db (cursor,  ensembl_db_name[species])
        cognate_stable_id = gene2stable(cursor, cognate_gene_id)
        print  >>of,  orthos_tabstring ([human_stable_id, cognate_stable_id, species, trivial_name[species]])


    of.close()
    
    cursor.close()
    db    .close()

#########################################
def main():

    no_threads = 1

    db = connect_to_mysql()
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    parallelize (no_threads, dump_orthos, all_species, [local_db, ensembl_db_name])



#########################################
if __name__ == '__main__':
    main()
