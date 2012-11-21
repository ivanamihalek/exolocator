#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re

from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update

from   el_utils.ensembl import  *
from   el_utils.map     import  Map, get_maps
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader

from bitstring import Bits



#########################################
def main():

    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    gene_id = 355489


    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
    human_exons = gene2exon_list(cursor, gene_id)
    print gene2stable(cursor, gene_id), get_description (cursor, gene_id)

    tetra = 'tetraodon_nigroviridis'
    tetra_genome_db_id = species2genome_db_id(cursor, tetra)
    print tetra, tetra_genome_db_id
    print

    for human_exon  in  human_exons:
        if (not human_exon.is_canonical):
            continue

        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        qry  = "select cognate_exon_id, cognate_exon_known, cigar_line, msa_bitmap  from exon_map where "
        qry += " cognate_genome_db_id = %d " % tetra_genome_db_id 
        qry += " and exon_id  = %d " % human_exon.exon_id
        qry += " and exon_known = %d " % human_exon.is_known

        rows = search_db (cursor, qry)
        print "----------------------------------"
        print human_exon.exon_id, human_exon.is_known
        if rows:
            for row in rows: 
                if ( len(row) == 4):
                    [exon_id_2, exon_known_2, cigar_line, msa_bitmap] = row
                    gene_2 = exon_id2gene_id(cursor, ensembl_db_name[tetra], exon_id_2, exon_known_2)
                    print  exon_id_2, exon_known_2, gene_2, cigar_line, "|||||",  msa_bitmap
                else:
                    print row
#########################################
if __name__ == '__main__':
    main()
