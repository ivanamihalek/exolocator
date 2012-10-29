#!/usr/bin/python

import MySQLdb
import commands

from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl import  get_species, get_gene_ids, get_compara_name
from   el_utils.ensembl import  gene2stable, gene2stable_canon_transl
from   el_utils.tree    import  Tree, Node
from   el_utils.threads import  parallelize


#########################################
def track_ancestry(cursor, leaf):
    
    # find member id for the species:
    qry = "select taxon_id from genome_db where name='%s'" % leaf.name
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose=True)
        exit(1)

    taxon_id = rows[0][0]
    print leaf.name, taxon_id

    

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    ensembl_compara_name           = get_compara_name(cursor)
    # cursor now points to the ensembl_compare database:

    tree   = Tree()
    for species in all_species:
        leaf = Node(species)
        tree.leafs.append(leaf)

    tree.build(cursor)

    print
    print tree.nhx_string()
    print

    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
