#!/usr/bin/python

import MySQLdb
import commands, operator
from   el_utils.mysql   import  connect_to_mysql, search_db, store_or_update
from   el_utils.ensembl import  *

#########################################
def gene_name2gene_id(cursor, gene_name):
    
    qry  = "select ensembl_id from object_xref, external_synonym "
    qry += "where object_xref.ensembl_object_type = 'Gene'  "
    qry += "and object_xref.xref_id= external_synonym.xref_id "
    qry += "and external_synonym.synonym = '%s' " % gene_name
    qry += "group by synonym"
    rows = search_db (cursor, qry)
    if not rows:
        return ""
    else:
        return rows[0][0]
#########################################
def search_description (cursor, gene_name):
    qry  = "select gene_id, description from gene "
    qry += "where description like '%"+gene_name+"%'"
    rows = search_db (cursor, qry)
    if not rows:
        return ["", ""]
    else:
        return rows[0]


#########################################
#########################################
def main():

    local_db   = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    switch_to_db(cursor, ensembl_db_name['homo_sapiens'])

    tot = 0
    number_of_patched_exons = {}
    for species in all_species:

        switch_to_db(cursor, ensembl_db_name[species])
        qry  =  "select count(1) from sw_exon"
        rows = search_db(cursor, qry)

        if rows and rows[0][0]:
            number_of_patched_exons[species] = int (rows[0][0])
            
    sorted_species = sorted(number_of_patched_exons.iteritems(), key=operator.itemgetter(1))
    sorted_species.reverse()
    #####################
    seen    = {}
    patches = {}
    for species, no_exons in sorted_species:
        
        switch_to_db(cursor, ensembl_db_name[species])
        qry  = "select sw_exon.gene_id, sw_exon.maps_to_human_exon_id, exon_seq.exon_id, exon_seq.protein_seq from exon_seq, sw_exon"
        qry += " where exon_seq.exon_id = sw_exon.exon_id and exon_seq.is_sw = 1"
        rows = search_db (cursor, qry, verbose=True)
        ct = 0
        for row in rows:
            [gene_id, human_exon_id, exon_id, pepseq] = row
            if '*' in pepseq:    continue
            if 'XXX' in pepseq:    continue
            if len(pepseq) < 10: continue
            ct += 1
            human_gene_id = exon_id2gene_id (cursor, ensembl_db_name['homo_sapiens'], human_exon_id, 1)
        
            #print ct, exon_id, pepseq
            #print "\t", gene2stable(cursor, human_gene_id, ensembl_db_name['homo_sapiens']),
            #print get_description (cursor, human_gene_id, ensembl_db_name['homo_sapiens'])

            if not seen.has_key(human_gene_id):
                seen[human_gene_id]    = 1
                patches[human_gene_id] = []
            else:
                seen[human_gene_id] += 1
            patches[human_gene_id].append([species, exon_id, pepseq])

    sorted_human_genes = sorted(seen.iteritems(), key=operator.itemgetter(1))
    sorted_human_genes.reverse()
    #####################
    for human_gene_id, ct in  sorted_human_genes[1:5]:
        print
        print "most seen: ",  human_gene_id, ct
        print "\t", gene2stable(cursor, human_gene_id, ensembl_db_name['homo_sapiens']),
        print get_description (cursor, human_gene_id, ensembl_db_name['homo_sapiens'])
        
        for [species, exon_id, pepseq] in patches[human_gene_id]:
            print "\t", species, exon_id, pepseq

    cursor.close()
    db    .close()
    



#########################################
if __name__ == '__main__':
    main()
