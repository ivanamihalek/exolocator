#########################################

import os
from mysql import *
from utils import erropen


#########################################
def get_theme_ids(cursor, ensembl_db_name, config_reader, theme_name):

    path = config_reader.get_path('resources')
    fnm  = path + '/' + theme_name + '.txt'
    if not os.path.exists(fnm):
        print fnm, "not found"
        exit(1)

    if not os.path.getsize(fnm) > 0:
        print fnm, "empty"
        exit(1)

    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        
    inf = erropen(fnm, "r")
    gene_ids = []
    for line in inf:
        line = line.rstrip()
        if "\t" in line:
            fields = line.split("\t")
            stable_id = fields[0]
        elif "\s" in line:
            fields = line.split("\t")
            stable_id = fields[0]
        else:
            stable_id = line
        qry  = "select gene_id, description from gene where stable_id='%s'" % stable_id
        rows = search_db (cursor, qry)
        if not rows: continue
        if 'ERROR' in rows[0]:
            print rows[0]
            exit(1)
        gene_ids.append(int(rows[0][0]))
    inf.close()

    return gene_ids

#########################################
def  human_genes_w_sw_sharp_annotation (cursor, ensembl_db_name):

    genes_with_patch = 0

    human_exons = []
    for db in ensembl_db_name.values():
        switch_to_db (cursor, db)
        qry  = "select distinct(maps_to_human_exon_id) from sw_exon"
        rows = search_db (cursor, qry)
        if not rows: continue
        for row in rows:
            if not row[0] in human_exons:
                human_exons.append(row[0])
 
    human_genes = []
    switch_to_db (cursor,ensembl_db_name['homo_sapiens'])
    for exon_id in human_exons:
        qry = "select gene_id from gene2exon where exon_id = %d " % exon_id
        rows = search_db (cursor, qry)
        if not rows: continue
        for row in rows:
            if not row[0] in human_genes:
                human_genes.append(row[0])
        
    
    print "human exons: ", len(human_exons),  "from genes: ", len(human_genes)


    return human_genes

