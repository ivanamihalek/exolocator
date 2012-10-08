#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  get_species, get_gene_ids
from   el_utils.objects import  Exon
from   el_utils.threads import  parallelize



#########################################
def get_gene_region (cursor, gene_id, is_known=None):

    qry     = "SELECT seq_region_id, seq_region_start, seq_region_end, "
    qry    += " seq_region_strand, canonical_transcript_id "
    qry    += " FROM gene WHERE gene_id=%d"  %  gene_id
    if (not is_known is None):
        qry    += " and is_known=%d"  % is_known
    rows    = search_db (cursor, qry, verbose=False)

    if (not rows):
        rows    = search_db (cursor, qry, verbose=True)
        return []
    elif ( 'Error' in rows[0]):
        print  rows[0]
        return []

    return rows[0]

#########################################
def get_known_exons (cursor, gene_id, species):

    exons = []

    # get the region on the gene
    ret = get_gene_region (cursor, gene_id)
    if  ret:
        [gene_seq_id, gene_region_start, gene_region_end, 
         gene_region_strand, canonical_transcript_id] = ret
    else:
        print "region not retrived for ", species, gene_id
        exit(1)

    return exons

#########################################
def get_predicted_exons (cursor, gene_id, species):

    exons = []

    # get the region on the gene
    ret = get_gene_region (cursor, gene_id)
    if  ret:
        [gene_seq_id, gene_region_start, gene_region_end, 
         gene_region_strand, canonical_transcript_id] = ret
    else:
        print "region not retrived for ", species, gene_id
        exit(1)


    qry    = "SELECT  * FROM prediction_exon WHERE seq_region_id = %d "  %  gene_seq_id
    qry   += " AND  seq_region_start >= %d AND seq_region_start <= %d " %  \
        (gene_region_start, gene_region_end)
    qry   += " AND  seq_region_end   >= %d AND seq_region_end   <= %d " %  \
        (gene_region_start, gene_region_end)
    rows   = search_db (cursor, qry)

    if (not rows):
        return []
    for row in rows:
        exon  = Exon()
        exon.load_from_ensembl_prediction (gene_region_start, row)
        exons.append(exon)

    return exons

#########################################
def find_exons (cursor, gene_id, species):

    coding_region_start = -1
    coding_region_end   = -1
    exons               = []
    
    # get all exons from the 'exon' table
    exons = get_known_exons (cursor, gene_id, species)
    # get all exons from the 'predicted_exon' table
    if (not species == 'homo_sapiens'):
        exons += get_predicted_exons (cursor, gene_id, species)

    return [coding_region_start, coding_region_end, exons]

#########################################
def gene2exon(species_list, ensembl_db_name):

    db     = connect_to_mysql()
    cursor = db.cursor()


    for species in species_list:
        print  species
        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        
        for gene_id in gene_ids:
            # find all exons associated with the gene id
            [coding_region_start, coding_region_end, exons] = \
                find_exons (cursor, gene_id, species)
            if (not exons):
                print  gene_id, " no exons found" 
                exit(1)

            print gene_id, " number of  exons: ", len(exons)
            for exon in exons:
                print
                print exon
            exit(1)

    cursor.close()
    db.close()

    return True



#########################################
def main():

    no_threads = 1

    db     = connect_to_mysql()
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    parallelize (no_threads, gene2exon, all_species, ensembl_db_name)



#########################################
if __name__ == '__main__':
    main()
