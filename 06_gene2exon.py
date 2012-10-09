#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  get_species, get_gene_ids, gene2stable
from   el_utils.objects import  Exon
from   el_utils.threads import  parallelize



#########################################
def get_gene_region (cursor, gene_id, is_known=None):

    qry     = "SELECT seq_region_id, seq_region_start, seq_region_end, "
    qry    += " seq_region_strand "
    qry    += " FROM gene WHERE gene_id=%d"  %  gene_id
    if (not is_known is None):
        qry  += " and is_known=%d"  % is_known
    rows    = search_db (cursor, qry, verbose=False)

    if (not rows):
        rows    = search_db (cursor, qry, verbose=True)
        return []
    elif ( 'Error' in rows[0]):
        print  rows[0]
        return []

    return rows[0]

#########################################
def get_canonical_transcript_id (cursor, gene_id):

    qry     = "SELECT canonical_transcript_id"
    qry    += " FROM gene WHERE gene_id=%d"  %  gene_id
    rows    = search_db (cursor, qry, verbose=False)

    if (not rows):
        rows    = search_db (cursor, qry, verbose=True)
        return ""
    elif ( 'Error' in rows[0]):
        print  rows[0]
        return ""


    return rows[0][0]

#########################################
def get_exons (cursor, gene_id, species, table):

    exons = []

    # get the region on the gene
    ret = get_gene_region (cursor, gene_id)
    if  ret:
        [gene_seq_id, gene_region_start, gene_region_end, 
         gene_region_strand] = ret
    else:
        print "region not retrived for ", species, gene_id
        return []

    qry    = "SELECT  * FROM " + table + " WHERE seq_region_id = %d "  %  gene_seq_id
    qry   += " AND  seq_region_start >= %d AND seq_region_start <= %d " %  \
        (gene_region_start, gene_region_end)
    qry   += " AND  seq_region_end   >= %d AND seq_region_end   <= %d " %  \
        (gene_region_start, gene_region_end)
    rows   = search_db (cursor, qry)

    if (not rows):
        return []
    for row in rows:
        exon  = Exon()
        if ( table == 'exon'):
            exon.load_from_ensembl_exon (gene_region_start, row)
        else:
            exon.load_from_ensembl_prediction (gene_region_start, row)
        exon.gene_id = gene_id
        exons.append(exon)
 
    return exons


#########################################
def get_canonical_exon_ids (cursor, canonical_transcript_id):

    canonical_exon_ids = []
    qry = "select exon_id from exon_transcript "
    qry += " where transcript_id = %d " % canonical_transcript_id
    rows   = search_db (cursor, qry)
    if (not rows):
        return []
    for row in rows:
        canonical_exon_ids.append(row[0])

    return canonical_exon_ids


#########################################
def  mark_canonical (cursor, gene_id, exons):

    canonical_transcript_id = get_canonical_transcript_id(cursor, gene_id)
    if not canonical_transcript_id:
        print "canonical_transcript_id  not retrived for ",  gene_id
        return False

    canonical_exon_ids = get_canonical_exon_ids (cursor, canonical_transcript_id)
    for exon in exons:
        if exon.is_known and (exon.exon_id in canonical_exon_ids):
            exon.is_canonical = 1
        else:
            exon.is_canonical = 0
            
#########################################
def fill_in_annotation_info (cursor, gene_id, exons):

    for exon in exons:
        if ( exon.is_known):
            qry  = "select transcript.analysis_id "
            qry += " from transcript, exon_transcript "
            qry += " where exon_transcript.transcript_id = transcript.transcript_id "
            qry += " and  exon_transcript.exon_id = %d " % exon.exon_id
        else:
            qry  = " select prediction_transcript.analysis_id "
            qry += " from prediction_transcript, prediction_exon "
            qry += " where prediction_exon.prediction_transcript_id "
            qry += " = prediction_transcript.prediction_transcript_id and  "
            qry += " prediction_exon.prediction_exon_id= %d " %  exon.exon_id

        rows   = search_db (cursor, qry)

        if (not rows):
            exon.analysis_id = 0
        else:
            exon.analysis_id = rows[0][0]

###################################
def get_logic_name(analysis_id, cursor):
        qry = "SELECT logic_name FROM analysis WHERE analysis_id = %d" % analysis_id
        rows    = search_db (cursor, qry)
        if (not rows):
            logic_name = ''
        else:
            logic_name = rows[0][0]
        return logic_name


#########################################
def sort_out_covering_exons (cursor, exons):

    # havana is manually curated an gets priority
    is_ensembl = {}
    is_havana  = {}
    for exon in exons:
        logic_name = get_logic_name(exon.analysis_id, cursor)
        is_ensembl[exon] = ('ensembl' in logic_name)
        is_havana [exon] = ('havana'  in logic_name)


    # find a representative for each overlapping cluster
    # if we can find havana, it will be havana
    # otherwise we move to less reliable annotation
    overlapping_cluster = {}
    for i in range(len(exons)):

        exon_1 = exons[i]

        already_seen = 0
        for exon_2 in overlapping_cluster.keys():
            if exon_1 in overlapping_cluster[exon_2]:
                already_seen = 1
                break
        if (already_seen):
            continue

        exon_start_1 = exon_1.start_in_gene
        exon_end_1   = exon_1.end_in_gene

        the_master_exon_is = -1

        for j in range(i+1, len(exons)):
     
            exon_2       = exons[j]
            exon_start_2 = exon_2.start_in_gene
            exon_end_2   = exon_2.end_in_gene

            the_master_exon_is = -1
            
            if (exon_start_1 <= exon_start_2 < exon_end_2 <= exon_end_1):
                the_master_exon_is = 1

            elif ( exon_start_2 <= exon_start_1 < exon_end_1 <= exon_end_2):
                the_master_exon_is = 2

                
            if ( the_master_exon_is > 0 ):
                # however, we may change our minds
                if  (exon_1.is_canonical and not exon_2.is_canonical):
                    the_master_exon_is = 1

                elif (exon_2.is_canonical and not exon_1.is_canonical):
                    the_master_exon_is = 2

                elif (is_havana[exon_1] and not is_havana[exon_2]):
                    the_master_exon_is = 1

                elif (is_havana[exon_2] and not is_havana[exon_1]):
                    the_master_exon_is = 2

                elif (is_ensembl[exon_1] and not is_ensembl[exon_2]):
                    the_master_exon_is = 1

                elif (is_ensembl[exon_2] and not is_ensembl[exon_1]):
                    the_master_exon_is = 2

                elif (exon_1.is_known and not exon_2.is_known):
                    the_master_exon_is = 1

                elif (exon_2.is_known and not exon_1.is_known):
                    the_master_exon_is = 2

            # the overlapping exons are all grouped
            # under the heading of the "master" version of the exon
            if   ( the_master_exon_is < 0):
                continue
            elif ( the_master_exon_is == 1):

                if (not overlapping_cluster.has_key(exon_1) ):
                    overlapping_cluster[exon_1] = []
                overlapping_cluster[exon_1].append(exon_2)
     
            elif ( the_master_exon_is == 2):

                if (not overlapping_cluster.has_key(exon_2) ):
                    overlapping_cluster[exon_2] = []
                overlapping_cluster[exon_2].append(exon_1)
        
        
        if   ( the_master_exon_is < 0):
            if (not overlapping_cluster.has_key(exon_1) ):
                    overlapping_cluster[exon_1] = []


    for exon_1 in overlapping_cluster.keys():
        exon_1.covering_exon       = -1 # nobody's covering this guy
        exon_1.covering_exon_known = -1 # formal

        for exon_2 in overlapping_cluster[exon_1]:
            exon_2.covering_exon        = exon_1.exon_id;
            exon_2.covering_exon_known  = exon_1.is_known;

            
#########################################
def mark_coding (cursor, gene_id, species, exons):

    for exon in exons:

        exon.is_coding = 0
        if ( exon.is_known):
            # find related transcripts
            # if there is a related trascript, the thing is coding
            # (the default is not coding)
            qry  = "select count(1) from exon_transcript "
            qry += " where exon_id = " % exon.exon_id
            rows = search_db (cursor, qry)
            if (not rows[0][0]):
                continue
            # now need to check that it is within the coding region
            # if yes
            exon.is_coding = 1
            # otherwise
            exon.is_coding = 0
          
        else: # exon belogs to a predicted transcript
            # get precition exon id
            # check we fall within the translated region



    return [coding_region_start, coding_region_end]


#########################################
def find_exons (cursor, gene_id, species):

    coding_region_start = -1
    coding_region_end   = -1
    exons               = []
    
    # get all exons from the 'exon' table
    exons = get_exons (cursor, gene_id, species, 'exon')
    # get all exons from the 'predicted_exon' table
    if (not species == 'homo_sapiens'):
        exons += get_exons (cursor, gene_id, species, 'prediction_exon')
    # mark the exons belonging to canonical transcript
    mark_canonical (cursor, gene_id, exons)
    # get annotation info
    fill_in_annotation_info (cursor, gene_id, exons)
    # find covering exons
    sort_out_covering_exons (cursor, exons)
    # mark coding exons
    ret = mark_coding (cursor, gene_id, species, exons)
    if ret:
        [coding_region_start, coding_region_end] = ret
    else:
        print "error retrieving coding transcripts for ", gene_id
        exit (1)


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
            print "\t", gene2stable(cursor, gene_id=gene_id)

            # find all exons associated with the gene id
            [coding_region_start, coding_region_end, exons] = \
                find_exons (cursor, gene_id, species)
            if (not exons):
                print  gene_id, " no exons found" 
                exit(1)

            print gene_id, " number of exons: ", len(exons)
            for exon in exons:
                print
                print exon
            break

            # store to gene2exon table


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
