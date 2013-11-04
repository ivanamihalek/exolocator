#!/usr/bin/python

import MySQLdb
import commands
from   random import choice
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  *
from   el_utils.exon    import  Exon
from   el_utils.threads import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.special_gene_sets  import  get_theme_ids
from   el_utils.config_reader      import ConfigurationReader

#########################################
def get_gene_region (cursor, gene_id, is_known=None):

    qry     = "SELECT seq_region_id, seq_region_start, seq_region_end, "
    qry    += " seq_region_strand "
    qry    += " FROM gene WHERE gene_id=%d"  %  gene_id
    if (not is_known is None and is_known):
        qry  += " and  status='known' "
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
def get_known_exons (cursor, gene_id, species):

    exons = []

    qry  = "select distinct exon_transcript.exon_id from  exon_transcript, transcript "
    qry += " where exon_transcript.transcript_id = transcript.transcript_id "
    qry += " and transcript.gene_id = %d " % gene_id

    rows = search_db (cursor, qry)
    
    if (not rows ):
        return []
    if ('Error' in rows[0]):
        search_db (cursor, qry, verbose = True)
        return []

    # get the region on the gene
    ret = get_gene_region (cursor, gene_id)
    if  ret:
        [gene_seq_id, gene_region_start, gene_region_end, 
         gene_region_strand] = ret
    else:
        print "region not retrived for ", species, gene_id
        return []

    exon_ids = []
    for row in rows:
        exon_ids.append(row[0])

    for exon_id in exon_ids:
        qry = "select * from exon where exon_id=%d" % exon_id
        rows = search_db (cursor, qry)
        if (not rows or 'Error' in rows[0]):
            search_db (cursor, qry, verbose = True)
            continue
        exon         = Exon()
        exon.gene_id = gene_id
        exon.load_from_ensembl_exon (gene_region_start, gene_region_end, rows[0])
        exons.append(exon)

    return exons


#########################################
def get_predicted_exons (cursor, gene_id, species):

    exons = []

    # get the region on the gene
    ret = get_gene_region (cursor, gene_id)
    if  ret:
        [gene_seq_id, gene_region_start, gene_region_end, 
         gene_region_strand] = ret
    else:
        print "region not retrived for ", species, gene_id
        return []

    qry    = "SELECT  * FROM  prediction_exon  WHERE seq_region_id = %d "  %  gene_seq_id
    qry   += " AND  seq_region_start >= %d AND seq_region_start <= %d " %  \
        (gene_region_start, gene_region_end)
    qry   += " AND  seq_region_end   >= %d AND seq_region_end   <= %d " %  \
        (gene_region_start, gene_region_end)
    rows   = search_db (cursor, qry)

    if (not rows):
        return []
    for row in rows:
        exon         = Exon()
        exon.gene_id = gene_id
        exon.load_from_ensembl_prediction (gene_region_start, gene_region_end, row)
        exons.append(exon)
 
    return exons


#########################################
def get_exons (cursor, gene_id, species, table):

    if ( table == 'exon'):
        return get_known_exons (cursor, gene_id, species)
    else:
        return get_predicted_exons (cursor, gene_id, species)



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
def get_canonical_coordinates (cursor, canonical_transcript_id):
    qry = "select seq_start, start_exon_id,  seq_end, end_exon_id "
    qry += " from translation where transcript_id = %d " % canonical_transcript_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         return []
    return rows[0]

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
    [canonical_start_in_exon, canonical_start_exon_id,
     canonical_end_in_exon, canonical_end_exon_id] = get_canonical_coordinates (cursor, canonical_transcript_id)
    
    start_found = False
    end_found   = False
    for exon in exons:
        if (exon.exon_id == canonical_start_exon_id):
            start_found = True
            exon.canon_transl_start = canonical_start_in_exon-1
            canonical_start_in_gene = exon.start_in_gene+exon.canon_transl_start
        if (exon.exon_id == canonical_end_exon_id):
            end_found = True
            exon.canon_transl_end = canonical_end_in_exon-1
            canonical_end_in_gene = exon.start_in_gene+exon.canon_transl_end
    if ( not start_found ):
        print "canonical translation start not found for ", gene_id
        exit(1)
    if ( not end_found ):
        print "canonical translation end not found for ", gene_id
        exit(1)

    # for each exon in canonical 
    qry = "select exon_id  from exon_transcript where transcript_id= %d" % canonical_transcript_id
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, veerbose=True)
        exit(1)

    canonical_ids = []
    for row in rows:
        canonical_ids.append(row[0])

    for exon in exons:
       exon.is_canonical = 0 # default
       if ( not exon.is_known):
           continue
       if (exon.exon_id in canonical_ids):
           exon.is_canonical = 1

         
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



#########################################
def get_transcript_ids(cursor, gene_id, species):

    rows = []
    if ( species=='homo_sapiens'):
        qry    = "SELECT transcript_id  FROM transcript "
        qry   += " WHERE gene_id=%d AND status = 'known' AND biotype='protein_coding' "  \
            %  gene_id
        rows   = search_db (cursor, qry, verbose=False)

    if ((not rows) or species!='homo_sapiens'): # try softer criteria
        qry    = "SELECT transcript_id  FROM transcript "
        qry   += " WHERE gene_id=%d AND biotype='protein_coding' "  %  gene_id
        rows   = search_db (cursor, qry, verbose=False)
  
    if (not rows):
        return []

    transcript_ids = []
    for row in rows:
        transcript_ids.append(row[0])

    return  transcript_ids

#########################################
def get_exon_start(cursor, exon_id):

    qry  = "select seq_region_start from exon "
    qry += "where exon_id = %d " % exon_id
    rows = search_db (cursor, qry)
    if (not rows or 'Error' in rows[0]):
        print "start not found gor ", exon_id
        return None

    return rows[0][0]


#########################################
def get_exon_end(cursor, exon_id):

 
    qry  = "select seq_region_end from exon "
    qry += "where exon_id = %d " % exon_id
    rows = search_db (cursor, qry)
    if (not rows or 'Error' in rows[0]):
        print "start not found gor ", exon_id
        return None

    return rows[0][0]


#########################################
def get_translated_region(cursor, gene_id, species):

    # get the region on the gene
    is_known = (species == 'homo_sapiens')
    ret = get_gene_region (cursor, gene_id, is_known)
    if  ret:
        [gene_seq_id,gene_region_start, gene_region_end, 
         gene_region_strand] = ret
    else:
        print "region not retrived for ", species, gene_id, species
        return []


    transcript_ids = get_transcript_ids(cursor, gene_id, species)


    transl_region_start = gene_region_end
    transl_region_end   = gene_region_start

    for transcript_id in transcript_ids:
   
        qry  = "SELECT seq_start, start_exon_id, seq_end, end_exon_id " 
        qry += " FROM translation WHERE transcript_id=%d"  %  transcript_id
        rows = search_db (cursor, qry)
        if (not rows):
            continue
        exon_seq_start = rows[0][0]
        start_exon_id  = rows[0][1]
        exon_seq_end   = rows[0][2]
        end_exon_id    = rows[0][3]

        if (gene_region_strand > 0):
            start = {}
 
            start[start_exon_id] = get_exon_start(cursor, start_exon_id)
            start[end_exon_id]   = get_exon_start(cursor, end_exon_id)

            this_translation_region_start = start[start_exon_id] + exon_seq_start - 1
            this_translation_region_end   = start[end_exon_id]   + exon_seq_end   - 1

        else: 
            end   = {}  

            end[start_exon_id] = get_exon_end (cursor, start_exon_id)
            end[end_exon_id]   = get_exon_end (cursor, end_exon_id)

            this_translation_region_start = end[end_exon_id]   - exon_seq_end   + 1
            this_translation_region_end   = end[start_exon_id] - exon_seq_start + 1

        if (this_translation_region_start <= transl_region_start):
            transl_region_start = this_translation_region_start
        
        if (this_translation_region_end >= transl_region_end):
            transl_region_end = this_translation_region_end

    return [transl_region_start, transl_region_end, gene_region_strand]

#########################################################
# resolve which exons are "master" among the ones
# which are havana, ensembl, predicted ....

from itertools  import combinations
from pygraph    import digraph
from pygraph.algorithms.cycles import find_cycle
from pygraph.algorithms.accessibility import accessibility

#########################################
def find_master (exon_1, exon_2, is_ensembl, is_havana):

    master_exon    = None
    covered_exon   = None

    havana_exon    = None
    ensembl_exon   = None
    canonical_exon = None
    superset_exon  = None
    known_exon     = None

    exon_start_1 = exon_1.start_in_gene
    exon_end_1   = exon_1.end_in_gene
    exon_start_2 = exon_2.start_in_gene
    exon_end_2   = exon_2.end_in_gene

    if (exon_start_1 > exon_end_2 or exon_start_2 > exon_end_1):
        return None, None # Fully disjoint exons

    if (exon_start_1 <= exon_start_2 < exon_end_2 <= exon_end_1):
        superset_exon = exon_1
    elif ( exon_start_2 <= exon_start_1 < exon_end_1 <= exon_end_2):
        superset_exon = exon_2

    if   (exon_1.is_canonical and not exon_2.is_canonical):
        canonical_exon = exon_1
    elif (exon_2.is_canonical and not exon_1.is_canonical):
        canonical_exon = exon_2

    if   (is_havana[exon_1] and not is_havana[exon_2]):
        havana_exon = exon_1
    elif (is_havana[exon_2] and not is_havana[exon_1]):
        havana_exon = exon_2

    if (is_ensembl[exon_1] and not is_ensembl[exon_2]):
        ensembl_exon = exon_1
    elif (is_ensembl[exon_2] and not is_ensembl[exon_1]):
        ensembl_exon = exon_2

    if (exon_1.is_known and not exon_2.is_known):
        known_exon = exon_1
    elif (exon_2.is_known and not exon_1.is_known):
        known_exon = exon_2

    if havana_exon is not None:
        master_exon = havana_exon
    elif canonical_exon is not None:
        master_exon = canonical_exon
    elif ensembl_exon is not None:
        master_exon = ensembl_exon
    elif known_exon is not None:
        master_exon = known_exon
    elif superset_exon is not None:
        master_exon = superset_exon

    if (master_exon == exon_1):
        covered_exon = exon_2
    else:
        covered_exon = exon_1

    return master_exon, covered_exon
    
#########################################
def sort_out_covering_exons (cursor, exons):

    # havana is manually curated and gets priority
    is_ensembl = {}
    is_havana  = {}
    for exon in exons:
        logic_name = get_logic_name(cursor, exon.analysis_id)
        is_ensembl[exon] = ('ensembl' in logic_name)
        is_havana [exon] = ('havana'  in logic_name)

    dg = digraph()
    dg.add_nodes(exons)
    for exon1, exon2 in combinations(dg.nodes(),2):
        master, covered = find_master(exon1,exon2,is_ensembl,is_havana)
        if master is not None and covered is not None:
            dg.add_edge(master,covered)
    assert not find_cycle(dg)
    clusters = dict(((k,v) for k,v in accessibility(dg).iteritems()
                     if not dg.incidents(k)))
    for k in clusters:
        clusters[k].remove(k)

    for master_exon, covered_list in clusters.iteritems():
        master_exon.covering_exon       = -1 # nobody's covering this guy
        master_exon.covering_exon_known = -1 # formal

        for covered_exon in covered_list:
            covered_exon.covering_exon       = master_exon.exon_id;
            covered_exon.covering_exon_known = master_exon.is_known;

           
#########################################
def mark_coding (cursor, gene_id, species, exons):

    ret = get_translated_region(cursor, gene_id, species)
    if ( not ret ):
        return False

    [transl_region_start, transl_region_end, strand] = ret
        
    translated_length = 0
    for exon in exons:

        exon.is_coding = 0
        if (exon.is_known):
            # inocent until proven guilty
            exon.is_coding = 0
            # find related transcripts
            # if there is a related trascript, the thing is coding
            # (the default is 'not coding')
            qry  = "select count(1) from exon_transcript "
            qry += " where  exon_id = %d " % exon.exon_id
            rows = search_db (cursor, qry)
            if (not rows[0][0]):
                continue

            # now need to check that it is within the coding region
            exon_start = get_exon_start (cursor, exon.exon_id)
            exon_end   = get_exon_end   (cursor, exon.exon_id)

            translated_length += exon_end-exon_start+1

            if ( exon_end < transl_region_start or
                 transl_region_end   < exon_start):
                exon.is_coding = 0

            else:
                # there is _a_ translation that covers this exon
                # otherwise I could have two disjunct translations from
                # the saem gene, and a couple of exons in the middle that 
                # are never translated - is that possible?
                exon.is_coding = 1

        else: # exons belongs to a  predicted transcript 
            # == we don't know if it is coding or not
            exon.is_coding = None
            

    return True


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
    mark_coding (cursor, gene_id, species, exons)


    return exons

#########################################
def gene2exon(species_list, db_info):

    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
        cfg    = ConfigurationReader()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    for species in species_list:

        print
        print "############################"
        print  species
        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        
        #gene_ids = get_theme_ids(cursor, ensembl_db_name, cfg, 'missing_seq')

        no_exons_found  = 0
        tot = 0
        status_not_known = 0
        biotype_not_protein_coding = 0
        seen = {}
        # check all genes (for one species, for example)
        #for gene_id in gene_ids:
        # check a random sample of genes
        for tot in range(50):
            gene_id = choice(gene_ids)

            if seen.has_key(gene_id):
                continue
            seen[gene_id] = True

            #tot += 1
            status  = get_status (cursor, gene_id)
            if not 'KNOWN' in status: 
                status_not_known += 1
                #print status
                continue

            biotype = get_biotype (cursor, gene_id)
            if not biotype == 'protein_coding': 
                #print biotype
                biotype_not_protein_coding += 1
                continue

            #print 
            #print gene_id, gene2stable (cursor, gene_id), biotype
            # find all exons associated with that gene id
            exons = find_exons (cursor, gene_id, species)
            if (not exons):
                no_exons_found +=1
                #print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
                
            continue
   
            length = 0
            for exon in exons:
                if (not exon.is_canonical or not exon.is_coding): 
                    # in some genomes an exon appears as canonical
                    # even though it is not coding
                    continue
 
                if (exon.canon_transl_end is None):
                    length += exon.end_in_gene - exon.start_in_gene + 1
                else:
                    length += exon.canon_transl_end + 1

                if (not exon.canon_transl_start is None):
                    length -= exon.canon_transl_start 
            
            if (not length):
                print gene2stable (cursor, gene_id = gene_id), " no exons marked as canonical"
                continue

                
            # what is the length of the canonical transcript according to Ensembl
            canonical_transl_id = gene2stable_canon_transl(cursor, gene_id)
            if ( not canonical_transl_id):
                print "no canonical transl id found for ", gene_id
                continue


            # get canonical transcript from ensembl fasta database
            cmd = acg.generate_fastacmd_protein_command (canonical_transl_id, species, 
                                                         "all", None)
            fasta = commands.getoutput(cmd)
            
            if (not fasta):
                print gene2stable (cursor, gene_id = gene_id), "fasta not found for ", canonical_transl_id
                continue
            
            translation_length = 0
            for line in fasta.split('\n'):
                if '>' in line:
                    continue
                translation_length += len(line)


            if ( abs(length/3-translation_length) > 3):
                ct +=1
                print gene2stable (cursor, gene_id = gene_id),
                print "exon length ", length/3, " transl len ", translation_length,
                print "(", ct, " out of ", tot, ")"
            

        print "total:", tot
        print "status_not_known: ",  status_not_known
        print "biotype_not_protein_coding: ", biotype_not_protein_coding
        print "no_exons_found: ", no_exons_found

    cursor.close()
    db.close()

    return True



#########################################
def main():

    no_threads = 1

    local_db = False

   
    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    parallelize (no_threads, gene2exon, all_species, [local_db, ensembl_db_name])



#########################################
if __name__ == '__main__':
    main()