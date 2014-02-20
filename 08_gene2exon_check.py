#!/usr/bin/python

import MySQLdb
import commands
from   random               import choice
from   el_utils.mysql       import  connect_to_mysql, search_db
from   el_utils.ensembl     import  *
from   el_utils.el_specific import  *
from   el_utils.exon        import  Exon
from   el_utils.threads     import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator

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
        print "start not found for ", exon_id
        return None

    return rows[0][0]

#########################################
def get_translated_region_talkative(cursor, gene_id, species):

    # get the region on the gene
    is_known = (species == 'homo_sapiens')
    ret = get_gene_region (cursor, gene_id, is_known)
    if  ret:
        [gene_seq_id,gene_region_start, gene_region_end, 
         gene_region_strand] = ret
    else:
        print "region not retrived for ", species, gene_id, species
        return []

    canonical_transcript_id = get_canonical_transcript_id (cursor, gene_id)
    transcript_ids = get_transcript_ids(cursor, gene_id)
    print transcript_ids
    print "canonical: ", canonical_transcript_id

    transl_region_start = gene_region_end
    transl_region_end   = gene_region_start

    print "transl region start:", transl_region_start
    print "transl region end:", transl_region_end

    for[ transcript_id, transcript_stable] in transcript_ids:
   
        qry  = "SELECT seq_start, start_exon_id, seq_end, end_exon_id " 
        qry += " FROM translation WHERE transcript_id=%d"  %  transcript_id
        rows = search_db (cursor, qry)
        if (not rows):
            continue
        exon_seq_start = rows[0][0]
        start_exon_id  = rows[0][1]
        exon_seq_end   = rows[0][2]
        end_exon_id    = rows[0][3]
        
        print
        if transcript_id == canonical_transcript_id:
            print "canonical: "
        print "transcript id: ", transcript_id
        print "start exon id:",  start_exon_id, "transl start (in the exon) ", exon_seq_start
        print "end exon id:",    end_exon_id,   "transl end (in the exon)", exon_seq_end
        

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

    return


#########################################
def inspect (exons):


    total_len = 0

    for exon in exons:

        exon_len = 0
        if  exon.canon_transl_end is None:
            exon_len = exon.end_in_gene - exon.start_in_gene + 1
        else:
            exon_len = exon.canon_transl_end

        if not exon.canon_transl_start is None:
            exon_len -= exon.canon_transl_start - 1

        total_len += exon_len

        print "*****"
        print "exon id: ",  exon.exon_id         
        print "canonical:", exon.is_canonical
        print "coding:",    exon.is_coding
        print "start in gene: ", exon.start_in_gene 
        print "end in gene: ",   exon.end_in_gene 
        print "canon transl start: ", exon.canon_transl_start
        print "canon transl end: ",   exon.canon_transl_end
        print "exon len %5d   total %d " % (exon_len, total_len)

    print

#########################################
def main():

    local_db = False

    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    if  len(sys.argv) > 1:
        species_list = sys.argv[1:]
    else:
        species_list = all_species
 
    species_list = ['homo_sapiens']
    ############################
    for species in species_list:
        print
        print "############################"
        print  species

        switch_to_db (cursor, ensembl_db_name[species])

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        
        ct = 0
        tot = 0

        #for tot in range(1000):
        #for gene_id in gene_ids:
        for gene_id in [727579]:
            tot += 1
            #gene_id = choice(gene_ids)
            # find all canonical coding exons associated with the gene id
            exons = get_canonical_coding_exons (cursor, gene_id)
            if (not exons):
                ct +=1
                print gene_id, gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
                
            if not tot%100:
                print species, tot, ct

            # add up the coding length of the canonical exons
            exons_sorted = exons.sort(key=lambda exon: exon.start_in_gene)

            inside_the_coding_range = False
            start_properly_marked   = False
            length = 0
            for exon in exons_sorted:

                if not exon.canon_transl_start is None:
                    start_properly_marked   = True # if it is not propermy marked, we'll never start reading
                    inside_the_coding_range = True
                    length -= exon.canon_transl_start - 1

                if not exon.canon_transl_start is None:
                    inside_the_coding_range = True
                    length += exon.canon_transl_end 

                if inside_the_coding_range:
                    length += exon.end_in_gene - exon.start_in_gene + 1

            # take that all exons are coding full length if there is no start and end annotation
            # (this I believe is the case for predicted transcripts)
            if not start_properly_marked:
                length = 0
                for exon in exons_sorted:
                    length += exon.end_in_gene - exon.start_in_gene + 1
                
            
            if (not length):
                print gene2stable (cursor, gene_id = gene_id), " no exons marked as canonical"
                continue

            # what is the length of the canonical transcript according to Ensembl
            canonical_translation = get_canonical_transl (acg, cursor, gene_id, species)
            if ( not canonical_translation):
                print "no canonical transl found for ", gene_id
                exit(1)
                continue


            if ( abs(length/3-len(canonical_translation)) > 3):
                ct +=1
                print gene_id, gene2stable (cursor, gene_id), get_description (cursor, gene_id)
                print "(length of all exons)/3 ", length/3, 
                print " does not match reported canonical transl len ", len(canonical_translation)
                # print out all exons
                print "exons:"
                inspect (exons)
                print re.sub("(.{50})", "\\1\n", canonical_translation)  # print canonical sequence with \n stuck in every 50 positions     
                print
                # print out exons more carefully filtered to belong to the canonical version of the translation
                print
                get_translated_region_talkative (cursor, gene_id, species)
                all_exons =  gene2exon_list (cursor, gene_id)
                print "all exons:"
                inspect (all_exons)

                exit(1)

        print species, "checked a sample of ", tot, "genes;  problematic:", ct

    cursor.close()
    db.close()

    return True




#########################################
if __name__ == '__main__':
    main()
