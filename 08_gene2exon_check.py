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
def inspect (exons, canonical_translation):


    total_len = 0

    for exon in exons:

        exon_len = 0
        if  exon.canon_transl_end is None:
            exon_len = exon.end_in_gene - exon.start_in_gene + 1
        else:
            exon_len = exon.canon_transl_end + 1
        if not exon.canon_transl_start is None:
            exon_len -= exon.canon_transl_start
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
        # print canonical sequence with the newline stuck in every 50 positions
        print re.sub("(.{50})", "\\1\n", canonical_translation)       

        # where did the information about translation start and end come at this point?

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
        for tot in range(1000):

            gene_id = choice(gene_ids)

            tot += 1
            # find all exons associated with the gene id
            exons = filter (lambda x: x.is_coding and x.is_canonical, gene2exon_list (cursor, gene_id))
            if (not exons):
                ct +=1
                print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
                
            if not tot%100:
                print species, tot, ct

            # add up the coding length of the canonical exons
            length = 0
            for exon in exons:

                if  exon.canon_transl_end is None:
                    length += exon.end_in_gene - exon.start_in_gene + 1
                else:
                    length += exon.canon_transl_end + 1

                if not exon.canon_transl_start is None:
                    length -= exon.canon_transl_start
            
            if (not length):
                print gene2stable (cursor, gene_id = gene_id), " no exons marked as canonical"
                continue

            # what is the length of the canonical transcript according to Ensembl
            canonical_translation = get_canonical_transl (acg, cursor, gene_id, species)
            if ( not canonical_translation):
                print "no canonical transl found for ", gene_id
                continue


            if ( abs(length/3-len(canonical_translation)) > 3):
                ct +=1
                print gene2stable (cursor, gene_id = gene_id),
                print "(length of all exons)/3 ", length/3, 
                print " does not match reported canonical transl len ", len(canonical_translation)
                inspect (exons, canonical_translation)
                exit(1)

        print species, "checked a sample of ", tot, "genes;  problematic:", ct

    cursor.close()
    db.close()

    return True




#########################################
if __name__ == '__main__':
    main()
