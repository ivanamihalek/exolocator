#!/usr/bin/python

import MySQLdb
import commands
from   random import choice
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.ensembl import  *
from   el_utils.exon    import  Exon
from   el_utils.threads import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator




            


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

    species_list = all_species
 
    for species in species_list:
        print
        print "############################"
        print  species
        #if (species in ['ailuropoda_melanoleuca', 'anolis_carolinensis', 
        #                'bos_taurus','danio_rerio']):
        #    continue
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
            exons = gene2exon_list (cursor, gene_id)
            if (not exons):
                ct +=1
                print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
                
            if not tot%100:
                print species, tot, ct

            length = 0
            for exon in exons:
                if (not exon.is_canonical): 
                    continue
                
                if (not exon.is_coding): 
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
            cmd = acg.generate_fastacmd_protein_command (canonical_transl_id, species,  "all", None)
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
                print "exon length ", length/3, 
                print " does not match reported canonical transl len ", translation_length
            
        print species, "checked a sample of ", tot, "genes;  problematic:", ct

    cursor.close()
    db.close()

    return True




#########################################
if __name__ == '__main__':
    main()
