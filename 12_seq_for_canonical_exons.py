#!/usr/bin/python

import MySQLdb
import commands
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.el_specific import  *
from   el_utils.exon    import  Exon
from   el_utils.threads import  parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.special_gene_sets  import get_theme_ids

# BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#def is_mitochondrial (cursor, seq_region_id):

    
########################################
def  get_canonical_exons (cursor, gene_id):

    exons = gene2exon_list (cursor, gene_id)
    if (not exons):
        #print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
        return []
        #exit(1) # shouldn't happen at this point

    # sorting exons in place by their start in gene:
    exons.sort(key=lambda exon: exon.start_in_gene)

    canonical_coding_exons = []
    reading = False
    for exon in exons:        
        if (not exon.is_canonical or  not exon.is_coding): 
             continue
        if (not exon.canon_transl_start is None):
            reading = True
        if (reading):
            canonical_coding_exons.append(exon)
        if (not exon.canon_transl_end is None):
            break  

    return canonical_coding_exons



#########################################
def reconstruct_exon_seqs (gene_seq, exons):

    """
    Return dna sequences belonging to each exon in the input list, according to Ensembl.
    """

    exon_seq    = {} 
    left_flank  = {}
    right_flank = {}

    for exon in exons:
        start   = exon.start_in_gene
        end     = exon.end_in_gene
        exon_id = exon.exon_id
        left_flank[exon_id]  = gene_seq[start-15:start]
        exon_seq[exon_id]    = gene_seq[start:end+1]
        right_flank[exon_id] = gene_seq[end+1:end+16]
    

    return [exon_seq, left_flank, right_flank]
            
#########################################
def store (cursor, exons, exon_seq, left_flank, right_flank, canonical_exon_pepseq):
    
    """
    Calls store_or_update() from el_utils.mysql.py:
    It sets the fixed fields to be exon_id, is_known, and is_sw;
    The rest of the exon info is updateable.
    """

    for exon in exons:
        exon_id = exon.exon_id
        #####
        fixed_fields  = {}
        update_fields = {}
        fixed_fields['exon_id']          = exon_id
        fixed_fields['is_known']         = exon.is_known
        fixed_fields['is_sw']            = 0
        if (exon_seq[exon_id]):
            update_fields['dna_seq']     = exon_seq[exon_id]
        if (left_flank[exon_id] ):
            update_fields['left_flank']  = left_flank[exon_id]
        if (right_flank[exon_id] ):
            update_fields['right_flank'] = right_flank[exon_id]
        if ( canonical_exon_pepseq.has_key(exon_id) and canonical_exon_pepseq[exon_id]):
            update_fields['protein_seq'] = canonical_exon_pepseq[exon_id]
        #####
        store_or_update (cursor, 'exon_seq', fixed_fields, update_fields)



#########################################
def store_exon_seqs(species_list, db_info):

    """
    The core of the operation: for each exon find and store dna sequence; for canonical exons also add the translation.
    For each species retrieves all protein coding genes, 
    and for each of the genes its full sequence and  the list of the known exons;  
    assigns dna  sequence to each exon (+ the translation for the canonical exons) and stores it in db.

    The reason why only the canonical exons get the translation at this point is that the canonical translation available in
    Ensembl is used to make sure we are looking at the right genome region. Namely, some sequences have a "patch"
    with a corrected or completed sequence. The info about that can be found in the assembly_exception table.
    Once we have the translation we can just as well store it.

    Note that the information about the translation start and end positions are added in the following script in the pipeline.

    """
    
    special    = ''

    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
        cfg    = ConfigurationReader()
    else:
        db     = connect_to_mysql          (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader       (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    for species in species_list[:len(species_list)/2]:
        print
        print "############################"
        print  species
        #continue
        if not switch_to_db(cursor, ensembl_db_name[species]):
            return False      
 
        if special:
            print "using", special, "set"
            gene_ids = get_theme_ids (cursor,  ensembl_db_name, cfg, special )
        elif (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        ###########################
        seqs_not_found = []
        ct  = 0
        tot = 0
        for gene_id in gene_ids:
            tot += 1
            if (not  tot%1000):
                print species, "tot genes:", tot, " fail:", ct
               
            # extract raw gene  region - bonus return from checking whether the 
            # sequence is correct: translation of canonical exons
            ret = get_gene_seq(acg, cursor, gene_id, species)
            [gene_seq, canonical_exon_pepseq, file_name, seq_name, seq_region_start, seq_region_end]  = ret

            if (not gene_seq or not canonical_exon_pepseq):
                ct += 1
                print 'no sequence found for ', gene_id, "   ",   ct, "out of ", tot
                seqs_not_found.append(gene_id)
                continue

            # get _all_ exons
            exons = gene2exon_list(cursor, gene_id, ensembl_db_name[species])
            if (not exons):
                #print 'no exons for ', gene_id
                #exit(1)
                ct += 1
                continue

            # get the sequence for each of the exons, as well as for the flanks
            # (the return are three dictionaries, with exon_ids as keys)
            [exon_seq, left_flank, right_flank] = reconstruct_exon_seqs (gene_seq, exons)
            store (cursor, exons, exon_seq, left_flank, right_flank, canonical_exon_pepseq)


        print species, "done; tot:", tot, " fail:", ct
        if (seqs_not_found):
            outf = open(species+".seqs_not_found", "w")
            for not_found in seqs_not_found:
                print >> outf, str(not_found)+", ",
            print >> outf, "\n"
            outf.close
    cursor.close()
    db    .close()



#########################################
def main():

    """
    Main entry point, but in reality does nothing except taking care of the parallelization.
    The parallelization here is per-species.
    """


    no_threads = 10

    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    cursor.close()
    db    .close()


    parallelize (no_threads, store_exon_seqs, all_species, [local_db, ensembl_db_name])





#########################################
if __name__ == '__main__':
    main()
'''
    # Ivana:
    # What is the proper way to find out whether the seq_region is mitochondrial?
    # EMily:
    # It's in the seq_region table, under name. Name will either be a chromosome
    # number, MT for mitochrondria or the name of a contig.
Hi Ivana

To find which contigs are mitochondrial, use the assembly table.

select s1.name from seq_region s1, assembly asm, seq_region s2 where
s1.seq_region_id = asm.cmp_seq_region_id and s2.seq_region_id =
asm.asm_seq_region_id and s2.name = "MT" ;
in human, returns
+----------------+
| name |
+----------------+
| NC_012920 |
| J01415.1.16569 |
+----------------+

Emily

Hi Ivana

The reason for this is that patches are not always the same length as the
genomic region they patch over. In most cases, a patch corrects sequencing
errors but the number of bp stays the same, but in some cases, a patch adds a
chunk of sequence. We anchor the 5' (wrt the chromosome orientation) end of the
patch to the identical reference coordinates, and allow the 3' end to differ
slightly from the reference coordinates. We don't change the complete genomic
coordinates every time a patch is added (this happens when we bring out a new
assembly) as this would be more hassle than it's worth and most people don't
notice it anyway. However the one thing it does affect is the coordinates of
genes that overlap the 3' end of the patch.

ENSG00000261899, and presumably the other genes you've had a similar issue
with, overlaps the 3' end of a patch, so its coordinates on the patch are
shifted compared to its coordinates on the reference genome.

Here's the data from the assembly_exception table for the particular patch that
affects ENSG00000261899:
+-----------------------+---------------+------------------+----------------+-------------+------------------+-----+
| assembly_exception_id | seq_region_id | seq_region_start | seq_region_end |exc_type | exc_seq_region_id | exc_seq_region_start | exc_seq_region_end | ori
|
+-----------------------+---------------+------------------+----------------+-------------+-------------------+----------------------+--------------------+-----+
| 67                    | 1000057054    | 36453102 | 36596491 | PATCH_NOVEL | 27508 | 36453102 |36590458 | 1 |
+-----------------------+---------------+------------------+----------------+-------------+-------------------+----------------------+--------------------+-----+

To get over this you need to link over the assembly_exception table. This will
give you the exact relationship between the patch you are looking at and the
chromosome it is linked to. If you then use the
seq_region_start/exc_seq_region_start relationship, this will cover all cases,
whether there is a shift or not.

Hope this helps,

Emily
'''
