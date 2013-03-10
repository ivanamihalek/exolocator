#!/usr/bin/python

import MySQLdb
import sys, commands
from   el_utils.mysql         import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl       import  *
from   el_utils.exon          import  Exon
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.threads       import  parallelize
from   el_utils.utils         import  erropen


#########################################
def exon_tabstring(exon, gene_stable_id, exon_stable_id, human_exon_stable_id,  species, analysis, exon_seqs):
    
    ret  = ""
    ret += str(exon.exon_id) +"\t"
    ret += gene_stable_id +"\t"
    ret += exon_stable_id +"\t"
    ret += str(exon.start_in_gene)   +"\t"
    ret += str(exon.end_in_gene)     +"\t"
    ret += str(exon.strand)          +"\t"
    ret += str(exon.is_known)        +"\t"
    ret += str(exon.is_coding)       +"\t"
    ret += str(exon.is_canonical)    +"\t"
    ret += str(exon.is_constitutive) +"\t"
    ret += species                   +"\t"
    ret += analysis +"\t"
    if analysis=='sw_sharp':
         ret += human_exon_stable_id +"\t" 

    ret += "\t".join( map((lambda token:  type(token) is str and token or str(token)), exon_seqs) )

    return ret

  
#########################################
def dump_exons (species_list, db_info):

    
    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        cfg      = ConfigurationReader()
    else:
        db     = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor   = db.cursor()

    out_path = cfg.get_path('afs_dumps')

    for species in species_list:

        if  species=='homo_sapiens':  continue

        outfile  = "{0}/{1}_sw_exons.txt".format(out_path, species)
        of       = erropen (outfile,"w")
        print outfile

        analysis = 'sw_sharp'

        switch_to_db (cursor,  ensembl_db_name[species])
        qry      = "select * from sw_exon"
        sw_exons = search_db(cursor, qry)

        if not sw_exons:
            continue


        for sw_exon in sw_exons:

            [sw_exon_id, gene_id, start_in_gene, end_in_gene, human_exon_id,
             exon_seq_id, template_exon_id, template_species, strand, phase, has_NNN, has_stop, has_3p_ss, has_5p_ss] = sw_exon

            if has_stop: continue

            human_coding = is_coding (cursor, human_exon_id, ensembl_db_name['homo_sapiens'])
            human_exon_stable_id  = exon2stable(cursor, human_exon_id)

            exon_seqs = get_exon_seq_by_db_id (cursor, exon_seq_id, ensembl_db_name[species])
            if (not exon_seqs):
                continue
            # exons seqs are its aa translation, left_flank, right_flank, and dna_seq
            [exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, 
             left_flank, right_flank, dna_seq] = exon_seqs

            # the first field return by get_exon_seqs is the exon_seq_id, so get rid of it
            gene_stable_id = gene2stable(cursor,gene_id)
            exon_stable_id = "anon"


            exon = Exon()

            exon.exon_id = sw_exon_id
            exon.start_in_gene = start_in_gene
            exon.end_in_gene   = end_in_gene
            exon.strand = strand
            exon.is_known        = 0
            exon.is_canonical    = 0
            exon.is_constitutive = 0
            exon.is_coding       = 1 if human_coding else 0
           

            print >> of, exon_tabstring (exon, gene_stable_id, exon_stable_id, human_exon_stable_id, species, analysis, exon_seqs[1:])


        of.close()
    
    cursor.close()
    db    .close()

#########################################
def main():

    no_threads = 1

    local_db = False

    if local_db:
        db = connect_to_mysql()
    else:
        db = connect_to_mysql (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    parallelize (no_threads, dump_exons, all_species, [local_db, ensembl_db_name])



#########################################
if __name__ == '__main__':
    main()
