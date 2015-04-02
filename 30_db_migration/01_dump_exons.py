#!/usr/bin/python -u

import MySQLdb
import sys, commands, os
from   el_utils.mysql         import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl       import  *
from   el_utils.exon          import  Exon
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.threads       import  parallelize
from   el_utils.utils         import  erropen


#########################################
def exon_tabstring(exon, gene_stable_id, exon_stable_id, species, analysis, exon_seqs):
    
    ret  = ""
    ret += str(exon.exon_id) +"\t"
    ret += gene_stable_id    +"\t"
    ret += exon_stable_id    +"\t"
    ret += str(exon.start_in_gene)   +"\t"
    ret += str(exon.end_in_gene)     +"\t"
    ret += str(exon.strand)          +"\t"
    ret += str(exon.is_known)        +"\t"
    ret += str(exon.is_coding)       +"\t"
    ret += str(exon.is_canonical)    +"\t"
    ret += str(exon.is_constitutive) +"\t"
    ret += species                   +"\t"
    ret += analysis +"\t"


    ret += "\t".join( map((lambda token:  type(token) is str and token or str(token)), exon_seqs) )

    return ret

  
#########################################
def dump_exons (species_list, db_info):

    
    [local_db, ensembl_db_name] = db_info
    db     = connect_to_mysql()
    cfg    = ConfigurationReader()
    cursor = db.cursor()

    out_path = "{0}/exons".format(cfg.get_path('afs_dumps'))
    if not os.path.exists(out_path):
        print out_path, "not found"
        exit (1) # exit on failed output dir check

    for species in species_list:
        #if (not species=='homo_sapiens'):
        #    continue
        outfile  = "{0}/{1}_exon_dump.txt".format(out_path, species)
        of       = erropen (outfile,"w")
        if not of:  continue
        switch_to_db (cursor,  ensembl_db_name[species])

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1, ref_only=True)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        source = get_analysis_dict(cursor)

        ct     = 0
        for gene_id in gene_ids:
            ct += 1
            if (not  ct%1000):
                print species, ct, len(gene_ids)

            # get _all_ exons
            exons = gene2exon_list(cursor, gene_id)
            if (not exons):
                print 'no exons for ', gene_id
                continue

            for exon in exons:

                if exon.covering_exon  > 0: continue
                # exons seqs are its aa translation, left_flank, right_flank, and dna_seq
                exon_seqs = get_exon_seqs(cursor, exon.exon_id, exon.is_known)
                if (not exon_seqs):
                    continue
                # human readable string describing the source of annotation for this exon
                if exon.is_known==2:
                    analysis = 'sw_sharp'
                elif exon.is_known==3:
                    analysis = 'usearch'
                else:
                    analysis = source[exon.analysis_id] 
                # the first field return by get_exon_seqs is the exon_seq_id, so get rid of it
                gene_stable_id = gene2stable(cursor,gene_id)
                if ( exon.is_known == 1):
                    exon_stable_id = exon2stable(cursor,exon.exon_id)
                elif ( exon.is_known == 2):
                    exon_stable_id = 'sw_sharp_'+str(exon.exon_id)
                elif ( exon.is_known == 3):
                    exon_stable_id = 'usearch_'+str(exon.exon_id)
                else:
                    exon_stable_id = "anon"

                print >> of, exon_tabstring (exon, gene_stable_id, exon_stable_id, species, analysis, exon_seqs[1:])


        of.close()
        print species, "done"
    
    cursor.close()
    db    .close()

#########################################
def main():

    no_threads = 10

    db = connect_to_mysql()
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    parallelize (no_threads, dump_exons, all_species, [local_db, ensembl_db_name])



#########################################
if __name__ == '__main__':
    main()
