#!/usr/bin/python

import MySQLdb
import sys, commands
from   el_utils.mysql         import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl       import  get_species, get_gene_ids, gene2stable, gene2exon_list, get_exon_seqs
from   el_utils.exon          import  Exon
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.threads       import  parallelize
from   el_utils.utils         import  erropen

#########################################
def  get_analysis_dict(cursor):
    source = {}
    qry  = "select analysis_id, logic_name  from analysis"
    rows = search_db (cursor, qry)
    if (not rows):
        print "blah?"
        return False
    for row in rows:
        source[row[0]] = row[1]
    return source

#########################################
def exon_tabstring(exon, stable_id, species, analysis, exon_seqs):
    
    ret  = ""
    ret += str(exon.exon_id) +"\t"
    ret += stable_id +"\t"
    ret += str(exon.start_in_gene)   +"\t"
    ret += str(exon.end_in_gene)     +"\t"
    ret += str(exon.strand)          +"\t"
    ret += str(exon.is_known)        +"\t"
    ret += str(exon.is_coding)       +"\t"
    ret += str(exon.is_canonical)    +"\t"
    ret += str(exon.is_constitutive) +"\t"
    ret += species                   +"\t"
    ret += analysis +"\t"
    ret += "\t".join(exon_seqs)

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
        #if (not species=='homo_sapiens'):
        #    continue
        outfile  = "{0}/{1}_exon_dump.txt".format(out_path, species)
        of       = erropen (outfile,"w")
        print of
        switch_to_db (cursor,  ensembl_db_name[species])

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        source = get_analysis_dict(cursor)

        tot    = 0
        ct     = 0
        for gene_id in gene_ids:
            tot += 1
            if (not  tot%1000):
                print species, ct, tot

            # get _all_ exons
            exons = gene2exon_list(cursor, gene_id)
            if (not exons):
                print 'no exons for ', gene_id
                sys.exit(1)

            for exon in exons:
                # exons seqs are its aa translation, left_flank, right_flank, and dna_seq
                exon_seqs = get_exon_seqs(cursor, exon.exon_id, exon.is_known)
                if (not exon_seqs):
                    ct += 1
                    continue
                # human readable string describing the source of annotation for this exon
                analysis  = source[exon.analysis_id] 
                # the first field return by get_exon_seqs is the exon_seq_id, so get rid of it
                print >> of, exon_tabstring (exon, gene2stable(cursor,gene_id), species, analysis, exon_seqs[1:])

        of.close()

    cursor.close()
    db    .close()

#########################################
def main():

    no_threads = 10

    local_db = False

    if local_db:
        db = connect_to_mysql()
    else:
        db = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    parallelize (no_threads, dump_exons, all_species, [local_db, ensembl_db_name])



#########################################
if __name__ == '__main__':
    main()
