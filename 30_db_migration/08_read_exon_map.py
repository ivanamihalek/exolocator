#!/usr/bin/python

import MySQLdb
import sys, commands, os
from   el_utils.mysql      import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl       import  get_species, get_gene_ids, gene2stable, gene2exon_list, get_exon_seqs
from   el_utils.exon          import  Exon
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.threads       import  parallelize

  
#########################################
def main ():

    
    db_name = "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)

    cfg     = ConfigurationReader (user="marioot", passwd="tooiram", check=False)

    inpath = cfg.get_path('afs_dumps')
    indir   = "%s/exon_map"     % inpath
    infile  = "%s/exon_map.sql" % indir
    if (not os.path.exists(infile)):
        print "not found: ", infile
        sys.exit(1)
    print "reading", infile

    #credentials = " -u marioot -ptooiram"
    #cmd = "mysql %s  exolocator_db  <  %s" % (credentials, infile)

    qry = "drop table exon_map"
    rows = search_db(cursor, qry)

    qry = "source %s" % infile
    print qry
    rows = search_db(cursor, qry, verbose=True)
    print rows

    return True


#########################################
if __name__ == '__main__':
    main()
