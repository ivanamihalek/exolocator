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

    
    local_db = False

    if local_db:
        db     = connect_to_mysql()
        cfg      = ConfigurationReader()
    else:
        db     = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor   = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    inpath = cfg.get_path('afs_dumps')
    indir   = "%s/exon_map"     % inpath
    infile  = "%s/exon_map.sql" % indir
    if (not os.path.exists(infile)):
        print "not found: ", infile
        sys.exit(1)

    credentials = " -u marioot -ptooiram"
    cmd = "mysql %s  exolocator_db  <  %s" % (credentials, infile)

    print cmd
    ret = commands.getoutput(cmd)
    
    print ret

    return True


#########################################
if __name__ == '__main__':
    main()
