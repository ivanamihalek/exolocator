#!/usr/bin/python -u

import MySQLdb
import sys, commands, os
from   el_utils.mysql      import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl       import  get_species, get_gene_ids, gene2stable, gene2exon_list, get_exon_seqs
from   el_utils.exon          import  Exon
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.threads       import  parallelize
from   el_utils.utils         import  erropen, mkdir_p

  
#########################################
def main ():

    
    db     = connect_to_mysql()
    cfg    = ConfigurationReader()
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    outpath = cfg.get_path('afs_dumps')
    outdir   = "{0}/exon_map".format(outpath)
    if (not os.path.exists(outdir)):
        mkdir_p(outdir)

    outfile  = "{0}/exon_map.sql".format(outdir)
    if os.path.exists('.creds'):
        [user, passwd, host, port] = read_creds()
    else:
        print "creds not found"
        exit(1)
    credentials = " -h {0} -P {1} -u {2}  -p{3}".format(host, port, user, password)
    cmd = "mysqldump {0} {1} exon_map > {2}".format (credentials, ensembl_db_name['homo_sapiens'], outfile)

    print cmd
    ret = commands.getoutput(cmd)
    
    print ret

    return True


#########################################
if __name__ == '__main__':
    main()
