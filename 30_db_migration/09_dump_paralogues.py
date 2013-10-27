#!/usr/bin/python


import MySQLdb, os
from   el_utils.mysql         import  *
from   el_utils.ensembl       import  *
from   el_utils.threads       import  parallelize
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.utils         import  erropen


########
def read_paralogues(cursor, gene_id):

    paralogues = []
    qry = "select gene.stable_id from gene, paralogue "
    qry += " where gene.gene_id = paralogue.cognate_gene_id "
    qry += " and paralogue.gene_id = %d " % gene_id
    rows = search_db (cursor, qry)


    if rows and not 'ERROR' in rows[0]:
        for row in rows:
            paralogues.append(row[0])
    
    return paralogues
        


#########################################
def dump_paralogues(species_list, db_info):
    
    [local_db, ensembl_db_name, outdir] = db_info
    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()


    for species in species_list:
        print
        print "############################"
        print  species
        qry = "use " + ensembl_db_name[species]
        search_db(cursor, qry)
        
        outfile  = "{0}/{1}_para_dump.txt".format(outdir, species)
        print outfile
        continue
        of       = erropen (outfile,"w")
        if not of: continue
        
        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        para_table = 'paralogue'

        ct   =  0
        seen = []
        for gene_id in gene_ids:
            ct += 1
            if not ct%100: print "\t", species, "   ", ct, "out of", len(gene_ids)

            if gene_id in seen: continue

            stable_id  = gene2stable(cursor, gene_id)
                
            paralogues = read_paralogues(cursor, gene_id)
            
            if ( paralogues):
                # dump
                for para in paralogues:
                    print >> of,  stable_id, para
                seen += paralogues
                
        of.close()
 
    cursor.close()
    db.close()

#########################################
def main():
    
    no_threads = 1
    local_db = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)


    outdir  = "{0}/para_dump".format(cfg.dir_path['afs_dumps'])
    print outdir
    if not os.path.exists(outdir):
        print outdir, "not found"
        exit(1) # exit after dir existence check
    cursor.close()
    db    .close()

    parallelize (no_threads, dump_paralogues, all_species, [local_db, ensembl_db_name, outdir])

    
    return True


#########################################
if __name__ == '__main__':
    main()
