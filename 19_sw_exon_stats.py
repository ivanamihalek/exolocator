#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl import  get_species, get_gene_ids


#########################################
def main():
    
    no_threads = 1
    local_db = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    tot = 0
    for species in all_species:

        switch_to_db(cursor, ensembl_db_name[species])
        qry  =  "select count(1) from sw_exon  join  exon_seq  "
        qry +=  "on sw_exon.exon_seq_id = exon_seq.exon_seq_id "
        qry +=  " where exon_seq.protein_seq  is not  null"
        rows = search_db(cursor, qry)
        
        if rows and rows[0][0]:
            print species,
            print "\t", rows[0][0]
            tot += int (rows[0][0])


    print
    print "total: ", tot


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
