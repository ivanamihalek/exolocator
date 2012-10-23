#!/usr/bin/python


import MySQLdb
from   el_utils.mysql         import  connect_to_mysql, search_db
from   el_utils.ensembl       import  get_species

####################################################
def get_seq_region_info(cursor, name):
    qry = "select * from seq_region where name = '%s'" % name
    rows = search_db (cursor, qry)
    if(len(rows) > 1):
        print "more than one entry associated with ", name
        exit (1)
    return rows[0]

####################################################
def main():
    db     = connect_to_mysql()
    cursor = db.cursor()


    [all_species, ensembl_db_name] = get_species (cursor)


    for species in ['danio_rerio']:
    #for species in all_species:
        print species

        qry = "use "+ensembl_db_name[species]
        search_db (cursor, qry)

        qry  = "select seq_region.name, seq_region.file_name from seq_region, gene "
        qry += " where gene.biotype='protein_coding' and gene.seq_region_id =  seq_region.seq_region_id "
            

        rows = search_db (cursor, qry)
        if (not rows):
            print "\t no seq region info found "
            continue
        tot = 0
        no_file = 0
        for row in rows:
            [name,  file_name] = row
            #print name, file_name
            tot += 1
            if (not file_name):
                no_file += 1
                print name, file_name
                #exit (1)

        print "\t tot seq_regions: ", tot, " no file: ", no_file
 
    cursor.close()
    db    .close()
    


####################################################
if __name__ == '__main__':
    main()
