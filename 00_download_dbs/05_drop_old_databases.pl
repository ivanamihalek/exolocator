#!/usr/bin/python
import MySQLdb
import os, sys
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.mysql   import  check_table_exists, store_or_update
from   el_utils.ensembl import  get_species, get_compara_name, species2taxid
from   el_utils.ncbi    import  get_ncbi_tax_name


#########################################
def main():

            
    db      = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()


    #######################################################
    # check if the config db exists -- if not, make it
    qry  = "show databases like '%s'" %  "%_74%"
    rows = search_db (cursor, qry, verbose=True)
    if rows:
        for row in rows:
            db_name = row[0]
            print db_name
            qry = "drop database " + db_name
            rows2 = search_db (cursor, qry, verbose=True)


    cursor.close()
    db.close()



#########################################
if __name__ == '__main__':
    main()
