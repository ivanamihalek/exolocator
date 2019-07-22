#!/usr/bin/python3

from el_utils.mysql  import connect_to_mysql, search_db
from config import Config

#########################################
def main():

    db      = connect_to_mysql(Config.mysql_conf_file)
    cursor = db.cursor()

    #######################################################
    # check if the config db exists -- if not, make it
    qry  = "show databases like '%s'" %  "%_94_%"
    rows = search_db (cursor, qry, verbose=True)
    if rows:
        for row in rows:
            db_name = row[0]
            print (db_name)
            qry = "drop database " + db_name
            rows2 = search_db (cursor, qry, verbose=True)


    cursor.close()
    db.close()



#########################################
if __name__ == '__main__':
    main()
