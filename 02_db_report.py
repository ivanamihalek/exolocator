#!/usr/bin/python


import MySQLdb
from   el_utils.mysql import connect_to_mysql


#########################################
def main():
    
    db = connect_to_mysql()
    cursor = db.cursor()

    qry = "show databases like '%core%'"
    cursor.execute(qry)

    rows = cursor.fetchall()
    if (not rows):
        print "No databases with 'core' in the name found"
        return 1

    ensembl_db_name = {}
    all_species = []
    for row in rows:
        db_name = row[0]
        name_tokens = db_name.split ('_')
        species = name_tokens[0]+'_'+ name_tokens[1]
        ensembl_db_name[species] = db_name
        all_species.append(species)

    for species in all_species:
        print species

    print "there are %d core dbs available " % len(all_species)

#########################################
if __name__ == '__main__':
    main()
