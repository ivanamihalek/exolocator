#!/usr/bin/python3

from el_utils.mysql  import *
from config import Config


def drop(cursor, db_name):
	print (f"dropping {db_name}")
	error_intolerant_search(cursor, f"drop database {db_name}")


#########################################
def main():
	version = 97

	db      = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#######################################################
	qry  = f"show databases like '_{version}_'"
	rows = error_intolerant_search(cursor, qry)
	if rows:
		for row in rows:
			db_name = row[0]
			drop(cursor, db_name)

	db_name = f"ensembl_compara_{version}"
	drop(cursor, db_name)

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
