import os
from wsgiref.simple_server import server_version

import MySQLdb, sys, warnings
from time import time

from dotenv import load_dotenv

load_dotenv()

#########################################
def error_intolerant_search(cursor, qry):
	ret =  search_db(cursor, qry)
	if not ret: return ret
	if type(ret[0][0])==str and 'error' in ret[0][0].lower():
		search_db(cursor, qry, verbose=True)
		exit()
	return ret


#########################################
def hard_landing_search(cursor, qry):
	ret = search_db(cursor, qry)
	if not ret or (type(ret[0][0])==str and 'error' in ret[0][0].lower()):
		search_db(cursor, qry, verbose=True)
		exit()
	return ret


########
def check_null(variable):
	if variable is None:
		return None
	if type(variable) is str and variable == "None":
		return None
	return variable


########
def switch_to_db(cursor, db_name):
	qry = "use %s" % db_name
	rows = search_db(cursor, qry, verbose=False)
	if (rows):
		print(rows)
		return False
	return True


########################################
def val2mysqlval(value):
	if value is None:
		return "null"
	elif type(value) is str:
		return "\'%s\'" % value.replace("'", "\\'")
	return f"{value}"


#########################################
def store_without_checking(cursor, table, fields, verbose=False, database=None, ignore=False):

	qry = "insert "
	# if we use INSERT IGNORE, the duplication attempt is ignored
	if ignore: qry += "ignore "
	qry += "%s.%s "%(database,table) if database else "%s "%table

	qry += "("
	qry += ",".join(fields.keys())
	qry += ")"

	qry += " values "
	qry += "("
	qry += ",".join([val2mysqlval(v) for v in fields.values()])
	qry += ")"

	rows  = search_db (cursor, qry, verbose)
	if verbose: print("qry:",qry,"\n", "rows:", rows)

	if rows:
		rows   = search_db (cursor, qry, verbose=True)
		print(rows)
		return -1

	rows = search_db (cursor, "select last_insert_id()" )
	try:
		row_id = int(rows[0][0])
	except:
		row_id = -1
	return row_id


#########################################
def store_or_update(cursor, table, fixed_fields, update_fields=None, verbose=False, primary_key='id'):
	# the update heuristics below expect that the type us either in or str

	primary_keys = []

	# check if the row exists
	conditions = " and ".join([f"{k}={val2mysqlval(v)}" for k, v in fixed_fields.items()])
	qry = "select %s from %s  where %s "  % (primary_key, table, conditions)
	rows   = error_intolerant_search(cursor, qry)
	exists = rows and len(rows)>0

	if exists:
		primary_keys = [row[0] for row in rows]
		row_id_string = ",".join([f"{row[0]}"for row in rows])

		if verbose: print(f"\nqry + exists? {exists} {row_id_string}")
		if not update_fields: return primary_keys

		if verbose: print("exists; updating")
		qry  = "update %s set " % table
		qry += ",".join([f"{k}={val2mysqlval(v)}" for k, v in update_fields.items()])
		qry += f" where {primary_key} in ({row_id_string})"

	else:  # if not, make a new one
		if verbose: print("does not exist; making new one")
		qry  = "insert into %s " % table
		keys = list(fixed_fields.keys())
		vals = list(fixed_fields.values())
		if update_fields:
			keys += list(update_fields.keys())
			vals += list(update_fields.values())
		qry += "(" + ",".join(keys) + ")"
		qry += " values "
		qry += "(" + ",".join([val2mysqlval(v) for v in vals]) + ")"

	rows   = search_db(cursor, qry, verbose)

	if verbose: print("qry:", qry, "\n", "rows:", rows)
	# if there is a return, it is an error msg
	if rows:
		rows   = search_db(cursor, qry, verbose=True)
		print(rows[0])
		return False

	if not exists:
		rows = search_db (cursor, "select last_insert_id()" )
		try:
			primary_keys = rows[0]
		except:
			return False

	return primary_keys


#########################################
def create_index(cursor, db_name, table, index_name, columns, verbose=False):
	time0 = 0
	# check whether this index exists already
	if verbose: print("checking existence of index {} on table {} in {}".format(index_name, table, db_name))
	qry = "show index from %s.%s where key_name like '%s'" % (db_name, table, index_name)
	rows = search_db(cursor, qry, verbose=False)
	if rows:
		if verbose: print("found index {} on table {} in {}.".format(index_name, table, db_name))
		return True

	if verbose:
		time0 = time()
		print("creating index {} on table {} in {} ...".format(index_name, table, db_name))
	# columns is a list of columns that we want to have indexed
	qry = "create index %s on %s.%s (%s)" % (index_name, db_name, table, ",".join(columns))
	error_intolerant_search(cursor, qry)
	if verbose:
		time_min = (time() - time0)/60
		print("\t\t  done in %.1f min." %  time_min)

	return True


#########################################
def column_exists(cursor, db_name, table_name, column_name):
	if not switch_to_db(cursor, db_name):
		return False

	qry = "show columns from " + table_name + " like '%s'" % column_name
	rows = search_db(cursor, qry, verbose=False)
	if (rows):
		if ( 'Error' in rows[0]):
			return False
		else:
			return True
	else:
		return False


def count_table_rows (cursor, db_name, table_name):
	# not my problem if the table does not exist
	qry = f"select count(*) from {db_name}.{table_name}"
	rows = search_db(cursor, qry, verbose=False)
	if rows:
		if 'Error' in rows[0]:
			return 0
		else:
			return rows[0][0]
	else:
		return 0


#########################################
def add_column(cursor, db_name, table_name, column_name, col_type, default=None, after_col=None):
	if not column_exists (cursor, db_name, table_name, column_name):
		qry = "alter table  %s.%s add  %s %s  " %(db_name, table_name, column_name, col_type)
		if default: qry += "default %s " % default
		if after_col: qry += "after %s" % after_col
		error_intolerant_search(cursor,qry)
	return


#########################################
def get_column_names (cursor, db_name, table_name):

	qry = f"show columns from {db_name}.{table_name}"
	rows = hard_landing_search (cursor, qry)
	return [row[0] for row in rows]


#########################################
def check_table_exists(cursor, db_name, table_name):
	if not switch_to_db(cursor, db_name):
		return False

	qry = "show tables like '%s'" % table_name
	rows = search_db(cursor, qry, verbose=True)
	if rows:
		if 'Error' in rows[0]:
			return False
		else:
			return True
	else:
		return False


############
def check_and_drop_table(cursor, db_name, table):
	search_db(cursor, "drop table if exists %s.%s"% (db_name, table))
	return


#########################################
def table_create_time(cursor, db_name, table_name):
	qry = "select create_time from information_schema.tables where "
	qry += "table_schema   = '%s' " % db_name
	qry += "and table_name = '%s' " % table_name

	rows = search_db(cursor, qry, verbose=False)
	if (not rows or 'Error' in rows[0]):
		search_db(cursor, qry, verbose=True)
		return ""
	else:
		return rows[0][0]


#######
def search_db(cursor, qry, verbose=False):
	warnings.filterwarnings('ignore', category=MySQLdb.Warning)
	try:
		cursor.execute(qry)
	except MySQLdb.Error as e:
		if verbose:
			print("Error running cursor.execute() for  qry:\n%s\n%s" % (qry, e.args[1]))
		return [["Error"], e.args]
	except MySQLdb.Warning as e: # this does not work for me - therefore filterwarnings
		if verbose:
			print("Warning running cursor.execute() for  qry:\n%s\n%s" % (qry, e.args[1]))
		return [["Warning"], e.args]

	try:
		rows = cursor.fetchall()
	except MySQLdb.Error as e:
		if verbose:
			print("Error running cursor.fetchall() for  qry:\n%s\n%s" % (qry, e.args[1]))
		return [["Error"], e.args]

	if len(rows) == 0:
		if verbose:
			print("No return for query:\n%s" % qry)
		return False

	# since python3 fetchall returns bytes inst of str in some  random fashion
	# not clear what's going on
	# here is a rather useless issue page on github
	# https://github.com/PyMySQL/mysqlclient-python/issues/145#issuecomment-283936456
	rows_clean = []
	for row in rows:
		rows_clean.append([r.decode('utf-8') if type(r)==bytes else r for r in row])
	return rows_clean


###################
import json

def read_creds ():

	[user, passwd, host, port] = [None, None, None, None]
	with open('.creds') as data_file:
		data = json.load(data_file)

	if 'user' in data:   user   = data['user']
	if 'passwd' in data: passwd = data['passwd']
	if 'host' in data:   host   = data['host']
	if 'port' in data:   port   = data['port']

	return [user, passwd, host, port]


########
def connect_to_mysql (conf_file):
	try:
		mysql_conn_handle = MySQLdb.connect(read_default_file=conf_file)
	except  MySQLdb.Error as e:
		print(("Error connecting to mysql (%s) " % (e.args[1])))
		sys.exit(1)
	return mysql_conn_handle


########
def mysql_server_connect (user=None, passwd=None, host='localhost', port: int = 3306):
	try:
		# data cannot be input from a local file without the last argument
		db = MySQLdb.connect(user=user, passwd=passwd, host=host, port=port, local_infile=1)
	except  MySQLdb.Error as e:
		print(("Error connecting to mysql server: %d %s" % (e.args[0], e.args[1])))
		exit(1)
	cursor = db.cursor()
	return cursor

def mysql_using_env_creds():
	cursor = mysql_server_connect(user=os.getenv('MYSQL_USER'),
                                  passwd=os.getenv('MYSQL_PASSWORD'),
                                  host=os.getenv('MYSQL_HOST', 'localhost'),
                                  port=int(os.getenv('MYSQL_PORT', 3306)))
	return cursor

def db_connect(db_name, user=None, passwd=None, host='localhost'):

	cursor = mysql_server_connect(user, passwd, host)
	switch_to_db(cursor, db_name)
	search_db(cursor, "set autocommit = 1")  # ... I thought it was a default
	return cursor


def mysql_server_conn_close(cursor):
	db = cursor.connection
	cursor.close()
	db.close()

