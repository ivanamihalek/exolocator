import MySQLdb


#######
def search_db (cursor, qry, verbose=False):


    try:
        cursor.execute(qry)
    except MySQLdb.Error, e:
        if verbose:
            print "Error running cursor.execute() for  qry: %s: %s " % (qry, e.args[1])
        return  ["ERROR: "+e.args[1]]


    try:
        rows = cursor.fetchall()
    except MySQLdb.Error, e:
        if verbose:
            print "Error running cursor.fetchall() for  qry: %s: %s " % (qry, e.args[1])
        return  ["ERROR: "+e.args[1]]
    
    if (len(rows) == 0):
        if verbose:
            print "No return for query %s"  % qry
        return False

    return rows

########
def connect_to_mysql ():
    try:
        db =MySQLdb.connect(user="root")
    except  MySQLdb.Error, e:
        print "Error connecting to mysql as root (%s) " % (e.args[1])
        exit(1)
 
    return db

########
def connect_to_db (db_name):

    try:
        db =MySQLdb.connect(user="root", db=db_name)
    except  MySQLdb.Error, e:
        print "Error connecting to %s: %d %s" % (db_name, e.args[0], e.args[1])
        exit(1)

    return db

