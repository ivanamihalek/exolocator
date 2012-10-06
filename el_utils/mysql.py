import MySQLdb


########
def store_or_update (cursor, table, fixed_fields, update_fields):


    conditions = ""
    first = True
    for [field, value] in fixed_fields.iteritems():
        if (not first):
            conditions += " and "
        if ( type (value) is int):
            conditions += " %s=%d "   % (field, value)
        else:
            conditions += " %s='%s' " % (field, value)
        first = False

    # check if the row exists
    qry = "select exists (select 1 from %s  where %s) "  % (table, conditions)
    rows   = search_db (cursor, qry)
    exists = rows and type(rows[0][0]) is int and rows[0][0]==1
   
    if exists: # if it exists, update
        qry = "update"

    else: # if not, make a new one

        qry = "insert into %s " % table
        qry += "("
        for fields in fixed_fields: # again will have to check for the type here



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

