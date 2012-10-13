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
    exists = rows and (type(rows[0][0]) is long) and (rows[0][0]==1)
  
 
    if exists: # if it exists, update
        qry = "update %s set" % table
        first = True
        for field, value in update_fields.iteritems():
            if (not first):
                qry += ", "
            qry += " %s = " % field
            if type(value) is int:
                qry += " %d" % value
            else:
                qry += " \'%s\'" % value

            first = False
        qry += " where %s " % conditions
        rows   = search_db (cursor, qry)
        if (rows):
            rows   = search_db (cursor, qry, verbose=True)
            return False

    else: # if not, make a new one

        qry = "insert into %s " % table
        qry += "("
        first = True
        for field in fixed_fields.keys()+update_fields.keys(): # again will have to check for the type here
            if (not first):
                qry += ", "
            qry += field
            first = False
        qry += ")"
       
        qry += " values "
        qry += "("
        first = True
        for value in fixed_fields.values()+update_fields.values(): # again will have to check for the type here
            if (not first):
                qry += ", "
            if type(value) is int:
                qry += " %d" % value
            else:
                qry += " \'%s\'" % value
            first = False
        qry += ")"
      
        rows   = search_db (cursor, qry)
        if (rows):
            rows   = search_db (cursor, qry, verbose=True)
            return False

    return True


#########################################
def create_index (cursor, db_name, index_name, table, columns):

    # columns is a list of columns that we want to have indexed
    qry = "use %s" % db_name
    rows = search_db (cursor, qry, verbose=False)
    if (rows):
        return False
    
    # check whether this index exists already
    qry = "show index from %s where key_name like '%s'" % ( table, index_name) 
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        print rows
        return True
   
    qry = "create index %s  on %s " % (index_name, table)
    qry += " ("
    first = True
    for column in columns:
        if (not first):
            qry += ", "
        qry += column
        first = False
    qry += ")"


    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        print rows
        return False
   
    return True




#########################################
def check_table_exists (cursor, db_name, table_name):
    
    qry = "use %s" % db_name
    rows = search_db (cursor, qry, verbose=False)
    if (rows):
        return False

    qry = "show tables like '%s'" % table_name
    rows = search_db (cursor, qry, verbose=False)
    if (rows):
        if ( 'Error' in rows[0]):
            return False
        else:
            return True
    else: 
        return False




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

