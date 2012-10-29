from mysql  import search_db



########
def taxid2name (cursor, tax_id):
    qry = "select name_txt from names where tax_id= %d " % tax_id
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""
    return rows[0][0]

########
def taxid2parentid (cursor, tax_id):
    qry = "select parent_tax_id from nodes where tax_id= %d " % tax_id
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""
    return rows[0][0]

########
def get_ncbi_tax_name (cursor):

    qry = "show databases like 'ncbi%tax%'"
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return rows[0][0]
