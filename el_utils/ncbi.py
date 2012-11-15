from mysql   import search_db, switch_to_db
from ensembl import species2taxid, get_compara_name


########
def taxid2name (cursor, tax_id):
    switch_to_db (cursor, get_ncbi_tax_name (cursor))
    qry  = "select name_txt from names where tax_id= %d " % tax_id
    qry += " and name_class = 'scientific name'";
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""
    return rows[0][0]

########
def taxid2trivial (cursor, tax_id):
    switch_to_db (cursor, get_ncbi_tax_name (cursor))
    qry  = "select name_txt from names where tax_id= %d " % tax_id
    qry += " and name_class = 'trivial'";
    rows = search_db (cursor, qry)
    if (not rows or 'ERROR' in rows[0]):
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

########
def get_common_name (cursor, species):
    switch_to_db(cursor, get_compara_name(cursor))
    tax_id = species2taxid (cursor, species)
    switch_to_db(cursor,get_ncbi_tax_name (cursor))
    qry   = "select name_txt from names where "
    qry  += "tax_id = %d and " % tax_id
    qry  += "name_class = 'genbank common name'"
    rows = search_db (cursor, qry)
    if rows:
        if ('ERROR' in rows[0]):
            search_db (cursor, qry, verbose = True)
            return ""
        else:
            return rows[0][0]

    return ""
