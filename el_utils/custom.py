
import os
from el_utils.utils   import  erropen
from el_utils.mysql   import  connect_to_mysql
from el_utils.mysql   import  switch_to_db,  search_db, store_or_update

#########################################
def get_theme_ids(cursor, cfg, theme_name):
    resources = cfg.dir_path['resources']
    fnm = resources + '/' + theme_name+'.txt'
    if not os.path.exists(fnm):
        print fnm, "not found"
        exit(1)

    if not os.path.getsize(fnm) > 0:
        print fnm, "empty"
        exit(1)
        
    inf = erropen(fnm, "r")
    gene_ids = []
    for line in inf:
        line.rstrip()
        [stable_id, name] = line.split("\t")
        qry = "select gene_id, description from gene where stable_id='%s'"% stable_id
        rows = search_db (cursor, qry)
        if not rows: continue
        gene_ids.append(rows[0][0])
    inf.close()

    return gene_ids
