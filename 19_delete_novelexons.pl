#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl import  get_species, get_gene_ids


#########################################
def main():
    
    special    = 'one'
    no_threads = 1
    local_db = False
    

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    print "using", special, "set"
    gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special)

    for human_gene_id in gene_list:
        
 	human_stable      = gene2stable    (cursor, human_gene_id)
        human_description = get_description(cursor, human_gene_id)
	if verbose:  print human_gene_id, human_stable, human_description
   
  	human_exons = [e for e in gene2exon_list(cursor, human_gene_id) 
                       if e.covering_exon < 0 and e.is_canonical and e.is_known]
        if not human_exons: 
            print "\t\t no exons found"
            continue

	human_exons.sort(key=lambda exon: exon.start_in_gene)
        for he in human_exons:
            he.stable_id = exon2stable (cursor, he.exon_id)
      
        


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
