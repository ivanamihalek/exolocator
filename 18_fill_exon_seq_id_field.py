#!/usr/bin/python3 -u

from config import Config
from el_utils.el_specific import  *

#########################################
def main():

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species (cursor)
	for species in all_species:
		print(species)
		print(f"\t making index")
		create_index(cursor,  ensembl_db_name[species],'gene2exon','exon_id_idx', ['exon_id'])
		print(f"\t index done")
		switch_to_db(cursor, ensembl_db_name[species])
		print(f"\t updating exon_seq_id")
		time0 = time()
		qry = "select exon_seq_id, exon_id from exon_seq"
		for exon_seq_id, exon_id in hard_landing_search(cursor, qry):
			qry = f"update gene2exon set exon_seq_id={exon_seq_id} where exon_id={exon_id}"
			error_intolerant_search(cursor, qry)
		print("\t done in %.1fmin"%((time()-time0)/60))
	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()

