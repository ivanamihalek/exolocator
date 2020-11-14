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

		print(f"\t setting exon_seq_id to null in gene2exon")
		qry = "update gene2exon set exon_seq_id=null"
		error_intolerant_search(cursor, qry)

		print(f"\t updating exon_seq_id")
		qry = "update gene2exon g, exon_seq e  set g.exon_seq_id=e.exon_seq_id where g.exon_id=e.exon_id"
		error_intolerant_search(cursor, qry)

	cursor.close()
	db.close()
	return True

#########################################
if __name__ == '__main__':
	main()

