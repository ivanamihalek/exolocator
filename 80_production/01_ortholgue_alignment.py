#!/usr/bin/python3 -u
from el_utils.mysql import  *
from el_utils.ensembl import *
from config import Config


#########################################
def main():

	if len(sys.argv)<2:
		print("usage: %s <gene symbol> " % sys.argv[0])
		exit()
	gene_name = sys.argv[1]
	species = 'homo_sapiens'  # the orthologue table is filled only here, for the moment

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	qry = "select ensembl_gene_id  from identifier_maps.hgnc where approved_symbol='%s'" % gene_name
	ensembl_id = hard_landing_search(cursor, qry)[0][0]

	[all_species, ensembl_db_name] = get_species(cursor)
	switch_to_db(cursor, ensembl_db_name[species])

	qry = "select gene_id from gene where stable_id='%s'" % ensembl_id
	gene_id = hard_landing_search(cursor, qry)[0][0]
	print(gene_name, ensembl_id, gene_id)

	for line in error_intolerant_search(cursor, "select * from orthologue where gene_id=%d" % 435660):
		print(line)

	cursor.close()
	db.close()

	return True


#########################################
if __name__ == '__main__':
	main()
