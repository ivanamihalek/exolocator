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

	qry = "select  cognate_gene_id, cognate_genome_db_id from orthologue where gene_id=%d" % gene_id
	for line in error_intolerant_search(cursor, qry):
		[cognate_gene_id, cognate_genome_db_id] = line
		qry = f"select db_name from exolocator_meta.db_names where genome_db_id={cognate_genome_db_id}"
		db_name = hard_landing_search(cursor, qry)[0][0]
		stable_canon_id = gene2stable_canon_transl_id(cursor, cognate_gene_id, db_name)
		print(db_name, cognate_gene_id, stable_canon_id)


	cursor.close()
	db.close()

	return True


#########################################
if __name__ == '__main__':
	main()
