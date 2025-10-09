#!/usr/bin/python3 -u

from el_utils.processes import *
from el_utils.special_gene_sets import *
from el_utils.known_exon_utils import exons_for_gene

from config import Config



#########################################
def exons_for_species(species_list, db_info):

	[ensembl_db_name, outdir] = db_info
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor, "set autocommit=1")
	logf = open("log.{}.txt".format(get_process_id()),"w")

	for species in species_list:
		# load this file later with
		# sudo mysqlimport  --local <species db>  outfir/species/gene2exon.tsv
		# mysqlimport strips any extension and uses what's left as a table name
		# before you begin, do
		# mysql> SET GLOBAL local_infile = 1;
		os.makedirs(f"{outdir}/{species}", exist_ok=True)
		outfile = open(f"{outdir}/{species}/gene2exon.tsv", "w")

		logf.write(species+" started\n")
		gene_ids = get_gene_ids(cursor, biotype='protein_coding', db_name= ensembl_db_name[species])
		# gene_ids = [8979]
		count = 0
		time0 = time()
		for gene_id in gene_ids:
			if gene_ids.index(gene_id)%1000==0:
				pct_of_genes_processed = float(int(gene_ids.index(gene_id)) + 1)/len(gene_ids)*100
				print("%50s:  %5.1f%%    %ds" % (species, pct_of_genes_processed, time()-time0))
				time0 = time()
			count = exons_for_gene(cursor, gene_id, ensembl_db_name[species], count, outfile, logf)
		outfile.close()
		print(f"{species} done", file=logf)
	logf.close()
	cursor.close()
	db.close()

	return True


def check_species_done( all_species,  outdir):
	unprocessed_species = []
	for species in all_species:
		gene2ex_file = f"{outdir}/{species}/gene2exon.tsv"
		if not os.path.exists(gene2ex_file):
			unprocessed_species.append(species)

	return unprocessed_species

#########################################
def main():
	outdir = "/home/ivana/scratch/raw_tables"
	os.makedirs(outdir, exist_ok=True)

	no_threads = 32
	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	#all_species = ["mus_musculus"]
	#all_species.remove('homo_sapiens')

	cursor.close()
	db    .close()

	unprocessed_species = check_species_done(all_species, "raw_tables")

	parallelize(no_threads, exons_for_species, unprocessed_species, [ensembl_db_name, outdir])


#########################################
if __name__ == '__main__':
	main()

'''
In v 101 a bunch of canonical transcript coordinates were missing for mus caroli
not sure if I should worry about that - there ar 75 cases like that
exmaple MGP_CAROLIEiJ_G0027698 Ppp1cc protein phosphatase 1 catalytic subunit gamma
The smae for mus pahari and mus spretus
'''