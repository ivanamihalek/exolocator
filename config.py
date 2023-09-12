
class Config:

	release_number = "101"

	fasta_repo = f"/storage/databases/ensembl-{release_number}/fasta"

	mysql_conf_file = "/home/ivana/.tcga_conf"
	ucsc_mysql_conf_file  = "/home/ivana/.ucsc_mysql_conf"
	blastp = "/home/ivana/third/blast/bin/blastp"
	blastdbcmd = "/home/ivana/third/blast/bin/blastdbcmd"
	muscle = "/home/ivana/third/muscle"

	# exon flankingin region - not clear how long is the splice length
	exon_flanking_region_length = 15
	min_accptbl_exon_sim = 0.7

	model_orgs = ['homo_sapiens', 'mus_musculus', 'danio_rerio', 'xenopus_tropicalis', 'gallus_gallus']
