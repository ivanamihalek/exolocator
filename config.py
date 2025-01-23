
class Config:

	release_number = "110"

	fasta_repo = f"/storage/databases/ensembl-{release_number}/fasta"

	mysql_conf_file = "/home/ivana/.tcga_conf"
	ucsc_mysql_conf_file  = "/home/ivana/.ucsc_mysql_conf"
	blastp = "/usr/bin/blastp"
	blastdbcmd = "/usr/bin/blastdbcmd"
	muscle = "/usr/third/muscle"
	seqtk = "/usr/third/seqtk/seqtk"
	# exon flankingin region - not clear how long is the splice length
	exon_flanking_region_length = 15
	min_accptbl_exon_sim = 0.7

	model_orgs = ['homo_sapiens', 'mus_musculus', 'danio_rerio', 'xenopus_tropicalis', 'gallus_gallus']
