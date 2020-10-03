


class Config:

	release_number = "101"

	fasta_repo = f"/storage/databases/ensembl-{release_number}/fasta"

	mysql_conf_file = "/home/ivana/.tcga_conf"
	ucsc_mysql_conf_file  = "/home/ivana/.ucsc_mysql_conf"
	blastp="/usr/bin/blastp"
	blastdbcmd = "/usr/bin/blastdbcmd"
	muscle ="/usr/bin/muscle"

	# exon flankingin region - not clear how long is the splice length
	exon_flanking_region_length = 15
