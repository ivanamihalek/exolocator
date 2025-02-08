
class Config:

	release_number = "113"
	fasta_repo = f"/media/ivana/portable/ensembl-{release_number}/fasta"
	mysql_repo = f"/media/ivana/portable/ensembl-{release_number}/mysql"
	# Skip species lists
	skip_species = [
		"ancestral_alleles", "bos_indicus_x_bos_taurus",
		"canis_lupus_familiarisbasenji", "canis_lupus_familiarisboxer",
		"canis_lupus_familiarisgreatdane", "canis_lupus_familiarisgsd",
		"cyprinus_carpio_germanmirror", "cyprinus_carpio_hebaored", "cyprinus_carpio_huanghe",
		"caenorhabditis_elegans", "ciona_intestinalis",
		"ciona_savignyi", "drosophila_melanogaster", "saccharomyces_cerevisiae"
	]
	# breeds to skip (when we have the corresponding generic genome)
	breed_species = [
		"anas_platyrhynchos_", "astyanax_mexicanus_", "bos_taurus_", "capra_hircus_",  "gallus_gallus_",
		"gasterosteus_aculeatus_", "hybrid", "mus_musculus_", "ovis_aries_",
		"oryzias_latipes_", "rattus_norvegicus_", "salmo_salar_", "sus_scrofa_",
	]

	ucsc_mysql_conf_file  = "/home/ivana/.ucsc_mysql_conf"
	blastp = "/usr/bin/blastp"
	blastdbcmd = "/usr/bin/blastdbcmd"
	muscle = "/usr/third/muscle"
	seqtk = "/usr/third/seqtk/seqtk"
	# exon flankingin region - not clear how long is the splice length
	exon_flanking_region_length = 15
	min_accptbl_exon_sim = 0.7

	model_orgs = ['homo_sapiens', 'mus_musculus', 'danio_rerio', 'xenopus_tropicalis', 'gallus_gallus']

	scratch = "/home/ivana/scratch"

# note from release 113:
#  ovis_aries dir contained ovis_aries_rambouillet_core_113_2.sql
# which I renamed to ovis_aries_core_113_2.sql, adn change the name of the db therein as well
# it was the only place where the naming was inconsistent in this way
