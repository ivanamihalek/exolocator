


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

	model_orgs = ['homo_sapiens', 'mus_musculus', 'danio_rerio', 'xenopus_tropicalis', 'gallus_gallus']
	species_subtrees = ['Archosauria',  # birds and crocs
	                    'Testudines',  # turtles
	                    'Lepidosauria',  # scaled reptiles (lizards and snakes)
	                    'Eutheria', 'Marsupialia',
	                    'ornithorhynchus_anatinus',  # platipus
	                    'Anura',  # frogs
	                    'latimeria_chalumnae',  # coelacanth
	                    'Euteleosteomorpha', 'Otomorpha', 'Osteoglossiformes',  # ray-finned fish
	                    'lepisosteus_oculatus', 'erpetoichthys_calabaricus',  # more ray-finned fish
	                    'callorhinchus milii',  # Australian ghostshark or elephant shark
	                    'Cyclostomata'  # nightmare stuff
	                    ]
