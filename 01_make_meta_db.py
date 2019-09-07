#!/usr/bin/python3

""" 
@package 01_make_meta_db.py
@brief  Creating ad filling the configuration database. Hardcoded stuff is here.
@author Ivana
@date   2012, 2019
@mainpage Overview
@main @section Purpose
@main  exolocator pipeline is intended to fill exolocator_db, the database behind the ExoLocator server. 
@main  ExoLocator server, in turn, is aimed at being a one-stop shopping point for comparative 
@main analysis of protein sequences.
@main @section  Design
@main  The main aim in desinging the pipeline has been intuitive use and reproducibility:  
@main by using the scripts in the order in which they are enumerated, the user should be able to
@main  recreate and/or update exolocator_db.
@main  @section Databases
@main exolocator pipeline functions by manipulating databases: everything except full genomic sequences is 
@main is stored in a database, from the configuration of the pipeline itself, to core anntotation data for each species
@main from Ensembl, and taxonomy information info from NCBI.
"""

from el_utils.mysql   import *
from el_utils.ensembl import get_species, get_compara_name, species2taxid
from el_utils.ncbi    import get_ncbi_tax_name
from config import Config

def set_trivial_names():
	trivial = {}
	trivial['ailuropoda_melanoleuca'] = 'panda'
	trivial['anas_platyrhynchos_platyrhynchos']     = 'duck'
	trivial['anolis_carolinensis']    = 'anole_lizard'
	trivial['astyanax_mexicanus']     = 'blind_cavefish'
	trivial['bison_bison_bison']      = 'bison'
	trivial['bos_taurus']             = 'cow'
	trivial['callithrix_jacchus']     = 'marmoset'
	trivial['canis_familiaris']       = 'dog'
	trivial['carlito_syrichta']       = 'tarsier'
	trivial['cebus_capucinus']        = 'capuchin_monkey'
	trivial['cavia_porcellus']        = 'guinea_pig'
	trivial['choloepus_hoffmanni']    = 'sloth'
	trivial['danio_rerio']            = 'zebrafish'
	trivial['dasypus_novemcinctus']   = 'armadillo'
	trivial['dipodomys_ordii']        = 'kangaroo_rat'
	trivial['echinops_telfairi']      = 'madagascar_hedgehog'
	trivial['equus_asinus_asinus']    = 'donkey'
	trivial['equus_caballus']         = 'horse'
	trivial['erinaceus_europaeus']    = 'european_hedgehog'
	trivial['felis_catus']            = 'cat'
	trivial['ficedula_albicollis']    = 'flycatcher'
	trivial['gadus_morhua']           = 'cod'
	trivial['gallus_gallus']          = 'chicken'
	trivial['gasterosteus_aculeatus'] = 'stickleback'
	trivial['gorilla_gorilla']        = 'gorilla'
	trivial['heterocephalus_glaber_female']  = 'naked_mole_rat_female'
	trivial['heterocephalus_glaber_male']  = 'naked_mole_rat_male'
	trivial['homo_sapiens']                = 'human'
	trivial['ictidomys_tridecemlineatus']  = 'squirrel'
	trivial['latimeria_chalumnae']         = 'coelacanth'
	trivial['lepisosteus_oculatus']        = 'spotted_gar'
	trivial['loxodonta_africana']          = 'elephant'
	trivial['macaca_mulatta']              = 'macaque'
	trivial['notamacropus_eugenii']        = 'wallaby'
	trivial['meleagris_gallopavo']         = 'turkey'
	trivial['microcebus_murinus']          = 'lemur'
	trivial['monodelphis_domestica']       = 'opossum'
	trivial['mus_musculus']                = 'mouse'
	trivial['mustela_putorius_furo']       = 'ferret'
	trivial['myotis_lucifugus']            = 'bat'
	trivial['nomascus_leucogenys']         = 'gibbon'
	trivial['ochotona_princeps']           = 'pika'
	trivial['oreochromis_niloticus']       = 'tilapia'
	trivial['ornithorhynchus_anatinus']    = 'platypus'
	trivial['oryctolagus_cuniculus']       = 'rabbit'
	trivial['oryzias_latipes']             = 'medaka'
	trivial['otolemur_garnettii']          = 'galago_lemur'
	trivial['ovis_aries']                  = 'sheep'
	trivial['pan_troglodytes']             = 'chimpanzee'
	trivial['papio_anubis']                = 'baboon'
	trivial['pelodiscus_sinensis']         = 'turtle'
	trivial['petromyzon_marinus']          = 'lamprey'
	trivial['poecilia_formosa']            = 'amazon_molly'
	trivial['pongo_abelii']                = 'orangutan'
	trivial['procavia_capensis']           = 'hyrax'
	trivial['pteropus_vampyrus']           = 'flying_fox'
	trivial['pygocentrus_nattereri']       = 'piranha'
	trivial['rattus_norvegicus']           = 'rat'
	trivial['sarcophilus_harrisii']        = 'tasmanian_devil'
	trivial['sorex_araneus']               = 'european_shrew'
	trivial['sus_scrofa']                  = 'pig'
	trivial['taeniopygia_guttata']         = 'zebra_finch'
	trivial['takifugu_rubripes']           = 'fugu'
	trivial['tetraodon_nigroviridis']      = 'pufferfish'
	trivial['tupaia_belangeri']            = 'tree_shrew'
	trivial['tursiops_truncatus']          = 'dolphin'
	trivial['vicugna_pacos']               = 'alpaca'
	trivial['xenopus_tropicalis']          = 'xenopus'
	trivial['xiphophorus_maculatus']       = 'platyfish'

	return trivial

#########################################
def  make_parameter_table (cursor, table):

	"""
	Creates parameter table in the config database.
	@param [cursor] db cursor, assumed top be pointing to the config database
	@retval True  on success
	@retval False on failure;  in that case the seach_db() call is repeated in verbose mode.
	"""

	print("making ", table)

	qry  = "create table " + table + "  (id int(10) primary key auto_increment)"
	rows = search_db (cursor, qry, verbose=True)
	if (rows):
		return False

	# make the columns
	for column  in  ['name', 'value']:
		qry = "alter table  %s  add  %s  varchar (50)" % (table, column)
		rows = search_db (cursor, qry, verbose=True)
		if (rows):
			return False
	return False

#########################################
def  make_flags_table (cursor, db_name, table):
	switch_to_db (cursor, db_name)
	table = 'flags'
	if check_table_exists(cursor, db_name, table):
		qry = "drop table " + table
		search_db(cursor, qry, verbose=True)

	qry = ""
	qry += "CREATE TABLE  %s (" % table
	qry += "     id INT NOT NULL, "
	qry += "  	 genome_db text NOT NULL, "
	qry += "  	 flag VARCHAR(30) NOT NULL, "
	qry += "  	 raised_by VARCHAR(50) NOT NULL, "
	qry += "  	 comment text NOT NULL, "
	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)




#########################################
def make_seqregion2file_table (cursor, db_name):

	switch_to_db (cursor, db_name)
	table_name = 'seqregion2file'
	if check_table_exists(cursor, db_name, table_name):
		qry = "drop table " + table_name
		search_db(cursor, qry, verbose=True)

	qry = ""
	qry += "  CREATE TABLE  %s (" % table_name
	qry += "     seqregion_id INT(10)  NOT NULL, "
	qry += "  	 seq_name VARCHAR (100) NOT NULL, "
	qry += "  	 file_name TEXT NOT NULL, "

	qry += "	 PRIMARY KEY (seqregion_id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)
	return


#########################################
def make_species_names_table(cursor, db_name):

	switch_to_db (cursor, db_name)
	table_name = 'species_names'
	if check_table_exists(cursor, db_name, table_name):
		qry = "drop table " + table_name
		search_db(cursor, qry, verbose=True)

	qry = ""
	qry += "  CREATE TABLE  %s (" % table_name
	qry += "     id INT(10)  NOT NULL AUTO_INCREMENT, "
	qry += "     tax_id INT(10)  NOT NULL, "
	qry += "  	 species VARCHAR (100) NOT NULL, "
	qry += "  	 shorthand VARCHAR (20) NOT NULL, "
	qry += "  	 trivial_name VARCHAR (50) NOT NULL, "

	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)
	return


#########################################
def feed_parameters(cursor, db_name):

	parameter = {}
	# in case I ever have to handle multiple versions of ensembl
	# (but for now I don't have enough space)
	# note though that there are functions in el_utils/mysql.py that assume
	# that whatever ensembl stuff is available to the mysql server corresponds to the same release
	parameter['ensembl_release_number'] = '97'
	parameter['blastp_e_value']         = "1.e-10" # it will be used as a string  when fmting the blastp cmd
	parameter['min_accptbl_exon_sim']   = 0.33333 #minimum acceptable exon similarity

	table = 'parameters'
	if check_table_exists(cursor, db_name, table):
		print(table, " found in ", db_name)
	else:
		print(table, " not found in ", db_name)
		make_parameter_table(cursor, table)

	fixed_fields  = {}
	update_fields = {}
	for [name, value] in parameter.items():
		fixed_fields['name']  = name
		update_fields['value'] = value
		store_or_update (cursor, table, fixed_fields, update_fields)


#########################################
def get_trivial_names(cursor, all_species):

	tax_id  = {}
	db_name = get_compara_name (cursor)
	if (not db_name):
		print("compara db not found")
		exit(1)
	switch_to_db(cursor, db_name)
	for species in all_species:
		tax_id[species] = species2taxid(cursor, species)

	# switch to ncbi taxonomy database
	tax_db_name = get_ncbi_tax_name (cursor)
	if not tax_db_name:
		print("ncbi taxonomy db not found")
		exit(1)
	switch_to_db(cursor, tax_db_name)

	common_name = {}
	for species in all_species:
		qry  = "select * from ncbi_taxa_name where taxon_id=%d " % tax_id[species]
		qry += "and name_class='genbank common name'"
		ret = search_db(cursor, qry)
		if not ret: continue
		common_name[species] = ret[0][1]

	trivial_names = set_trivial_names()
	for species in all_species:
		trivial = trivial_names.get(species, "not found")
		common = common_name.get(species, "not found")
		if trivial=="not found":
			if common == "not found":
				trivial = species.lower().replace(" ", "_").replace("'", "_")
			else:
				trivial = common.replace("common ","").lower().replace(" ", "_").replace("'", "_")
			trivial_names[species] = trivial

	return trivial_names, tax_id


# dingo and dog both have ENSCAF prefix in Ensembl:
# https://www.ensembl.org/Canis_lupus_dingo/Info/Annotation says:
# The genetic evidence indicates that the dingo clade originated from East Asian domestic dogs
# and was introduced through the Malay Archipelago into Australia. Based on a comparison with early fossils,
# dingo morphology has not changed over thousands of years.
# This suggests that there has been no artificial selection over this period
# and that the dingo represents an early form of dog.
#
# for chinese hamster, we have three different assmeblies
#########################################
def get_shorthand(cursor, species, ensembl_db):
	if species=='homo_sapiens':
		shorthand='HSA'
	elif species=='canis_lupus_dingo':
		shorthand='CLD'
	elif 'cricetulus_griseus' in species: # (three assemblies)
		shorthand = "CGR_" + species.split("_")[2].upper()
	else:
		qry="select meta_value  from %s.meta where meta_key='species.stable_id_prefix'" % ensembl_db
		stable_id_prefix = hard_landing_search(cursor, qry)[0][0]
		shorthand = stable_id_prefix.replace("ENS","").replace("EiJ_","").upper()
	return shorthand

#########################################
def tax_id_sanity_check(cursor, species, ensembl_db, ncbi_taxonomy_id):
	ncbi_taxonomy_id = int(ncbi_taxonomy_id)
	qry = "select meta_value from %s.meta where meta_key='species.taxonomy_id'" % ensembl_db
	tax_id_ensembl =  int(hard_landing_search(cursor, qry)[0][0])
	if tax_id_ensembl != ncbi_taxonomy_id:
		print("tax_id mismatch for " + species)
		print(tax_id_ensembl, ncbi_taxonomy_id)
		exit()
	return ncbi_taxonomy_id


#########################################
def feed_species_names (cursor, exolocator_meta_db_name, all_species, ensembl_db_name, trivial_name, ncbi_tax_id):

	make_species_names_table(cursor, exolocator_meta_db_name)

	dbs_with_shorthand = {}
	for species in all_species:
		shorthand = get_shorthand(cursor, species, ensembl_db_name[species])
		if not shorthand in dbs_with_shorthand: dbs_with_shorthand[shorthand]=[]
		dbs_with_shorthand[shorthand].append(ensembl_db_name[species])
		tax_id = tax_id_sanity_check(cursor, species, ensembl_db_name[species], ncbi_tax_id[species])
		print("%8d %35s %45s %20s   %s (%d)" % (tax_id, species, ensembl_db_name[species],
												shorthand, trivial_name[species], len(trivial_name[species])))
		fixed_fields  = {}
		update_fields = {}
		fixed_fields  ['species']   = species
		update_fields ['tax_id']    = tax_id
		update_fields ['shorthand'] = shorthand
		update_fields ['trivial_name'] = trivial_name[species]
		switch_to_db(cursor, exolocator_meta_db_name)
		store_or_update (cursor, 'species_names', fixed_fields, update_fields)

	# sanity warning - should not happen - redo with fiexd name shorthands if it does
	for shorthand, dbs in dbs_with_shorthand.items():
		if len(dbs)==1: continue
		print("\n dbs with same shorthand:")
		print(shorthand, dbs)
		print("This should not have happened. Redo with fiexd name shorthands.")
		exit()


#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	search_db(cursor,"set autocommit=1")


	#######################################################
	# check if the config db exists -- if not, make it
	exolocator_meta_db_name   = "exolocator_meta"
	qry  = "show databases like'%s'" % exolocator_meta_db_name
	rows = search_db (cursor, qry)
	if (not rows):
		print(exolocator_meta_db_name, "database not found")
		qry = "create database %s " % exolocator_meta_db_name
		rows = search_db (cursor, qry)
		if (rows):
			print("some problem creating the database ...")
			search_db (cursor, qry, verbose=True)
		print(exolocator_meta_db_name, "database found")

	switch_to_db(cursor, exolocator_meta_db_name)

	#######################################################
	# create flags database (for flagging arbitrary problems - to be filled as we go)
	make_flags_table(cursor, exolocator_meta_db_name, 'flags')

	#######################################################
	# store parameters, to make sure we are consistent and reproducible
	feed_parameters(cursor, exolocator_meta_db_name)

	#######################################################
	[all_species, ensembl_db_name] = get_species(cursor)

	#######################################################
	# get trivial names
	trivial_names, ncbi_tax_id = get_trivial_names (cursor, all_species)

	#######################################################
	# add species shorthands (used in ENS* names formation)
	# though we will not needed it until the paralogue alignment reconstruction point)
	feed_species_names (cursor, exolocator_meta_db_name, all_species, ensembl_db_name, trivial_names, ncbi_tax_id)

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
