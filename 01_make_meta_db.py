#!/usr/bin/python3

""" 
@package 00_make_meta_db.py
@brief  Creating ad filling the configuration database. Hardcoded stuff is here.
@author Ivana
@date   2012, 2025
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
import os

from dotenv import load_dotenv

from el_utils.mysql import *
from el_utils.ensembl import get_ensembl_species, get_compara_name, ncbi_species2taxid
from config import Config

# Load environment variables
load_dotenv()


def set_trivial_names():
    trivial = {}
    trivial['ailuropoda_melanoleuca'] = 'panda'
    trivial['anas_platyrhynchos_platyrhynchos'] = 'duck'
    trivial['anas_zonorhyncha'] = 'spot-billed_duck'
    trivial['anolis_carolinensis'] = 'anole_lizard'
    trivial['astyanax_mexicanus'] = 'blind_cavefish'
    trivial['bison_bison_bison'] = 'bison'
    trivial['bos_taurus'] = 'cow'
    trivial['callithrix_jacchus'] = 'marmoset'
    trivial['canis_familiaris'] = 'dog'
    trivial['cairina_moschata_domestica'] = 'muscovy_duck'
    trivial['carlito_syrichta'] = 'tarsier'
    trivial['cebus_capucinus'] = 'capuchin_monkey'
    trivial['cavia_porcellus'] = 'guinea_pig'
    trivial['choloepus_hoffmanni'] = 'sloth'
    trivial['danio_rerio'] = 'zebrafish'
    trivial['dasypus_novemcinctus'] = 'armadillo'
    trivial['dipodomys_ordii'] = 'kangaroo_rat'
    trivial['echinops_telfairi'] = 'madagascar_hedgehog'
    trivial['equus_asinus_asinus'] = 'donkey'
    trivial['equus_caballus'] = 'horse'
    trivial['erinaceus_europaeus'] = 'european_hedgehog'
    trivial['felis_catus'] = 'cat'
    trivial['ficedula_albicollis'] = 'flycatcher'
    trivial['gadus_morhua'] = 'cod'
    trivial['gallus_gallus'] = 'chicken'
    trivial['gasterosteus_aculeatus'] = 'stickleback'
    trivial['gorilla_gorilla'] = 'gorilla'
    trivial['heterocephalus_glaber_female'] = 'naked_mole_rat_female'
    trivial['heterocephalus_glaber_male'] = 'naked_mole_rat_male'
    trivial['homo_sapiens'] = 'human'
    trivial['ictidomys_tridecemlineatus'] = 'squirrel'
    trivial['latimeria_chalumnae'] = 'coelacanth'
    trivial['lepisosteus_oculatus'] = 'spotted_gar'
    trivial['loxodonta_africana'] = 'elephant'
    trivial['macaca_mulatta'] = 'macaque'
    trivial['nannospalax_galili'] = 'blind_mole_rat'
    trivial['notamacropus_eugenii'] = 'wallaby'
    trivial['meleagris_gallopavo'] = 'turkey'
    trivial['microcebus_murinus'] = 'lemur'
    trivial['monodelphis_domestica'] = 'opossum'
    trivial['mus_musculus'] = 'mouse'
    trivial['mustela_putorius_furo'] = 'ferret'
    trivial['myotis_lucifugus'] = 'bat'
    trivial['nomascus_leucogenys'] = 'gibbon'
    trivial['ochotona_princeps'] = 'pika'
    trivial['oreochromis_niloticus'] = 'tilapia'
    trivial['ornithorhynchus_anatinus'] = 'platypus'
    trivial['oryctolagus_cuniculus'] = 'rabbit'
    trivial['oryzias_latipes'] = 'medaka'
    trivial['otolemur_garnettii'] = 'galago_lemur'
    trivial['ovis_aries'] = 'sheep'
    trivial['pan_troglodytes'] = 'chimpanzee'
    trivial['papio_anubis'] = 'baboon'
    trivial['pelodiscus_sinensis'] = 'turtle'
    trivial['petromyzon_marinus'] = 'lamprey'
    trivial['poecilia_formosa'] = 'amazon_molly'
    trivial['pongo_abelii'] = 'orangutan'
    trivial['procavia_capensis'] = 'hyrax'
    trivial['pteropus_vampyrus'] = 'flying_fox'
    trivial['pygocentrus_nattereri'] = 'piranha'
    trivial['rattus_norvegicus'] = 'rat'
    trivial['salvator_merianae'] = 'argentine_giant_tegu'
    trivial['sarcophilus_harrisii'] = 'tasmanian_devil'
    trivial['serinus_canaria'] = 'canary'
    trivial['sorex_araneus'] = 'european_shrew'
    trivial['sus_scrofa'] = 'pig'
    trivial['taeniopygia_guttata'] = 'zebra_finch'
    trivial['takifugu_rubripes'] = 'fugu'
    trivial['tetraodon_nigroviridis'] = 'pufferfish'
    trivial['tupaia_belangeri'] = 'tree_shrew'
    trivial['tursiops_truncatus'] = 'dolphin'
    trivial['vicugna_pacos'] = 'alpaca'
    trivial['xenopus_tropicalis'] = 'xenopus'
    trivial['xiphophorus_maculatus'] = 'platyfish'

    return trivial


#########################################
def make_parameter_table(cursor, table):
    print("making ", table)

    qry = "create table " + table + "  (id int(10) primary key auto_increment)"
    rows = search_db(cursor, qry, verbose=True)
    if rows:
        return False

    # make the columns
    for column in ['name', 'value']:
        qry = "alter table  %s  add  %s  varchar (50)" % (table, column)
        rows = search_db(cursor, qry, verbose=True)
        if rows:
            return False
    return False


#########################################
def make_flags_table(cursor, db_name, table):
    switch_to_db(cursor, db_name)
    table = 'flags'
    if check_table_exists(cursor, db_name, table):
        qry = "drop table " + table
        search_db(cursor, qry, verbose=True)

    qry = ""
    qry += "CREATE TABLE  %s (" % table
    qry += "     id INT NOT NULL  AUTO_INCREMENT, "
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
def make_seqregion2file_table(cursor, db_name):
    switch_to_db(cursor, db_name)
    table_name = 'seqregion2file'
    if check_table_exists(cursor, db_name, table_name):
        qry = "drop table " + table_name
        search_db(cursor, qry, verbose=False)

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
def make_taxonomy_table(cursor, db_name):
    switch_to_db(cursor, db_name)
    table_name = 'taxonomy'
    if check_table_exists(cursor, db_name, table_name):
        qry = "drop table " + table_name
        search_db(cursor, qry, verbose=False)

    qry = ""
    qry += "  CREATE TABLE  %s (" % table_name
    qry += "     id INT NOT NULL  AUTO_INCREMENT, "
    qry += "  	 name VARCHAR (100) NOT NULL, "
    qry += "  	 value TEXT NOT NULL, "

    qry += "	 PRIMARY KEY (id) "
    qry += ") ENGINE=MyISAM"

    rows = search_db(cursor, qry)
    print(qry)
    print(rows)
    return


########################################
def make_taxonomy_groups_table(cursor, db_name):
    switch_to_db(cursor, db_name)
    table_name = 'taxonomy_groups'
    if check_table_exists(cursor, db_name, table_name):
        qry = "drop table " + table_name
        search_db(cursor, qry, verbose=False)

    qry = ""
    qry += "  CREATE TABLE  %s (" % table_name
    qry += "     id INT NOT NULL  AUTO_INCREMENT, "
    qry += "  	 name VARCHAR (50) NOT NULL, "
    qry += "  	 trivial_name VARCHAR (50) NOT NULL, "
    qry += "  	 representative_species  VARCHAR (50) NOT NULL, "
    qry += "  	 members  TEXT, "
    qry += "	 PRIMARY KEY (id) "
    qry += ") ENGINE=MyISAM"

    rows = search_db(cursor, qry)
    print(qry)
    print(rows)
    return


#########################################
def feed_parameters(cursor, db_name, version):
    parameter = {}
    # in case I ever have to handle multiple versions of ensembl
    # (but for now I don't have enough space)
    # note though that there are functions in el_utils/mysql.py that assume
    # that whatever ensembl stuff is available to the mysql server corresponds to the same release
    parameter['ensembl_release_number'] = version
    parameter['blastp_e_value'] = "1.e-10"  # it will be used as a string  when fmting the blastp cmd
    parameter['min_accptbl_exon_sim'] = 0.33333  #minimum acceptable exon similarity

    table = 'parameters'
    if check_table_exists(cursor, db_name, table):
        print(table, " found in ", db_name)
    else:
        print(table, " not found in ", db_name)
        make_parameter_table(cursor, table)

    fixed_fields = {}
    update_fields = {}
    for [name, value] in parameter.items():
        fixed_fields['name'] = name
        update_fields['value'] = value
        store_or_update(cursor, table, fixed_fields, update_fields)


########################
def make_db_names_table(cursor, db_name):
    switch_to_db(cursor, db_name)
    table_name = 'db_names'
    if check_table_exists(cursor, db_name, table_name):
        qry = "drop table " + table_name
        search_db(cursor, qry, verbose=True)

    qry = ""
    qry += "  CREATE TABLE  %s (" % table_name
    qry += "     genome_db_id INT(10)  NOT NULL AUTO_INCREMENT, "
    qry += "  	 species_name VARCHAR (100) NOT NULL, "
    qry += "  	 db_name VARCHAR (100) NOT NULL, "
    qry += "	 PRIMARY KEY (genome_db_id) "
    qry += ") ENGINE=MyISAM"

    rows = search_db(cursor, qry)
    print(qry)
    print(rows)
    return


#########
def fill_db_names_table(cursor, ensembl_db_name):
    compara_db_name = get_compara_name(cursor)
    # how is this suppesd to work if they keep pulling name out of their arses like this
    exceptions = {'ovis_aries': 'ovis_aries_rambouillet',
                  'cricetulus_griseus': 'cricetulus_griseus_crigri',
                    # some ornythologists are confused
                  'cyanoderma_ruficeps': 'stachyris_ruficeps'}
    # I am not sure what use if any I have from haveing both male and female H.glaber
    # male annotation is newer
    for species_name, db_name in ensembl_db_name.items():
        name = exceptions[species_name] if species_name  in exceptions else species_name
        qry = f"select genome_db_id from {compara_db_name}.genome_db where name='{name}'"
        ret = error_intolerant_search(cursor, qry)
        if not ret:
            # Fae 2025 -  genome_db table in the ensembl ftp website was old
            # some genomes were not listed
            print(f"genome_db_id not found for {species_name}")
            continue
        genome_db_id = ret[0][0]
        print(genome_db_id, species_name, db_name)
        fixed_fields = {"genome_db_id": genome_db_id, "species_name": species_name,
                        "db_name": ensembl_db_name[species_name]}
        store_or_update(cursor, "db_names", fixed_fields=fixed_fields, update_fields={}, primary_key='genome_db_id')
        print(f"{species_name} done")
    return


#########################################
def get_trivial_names(cursor, all_species, common_name):
    tax_id = {}
    db_name = get_compara_name(cursor)
    if not db_name:
        print("compara db not found")
        exit(1)
    switch_to_db(cursor, db_name)
    for species in all_species:
        tax_id[species] = ncbi_species2taxid(cursor, species)

    trivial_names = set_trivial_names()
    for species in all_species:
        trivial = trivial_names.get(species, "not found")
        common = common_name.get(species, "not found")
        if trivial == "not found":
            if common == "not found":
                trivial = species.lower().replace(" ", "_").replace("'", "_")
            else:
                trivial = common.replace("common ", "").lower().replace(" ", "_").replace("'", "_")
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
    if species == 'homo_sapiens':
        shorthand = 'HSA'
    elif species == 'canis_lupus_dingo':
        shorthand = 'CLD'
    else:
        qry = "select meta_value  from %s.meta where meta_key='species.stable_id_prefix'" % ensembl_db
        stable_id_prefix = hard_landing_search(cursor, qry)[0][0]
        shorthand = stable_id_prefix.replace("ENS", "").replace("EiJ_", "").upper()
    return shorthand


#########################################
def tax_id_sanity_check(cursor, species, ensembl_db, ncbi_taxonomy_id):
    ncbi_taxonomy_id = int(ncbi_taxonomy_id)
    qry = "select meta_value from %s.meta where meta_key='species.taxonomy_id'" % ensembl_db
    tax_id_ensembl = int(hard_landing_search(cursor, qry)[0][0])
    if tax_id_ensembl != ncbi_taxonomy_id:
        print("Warning: tax_id mismatch for " + species)
        print(f"ensembl:{tax_id_ensembl},  ncbi:{ncbi_taxonomy_id} - will use ncbi")

    return ncbi_taxonomy_id


#########################################
def make_species_names_table(cursor, db_name):
    switch_to_db(cursor, db_name)
    table_name = 'species_names'
    if check_table_exists(cursor, db_name, table_name):
        return
    # qry = "drop table " + table_name
    # search_db(cursor, qry, verbose=True)

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


#########
def feed_species_names(cursor, exolocator_meta_db_name, all_species, ensembl_db_name, trivial_name, ncbi_tax_id):
    make_species_names_table(cursor, exolocator_meta_db_name)

    dbs_with_shorthand = {}
    for species in all_species:
        shorthand = get_shorthand(cursor, species, ensembl_db_name[species])
        tax_id = tax_id_sanity_check(cursor, species, ensembl_db_name[species], ncbi_tax_id[species])
        print("%8d %35s %45s %20s   %s (%d)" % (tax_id, species, ensembl_db_name[species],
                                                shorthand, trivial_name[species], len(trivial_name[species])))
        fixed_fields = {}
        update_fields = {}
        fixed_fields['species'] = species
        update_fields['tax_id'] = tax_id
        update_fields['shorthand'] = shorthand
        update_fields['trivial_name'] = trivial_name[species]
        switch_to_db(cursor, exolocator_meta_db_name)
        store_or_update(cursor, 'species_names', fixed_fields, update_fields)

    # sanity warning - should not happen - redo with fiexd name shorthands if it does
    for shorthand, dbs in dbs_with_shorthand.items():
        if len(dbs) == 1: continue
        print("\n dbs with same shorthand:")
        print(shorthand, dbs)
        print("This should not have happened. Redo with fiexd name shorthands.")
        exit()


#########################################

#########################################
def main():
    # MySQL connection parameters from .env
    cursor = mysql_server_connect(user=os.getenv('MYSQL_USER'),
                                  passwd=os.getenv('MYSQL_PASSWORD'),
                                  host=os.getenv('MYSQL_HOST', 'localhost'),
                                  port=int(os.getenv('MYSQL_PORT', 3306)))

    error_intolerant_search(cursor, "set autocommit=1")

    #######################################################
    # check if the config db exists -- if not, make it
    exolocator_meta_db_name = "exolocator_meta"
    qry = "show databases like '%s'" % exolocator_meta_db_name
    rows = search_db(cursor, qry)
    if not rows:
        print(exolocator_meta_db_name, "database not found - creating a new one")
        qry = f"create database {exolocator_meta_db_name}"
        rows = search_db(cursor, qry)
        if rows:
            print("some problem creating the database ...")
            search_db(cursor, qry, verbose=True)
        print(exolocator_meta_db_name, "database ok")

    switch_to_db(cursor, exolocator_meta_db_name)
    #######################################################
    # store parameters, to make sure we are consistent and reproducible
    # this has to be stored first, because from this point on we read version from the parameter table
    feed_parameters(cursor, exolocator_meta_db_name, Config.release_number)
    #######################################################

    [all_species, ensembl_db_name, common_name] = get_ensembl_species(cursor)
    # for species in all_species:
    #     print(species, "     ", common_name[species], "     ", ensembl_db_name[species])


    #######################################################
    # create db_names table , mapping species id to db_name
    make_db_names_table(cursor, exolocator_meta_db_name)
    fill_db_names_table(cursor, ensembl_db_name)

    ########################################################
    # create flags table (for flagging arbitrary problems - to be filled as we go)
    make_flags_table(cursor, exolocator_meta_db_name, 'flags')

    ########################################################
    # create table to store tax info that we might need and don't want to recalculate each time
    make_taxonomy_table(cursor, exolocator_meta_db_name)
    make_taxonomy_groups_table(cursor, exolocator_meta_db_name)  # fish, frogs, mammals etc


    #######################################################
    # get trivial names
    [trivial_names, ncbi_tax_id] = get_trivial_names(cursor, all_species, common_name)

    #######################################################
    # add species shorthands (used in ENS* names formation)
    # though we will not need it until the paralogue alignment reconstruction point)
    feed_species_names(cursor, exolocator_meta_db_name, all_species, ensembl_db_name, trivial_names, ncbi_tax_id)

    mysql_server_conn_close(cursor)


#########################################
if __name__ == '__main__':
    main()
