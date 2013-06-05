#!/usr/bin/python

""" 
@package 01_make_config_db.py
@brief  Creating ad filling the configuration database. Hardcoded stuff is here.
@author Ivana
@date   2012
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
import MySQLdb
import os, sys
from   el_utils.mysql   import  connect_to_mysql, search_db
from   el_utils.mysql   import  check_table_exists, store_or_update
from   el_utils.ensembl import  get_species, get_compara_name, species2taxid
from   el_utils.ncbi    import  get_ncbi_tax_name


#########################################
def  make_parameter_table (cursor):

    """
    Creates parameter table in the config database.
    @param [cursor] db cursor, assumed top be pointing to the config database
    @retval True  on success
    @retval False on failure;  in that case the seach_db() call is repeated in verbose mode.
    """

    table = 'parameter'

    print "making ", table

    qry  = "create table " + table + "  (id int(10) primary key auto_increment)"
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        return False

    # make the columns
    for column  in  ['name', 'value']:
        qry = "alter table  %s  add  %s  varchar (20)" % (table, column)
        rows = search_db (cursor, qry, verbose=True)
        if (rows):
            return False
    return False
 
#########################################
def  make_path_table (cursor, table):

    print "making ", table

    qry  = "create table " + table + "  (id int(10) primary key auto_increment)"
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        return False

    # make the columns
    column  = 'name'
    qry = "alter table  %s  add  %s  varchar (20)" % (table, column)
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        return False

    column  = 'path'
    qry = "alter table  %s  add  %s  blob" % (table, column)
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        return False


#########################################
def  make_seqregion2file_table (cursor):

    table = 'seqregion2file'

    
    qry  = "create table " + table + "  (seqregion_id int(10) primary key)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    # make the columns
    column  = 'seq_name'
    qry = "alter table  %s  add  %s  varchar (100)" % (table, column)
    rows = search_db (cursor, qry)
    if (rows):
        return False

    column  = 'file_name'
    qry = "alter table  %s  add  %s  blob" % (table, column)
    rows = search_db (cursor, qry)
    if (rows):
        return False


#########################################
def make_table (cursor, table):
    
    if   table == 'util_path':
        make_path_table (cursor, table)
    elif table == 'dir_path':
        make_path_table (cursor, table)
    elif table == 'parameter':
        make_parameter_table (cursor)
    elif table == 'seqregion2file':
        make_seqregion2file_table (cursor)
       
    else:
        print "I don't know how to make table '%s'" % table



#########################################
def  feed_trivial_names (cursor, all_species):

    tax_id  = {}
    trivial = {}

    trivial['ailuropoda_melanoleuca']='panda' 
    trivial['anolis_carolinensis']='anole_lizard' 
    trivial['bos_taurus']='cow' 
    trivial['callithrix_jacchus']='marmoset' 
    trivial['canis_familiaris']='dog' 
    trivial['cavia_porcellus']='guinea_pig' 
    trivial['choloepus_hoffmanni']='sloth' 
    trivial['danio_rerio']='zebrafish' 
    trivial['dasypus_novemcinctus']='armadillo' 
    trivial['dipodomys_ordii']='kangaroo_rat' 
    trivial['echinops_telfairi']='madagascar_hedgehog' 
    trivial['equus_caballus']='horse' 
    trivial['erinaceus_europaeus']='european_hedgehog' 
    trivial['felis_catus']='cat' 
    trivial['gadus_morhua']='cod' 
    trivial['gallus_gallus']='chicken' 
    trivial['gasterosteus_aculeatus']='stickleback' 
    trivial['gorilla_gorilla']='gorilla' 
    trivial['homo_sapiens']='human' 
    trivial['ictidomys_tridecemlineatus']  ='squirrel' 
    trivial['latimeria_chalumnae']         ='coelacanth' 
    trivial['loxodonta_africana']          ='elephant' 
    trivial['macaca_mulatta']              ='macaque' 
    trivial['macropus_eugenii']='wallaby' 
    trivial['meleagris_gallopavo']='turkey' 
    trivial['microcebus_murinus']='lemur' 
    trivial['monodelphis_domestica']='opossum' 
    trivial['mus_musculus']='mouse' 
    trivial['mustela_putorius_furo']='ferret' 
    trivial['myotis_lucifugus']='bat' 
    trivial['nomascus_leucogenys']='gibbon' 
    trivial['ochotona_princeps']='pika' 
    trivial['oreochromis_niloticus']='tilapia' 
    trivial['ornithorhynchus_anatinus']='platypus' 
    trivial['oryctolagus_cuniculus']='rabbit' 
    trivial['oryzias_latipes']='medaka' 
    trivial['otolemur_garnettii']='galago_lemur' 
    trivial['pan_troglodytes']='chimpanzee' 
    trivial['pelodiscus_sinensis']='turtle' 
    trivial['petromyzon_marinus']='lamprey' 
    trivial['pongo_abelii']='orangutan' 
    trivial['procavia_capensis']='hyrax' 
    trivial['pteropus_vampyrus']='flying_fox' 
    trivial['rattus_norvegicus']='rat' 
    trivial['sarcophilus_harrisii']='tasmanian_devil' 
    trivial['sorex_araneus']='european_shrew' 
    trivial['sus_scrofa']='pig' 
    trivial['taeniopygia_guttata']='zebra_finch' 
    trivial['takifugu_rubripes']='fugu' 
    trivial['tarsius_syrichta']='tarsier' 
    trivial['tetraodon_nigroviridis']='pufferfish' 
    trivial['tupaia_belangeri']='tree_shrew' 
    trivial['tursiops_truncatus']='dolphin' 
    trivial['vicugna_pacos']='alpaca' 
    trivial['xenopus_tropicalis']='xenopus' 
    trivial['xiphophorus_maculatus']='platyfish' 


    db_name = get_compara_name (cursor)
    if (not db_name):
        print "compara db not found"
        exit(1)
    qry = "use %s " % db_name
    search_db (cursor, qry)
    for species in all_species:
        tax_id[species] = species2taxid (cursor, species)

    # switch to ncbi taxonomy database
    db_name = get_ncbi_tax_name (cursor)
    if (not db_name):
        print "ncbi taxonomy db not found"
        exit(1)

    qry = "use %s " % db_name
    search_db (cursor, qry)
    for species in all_species:
        if trivial.has_key(species):
            fixed_fields  = {}
            update_fields = {}
            fixed_fields ['tax_id']     = tax_id[species]
            fixed_fields ['name_class'] = 'trivial'
            update_fields['name_txt']   = trivial[species]
            store_or_update (cursor, 'names', fixed_fields, update_fields)
        else:
            print "trivial for ", species, " not found "
            trivial[species] = ""




    return True

#########################################
def main():

    

    dir_path = {}
    dir_path['ensembl_fasta'] = '/mnt/ensembl-mirror/release-69/fasta'
    # local juggling of data from one database base to the other
    dir_path['afs_dumps']     = '/afs/bii.a-star.edu.sg/dept/biomodel_design/Group/ivana/'
    dir_path['afs_dumps']    += 'ExoLocator/results/dumpster'
    dir_path['resources']     = '/afs/bii.a-star.edu.sg/dept/biomodel_design/Group/ivana/'
    dir_path['resources']    += 'pypeworks/exolocator/resources'
    dir_path['scratch']       = '/tmp'
    dir_path['maxentscan']    = '/afs/bii.a-star.edu.sg/dept/biomodel_design/Group/ivana/'
    dir_path['maxentscan']   += 'pypeworks/exolocator/pl_utils/maxentscan'

    util_path = {}
    util_path['mafft']    = '/usr/bin/mafft'
    util_path['blastall'] = '/usr/bin/blastall'
    util_path['fastacmd'] = '/usr/bin/fastacmd'
    util_path['sw#']      = '/usr/bin/swsharp'
    util_path['usearch']  = '/usr/bin/usearch'
    util_path['score3']   = dir_path['maxentscan'] + '/score3.pl'
    util_path['score5']   = dir_path['maxentscan'] + '/score5.pl'

    parameter = {}
    parameter['blastp_e_value'] = "1.e-10" # it will be used as a string  when fmting the blastp cmd
    parameter['blosum_hacked']  = "blosum_hacked.txt" # filename, to be found in resources
    parameter['min_accptbl_exon_sim'] = 0.33333 #minimum acceptable exon similarity

    # check if the paths are functioning (at this point at least)
    for util in util_path.values():
        if (not os.path.exists(util)):
            print util, " not found "
            sys.exit (1)

    for dir in dir_path.values():
        if (not os.path.exists(dir)):
            print dir, " not found "
            sys.exit (1)
        if (not os.path.isdir (dir)):
            print dir, " is not a directory "
            sys.exit (1)
            
    #db     = connect_to_mysql()
    #db     = connect_to_mysql(user="marioot", passwd="tooiram")
    db      = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()


    #######################################################
    # check if the config db exists -- if not, make it
    db_name   = "exolocator_config"
    qry  = "show databases like'%s'" % db_name
    rows = search_db (cursor, qry)
    if (not rows):
        print db_name, "database not found"
        qry = "create database %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            print "some problem creating the database ..."
            rows = search_db (cursor, qry, verbose = True)
    else:
        print db_name, "database found"

    qry = "use %s " % db_name
    search_db (cursor, qry)
        
    # make tables
    for table in ['util_path', 'dir_path', 'parameter']:
        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
        else:
            print table, " not found in ", db_name
            make_table (cursor, table)
   
    # fill util, dir and path tables 
    fixed_fields  = {}
    update_fields = {}
    for [name, path] in util_path.iteritems():
        fixed_fields['name']  = name
        update_fields['path'] = path
        store_or_update (cursor, 'util_path', fixed_fields, update_fields)

    fixed_fields  = {}
    update_fields = {}
    for [name, path] in dir_path.iteritems():
        fixed_fields['name'] = name
        update_fields['path'] = path
        store_or_update (cursor, 'dir_path', fixed_fields, update_fields)

    fixed_fields  = {}
    update_fields = {}
    for [name, value] in parameter.iteritems():
        fixed_fields['name']  = name
        update_fields['value'] = value
        store_or_update (cursor, 'parameter', fixed_fields, update_fields)

    #######################################################
    # add trivial names to ncbi_taxonomy.names
    [all_species, ensembl_db_name] = get_species (cursor)
    feed_trivial_names (cursor, all_species)


    cursor.close()
    db.close()



#########################################
if __name__ == '__main__':
    main()
