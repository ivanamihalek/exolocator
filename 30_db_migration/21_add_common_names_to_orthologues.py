#!/usr/bin/python

import MySQLdb, glob
import os, commands, sys
from   time import clock, time 
from   el_utils.mysql         import *
from   el_utils.config_reader import ConfigurationReader
from   el_utils.utils         import erropen
from   el_utils.threads       import  parallelize

##########################################
#########################################
def  trivial_names ():
    trivial = {}

    trivial['ailuropoda_melanoleuca']= 'panda' 
    trivial['anas_platyrhynchos']    = 'duck'
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
    trivial['erinaceus_europaeus'] = 'european_hedgehog' 
    trivial['felis_catus'] = 'cat' 
    trivial['ficedula_albicollis'] ='flycatcher'
    trivial['gadus_morhua']='cod' 
    trivial['gallus_gallus']='chicken' 
    trivial['gasterosteus_aculeatus']='stickleback' 
    trivial['gorilla_gorilla']='gorilla' 
    trivial['homo_sapiens']='human' 
    trivial['ictidomys_tridecemlineatus']  ='squirrel' 
    trivial['spermophilus_tridecemlineatus']  ='squirrel' 
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
    trivial['sorex_araneus']          = 'european_shrew' 
    trivial['sus_scrofa']             = 'pig' 
    trivial['taeniopygia_guttata']    = 'zebra_finch' 
    trivial['takifugu_rubripes']      = 'fugu' 
    trivial['tarsius_syrichta']       = 'tarsier' 
    trivial['tetraodon_nigroviridis'] = 'pufferfish' 
    trivial['tupaia_belangeri']       = 'tree_shrew' 
    trivial['tursiops_truncatus']     = 'dolphin' 
    trivial['vicugna_pacos']          = 'alpaca' 
    trivial['xenopus_tropicalis']     = 'xenopus' 
    trivial['xiphophorus_maculatus']  = 'platyfish' 

    return trivial

##########################################
def add_common_name_column (cursor):
    
    qry  = "ALTER TABLE ortholog  ADD common_name VARCHAR(50)"
    rows = search_db (cursor, qry)
    if (rows):
        return False

    return True


##########################################
def add_common_names (cursor, trivial):

    
    qry  = "select id, species from ortholog"
    rows = search_db (cursor, qry)
    for row in rows:
        if len(row) <2:
            print row
            exit(1)
        [id, species] = row
        qry  =  "update ortholog set common_name = '%s'" % trivial[species]
        qry +=  " where id = %s " % id
        ret  = search_db (cursor, qry)
        if (ret):
           search_db (cursor, qry, verbose=True) 
           exit (1)



#########################################
def main():

    
    no_threads = 1
    
    db_name =  "exolocator_db"
    db      = connect_to_mysql(user="marioot", passwd="tooiram")
    cursor  = db.cursor()
    switch_to_db (cursor, db_name)
    ###############
    trivial = trivial_names()
    #add_common_name_column (cursor)   
    add_common_names (cursor, trivial)

    ###############
    cursor.close()
    db    .close()


#########################################
if __name__ == '__main__':
    main()

