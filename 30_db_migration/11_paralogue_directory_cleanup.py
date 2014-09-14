#!/usr/bin/python -u
# make the best alignment we can using the maps
# we currently have at hand

import pdb
#pdb.set_trace()

import MySQLdb, commands, re, os, sys

from el_utils.mysql   import  connect_to_mysql, connect_to_db
from el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from el_utils.ensembl import  *
from el_utils.el_specific   import  *

from el_utils.utils   import  erropen, output_fasta, input_fasta, parse_aln_name
from el_utils.map     import  Map, get_maps
from el_utils.tree    import  species_sort
from el_utils.ncbi    import  taxid2trivial
from el_utils.almt_cmd_generator import AlignmentCommandGenerator
from el_utils.config_reader      import ConfigurationReader
from el_utils.translation        import phase2offset, translation_bounds, crop_dna, translate
from el_utils.processes import parallelize
from el_utils.custom  import get_theme_ids
from bitstring import Bits
from alignment import * # C implementation of smith waterman
from   random  import choice
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna


#########################################
def check_directory (cfg, species, species_id, pep_or_dna):
    
    directory = "{0}/para/{1}/{2}".format(cfg.dir_path['afs_dumps'], species_id, pep_or_dna)
    if not os.path.exists(directory):
        print "warning: ", directory, "not found"
        return ""
    return directory

#########################################
#########################################
#########################################
#########################################
def make_alignments (species_list, db_info):

    [local_db, ensembl_db_name] = db_info

    verbose      = False
    flank_length = 10

    if local_db:
        db     = connect_to_mysql()
        cfg    = ConfigurationReader()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    # find db ids adn common names for each species db
    [all_species, ensembl_db_name] = get_species (cursor)

    max_days = 60

    #species_list.reverse()
    for species in species_list:

        pep_produced = 0
        dna_produced = 0
        has_paralogues = 0

        species_shorthand = get_species_shorthand(cursor, species)
        print species, species_shorthand

        directory = check_directory (cfg, species, species_shorthand, "pep")
        if not directory: continue

        removed   = 0
        remaining = 0
        for dirname, dirnames, filenames in os.walk(directory):
            for filename in filenames:
                full_name =  os.path.join(dirname, filename)
                time_modified = os.path.getmtime(full_name)
                number_of_days_since_modified = (time.time() - time_modified)/(60*60*24)
                if number_of_days_since_modified > max_days:
                    #print "removing", filename, "made", number_of_days_since_modified, "ago"
                    os.remove(full_name)
                else:
                    remaining += 1
        print species, "done, removed", removed, "files, remaining", remaining
        
    

#########################################
def main():
    
    no_threads = 1

    local_db = False

    if local_db:
        db = connect_to_mysql()
    else:
        db = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    parallelize (no_threads, make_alignments, all_species, [local_db, ensembl_db_name])
    
    return True


#########################################
if __name__ == '__main__':
    main()

