#!/usr/bin/python


import MySQLdb
import commands, os
from   el_utils.mysql         import  connect_to_mysql, search_db, store_or_update
from   el_utils.config_reader import  ConfigurationReader
from   el_utils.ensembl       import  *


####################################################
def store_seq_filenames (cursor, name, file_names):
    fixed_fields  = {}
    update_fields = {}
    fixed_fields ['name']      = name
    update_fields['file_name'] = file_names
    retval = store_or_update (cursor, "seq_region", fixed_fields, update_fields)
    return retval

####################################################
def main():

    db = connect_to_mysql()
    cr = ConfigurationReader()

    cursor = db.cursor()
    fasta_path = cr.get_path('ensembl_fasta')

    [all_species, ensembl_db_name] = get_species (cursor)

    for species in all_species:
    #for species in ['danio_rerio']:
        print species
        dna_path = "{0}/{1}/dna".format(fasta_path, species)
        if (not os.path.exists(dna_path)):
            print "problem:", dna_path, "not found"
            exit(1)

        fasta_files = []
        for r,d,files in os.walk(dna_path):
            for file in files:
                if (not file[-3:] == ".fa"):
                    continue
                fasta_files.append(file)
        
        name2file = {}
        for file in fasta_files:
            print dna_path, file
            cmd = "grep '>' {0}/{1}".format(dna_path, file)
            ret = commands.getoutput(cmd)
            headers = ret.split("\n")
            print "number of headers: ", len(headers)
            for hdr in headers:
                fields = hdr.split(" ")
                name = fields[0].replace (">", "")
                #print name
                if (not name2file.has_key(name)):
                    name2file[name] = []
                name2file[name].append(file)

        qry = "use "+ensembl_db_name[species]
        search_db (cursor, qry)

        for name in name2file.keys():
            file_names = ""
            for file in  name2file[name]:
                if file_names:
                    file_names += " "
                file_names += file
            store_seq_filenames (cursor, name, file_names)
 
    cursor.close()
    db    .close()
    


####################################################
if __name__ == '__main__':
    main()
