#!/usr/bin/python3

import os

from config import Config

from el_utils.ensembl       import  *


####################################################
def store_seq_filenames (cursor, name, file_names):
    fixed_fields  = {}
    update_fields = {}
    fixed_fields ['name']      = name
    update_fields['file_name'] = file_names
    retval = store_or_update (cursor, "seq_region", fixed_fields, update_fields, primary_key='seq_region_id')
    return retval

####################################################
def main():

    db = connect_to_mysql(Config.mysql_conf_file)

    cursor = db.cursor()
    fasta_path = "/storage/databases/ensembl-{}/fasta".format(Config.release_number)

    [all_species, ensembl_db_name] = get_species (cursor)

    for species in all_species:
        print(species)
        dna_path = "{0}/{1}/dna".format(fasta_path, species)
        if (not os.path.exists(dna_path)):
            print("problem:", dna_path, "not found")
            exit(1)

        fasta_files = []
        for r,d,files in os.walk(dna_path):
            for file in files:
                if (not file[-3:] == ".fa"):
                    continue
                fasta_files.append(file)
        
        name2file = {}
        for file in fasta_files:
            #print(dna_path, file)
            cmd = "grep '>' {0}/{1}".format(dna_path, file)
            ret = subprocess.getoutput(cmd)
            headers = ret.split("\n")
            #print("number of headers: ", len(headers))
            for hdr in headers:
                fields = hdr.split(" ")
                name = fields[0].replace (">", "")
                #print name
                if (name not in name2file):
                    name2file[name] = []
                name2file[name].append(file)

        qry = "use "+ensembl_db_name[species]
        search_db (cursor, qry)

        for name in name2file.keys():
            file_names = " ".join(name2file[name])
            #print(name, file_names)
            store_seq_filenames (cursor, name, file_names)

    cursor.close()
    db    .close()
    


####################################################
if __name__ == '__main__':
    main()
