#!/usr/bin/python3 -u
import os.path

from dotenv import load_dotenv

from config import Config

from el_utils.ensembl import  *
from el_utils.processes import *


# Load environment variables
load_dotenv()


####################################################
def prioritize_versioning(fastasfile):
    headerct = 0
    index_versioning = 0
    with open(fastasfile) as inf:
        for line in inf:
            if not '>' in line: continue
            headerct +=1
            name = line.split()[0]
            if "." in name:
                index_versioning += 1
            if headerct==100: break
    if index_versioning/headerct>0.5:
        return [f".{v}" for v in range(1,13)] + [""]
    else:
        return [""] + [f".{v}" for v in range(1,13)]


###########
def canon_seqs_for_species(cursor, species, ensembl_db_name):

    print(f"#######################\n{species}")

    # we'll do this for human only, for now
    switch_to_db (cursor, ensembl_db_name[species])

    fasta_path = f"{Config.fasta_repo}/{species}/pep"
    ext = "pep.all.fa.gz"
    in_fastas = [fnm for fnm in os.listdir(fasta_path) if fnm[-(len(ext)):] == ext]
    if len(in_fastas) == 0:
        print(f"*.{ext} not found")
        exit()
    if len(in_fastas) > 1:
        print(f"multiple *.{ext}  found:", in_fastas)
        exit()
    gzipped_in_fasta = f"{fasta_path}/{in_fastas[0]}"
    in_fasta = gzipped_in_fasta[:-3]
    run_subprocess(f"gunzip {gzipped_in_fasta}")

    # this is very unsystematic, but some species such as Mus spretus do not have the version number
    # so in that case look for that name format first
    versioning = prioritize_versioning(in_fasta)

    out_fasta = f"{fasta_path}/canonical_peptides.fa"
    if os.path.exists(out_fasta): os.remove(out_fasta)

    #######################
    print(f"{species} collecting peptide ids")
    time0 = time()
    stable_translation_ids = []
    # get_gene_ids(cursor, db_name=None, biotype=None,  stable=False, ref_only=False):
    for gene_id in get_gene_ids(cursor, ensembl_db_name[species], biotype="protein_coding"):
        canonical_translation_stable_id = gene2stable_canon_transl_id(cursor, gene_id)
        if canonical_translation_stable_id: stable_translation_ids.append(canonical_translation_stable_id)
    mins = (time()-time0)/60
    print(f"{species} collecting peptide ids done in %.1f min" % mins)

    #######################
    print(f"{species} extracting sequences")
    time0 = time()
    names_file = f"{species}.names.txt"
    outf = open(names_file, "w")
    for stable_translation_id in stable_translation_ids:
        for version in versioning:
            print(f"{stable_translation_id}{version}", file=outf)
            break  # one version is enough

    cmd = f"{Config.seqtk} subseq {in_fasta} {names_file}"
    run_subprocess(cmd, stdoutfnm=out_fasta)

    secs = (time()-time0)
    print(f"{species} extracting seqs done in {secs:.1f} secS")

    run_subprocess(f'gzip {in_fasta}')
    if os.path.exists(names_file): os.remove(names_file)


###########
def canon_seqs(species_chunk, other_args):
    [ensembl_db_name] = other_args
    cursor = mysql_server_connect(user=os.getenv('MYSQL_USER'),
                                  passwd=os.getenv('MYSQL_PASSWORD'),
                                  host=os.getenv('MYSQL_HOST', 'localhost'),
                                  port=int(os.getenv('MYSQL_PORT', 3306)))

    for species in species_chunk:
        canon_seqs_for_species(cursor, species, ensembl_db_name)

    mysql_server_conn_close(cursor)


####################################################
def main():
    # MySQL connection parameters from .env
    cursor = mysql_server_connect(user=os.getenv('MYSQL_USER'),
                                  passwd=os.getenv('MYSQL_PASSWORD'),
                                  host=os.getenv('MYSQL_HOST', 'localhost'),
                                  port=int(os.getenv('MYSQL_PORT', 3306)))

    [all_species, ensembl_db_name] = get_species(cursor)
    mysql_server_conn_close(cursor)

    ok_species = []
    for species in all_species:
        species_path = f"{Config.fasta_repo}/{species}"
        if not os.path.exists(species_path):
            print(f"{species_path} not found")
        elif not os.path.exists(pep_path := f"{species_path}/pep"):
            print(f"{pep_path} not found")
        else:
            ok_species.append(species)

    number_of_chunks = 8

    print(f"number of species: {len(ok_species)}, number of chunks: {number_of_chunks}")
    parallelize(number_of_chunks, canon_seqs, ok_species, [ensembl_db_name])


#####################################################
if __name__=="__main__":
    main()