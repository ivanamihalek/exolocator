#!/usr/bin/python3 -u
import os.path
from glob import glob

from dotenv import load_dotenv

from config import Config

from el_utils.ensembl import  *
from el_utils.processes import *


# Load environment variables


####################################################
# why am I doing this?
def extract_peptide_ids(fastafile) -> dict:
    # for some reason, Ensembl is tacking the versionin to the stable id in the pep.fa file
    # e.g ENSP00000452345.1, not just ENSP00000452345
    # the search for stable id, however will return just the base part ENSP00000452345
    versions_in_fastafile = {}
    with open(fastafile) as inf:
        for line in inf:
            if line[0] != '>': continue
            name = line.split()[0][1:]  # get rid og the ">" thing
            if "." in name:
                name, version = name.split(".")
                version = "." + version
            else:
                version = ""
            if name not in versions_in_fastafile: versions_in_fastafile[name] = []
            versions_in_fastafile[name].append(version)

    if not versions_in_fastafile:
        print(f"no seq names in {fastafile} (?!)")
        exit()

    return versions_in_fastafile


def header_cleanup(complicated_header_fasta):
    tmp_out = complicated_header_fasta + ".tmp"
    with open(tmp_out, "w") as outf:
        with open(complicated_header_fasta) as inf:
            for line in inf:
                if line[0] == '>':
                    name = line.split()[0][1:]
                    if "." in name: name = name.split(".")[0]
                    outf.write(">" + name + "\n")
                else:
                    outf.write(line)
    print(f"header cleanup for {complicated_header_fasta} done")
    os.rename(tmp_out, complicated_header_fasta)


###########
def canon_seqs_for_species(cursor, species, ensembl_db_name):

    print(f"#######################\n{species}")

    # we'll do this for human only, for now
    switch_to_db (cursor, ensembl_db_name[species])

    fasta_path = f"{Config.fasta_repo}/{species}/pep"
    out_fasta = f"{fasta_path}/canonical_peptides.fa"
    print(fasta_path)
    # clean blast indexing files
    camefrom = os.getcwd()
    os.chdir(fasta_path)
    for fnm in glob("canonical_peptides.*"):
        print(f"removing {fnm}")
        os.remove(fnm)
    os.chdir(camefrom)

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
    name_versions_in_fasta_file  = extract_peptide_ids(in_fasta)

    #######################
    print(f"{species} collecting peptide ids")
    time0 = time()
    stable_translation_ids = []
    # get_gene_ids(cursor, db_name=None, biotype=None,  stable=False, ref_only=False):
    for gene_id in get_gene_ids(cursor, ensembl_db_name[species], biotype="protein_coding"):
        canonical_translation_stable_id = gene2stable_canon_transl_id(cursor, gene_id)
        if not canonical_translation_stable_id: continue
        if canonical_translation_stable_id not in name_versions_in_fasta_file:
            print(canonical_translation_stable_id, f"not found in {fasta_path} (?!)")
            exit()
        # e.g ENSP00000452040.2
        # just take the fist version that we bump into - I am not sure what to do with it anyway
        can_trsl_stable_id_w_version = canonical_translation_stable_id + name_versions_in_fasta_file[canonical_translation_stable_id][0]
        stable_translation_ids.append(can_trsl_stable_id_w_version)
    mins = (time()-time0)/60
    print(len(stable_translation_ids))
    print(f"{species} collecting canonical peptide ids done in {mins:.1f} min. Collected {len(stable_translation_ids)} ids." )
    print("some examples:", stable_translation_ids[:5])

    #######################
    print(f"{species} extracting sequences")
    time0 = time()
    names_file = f"{species}.names.txt"
    with open(names_file, "w") as outf:
        [print(stable_translation_id, file=outf) for stable_translation_id in stable_translation_ids]

    cmd = f"{Config.seqtk} subseq {in_fasta} {names_file}"
    run_subprocess(cmd, stdoutfnm=out_fasta)

    secs = (time()-time0)
    print(f"{species} extracting seqs done in {secs:.1f} secs")

    run_subprocess(f'gzip {in_fasta}')

    # now we've ended up with the versioning numbers in the canonical peptide file
    header_cleanup(out_fasta)

    if os.path.exists(names_file): os.remove(names_file)

###########
def canon_seqs(species_chunk, other_args):
    [ensembl_db_name] = other_args
    cursor = mysql_using_env_creds()

    for species in species_chunk:
        canon_seqs_for_species(cursor, species, ensembl_db_name)

    mysql_server_conn_close(cursor)


####################################################
def main():
    # MySQL connection parameters from .env
    cursor = mysql_using_env_creds()

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