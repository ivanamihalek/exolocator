#!/usr/bin/python3 -u
import os.path

from dotenv import load_dotenv

from el_utils import mysql
from el_utils.mysql import *
from el_utils.ensembl import *
from config import Config
from el_utils.processes import run_subprocess
from el_utils.tree import species_sort
from el_utils.utils import die_if_not_dir, die_if_not_nonzero_file


def write_to_fasta(home, species, stable_transl_id,  out_fasta):
    fasta_path = f"{Config.fasta_repo}/{species}/pep"
    if not os.path.exists(fasta_path):
        print(f"{fasta_path} not found")
        return False
    os.chdir(fasta_path)

    canonical_fasta = f"canonical_peptides.fa"
    if not os.path.exists(canonical_fasta):
        print(f"{canonical_fasta} not found in {fasta_path}")
        return False
    if os.path.getsize(canonical_fasta) == 0:
        print(f"{canonical_fasta} is empty in  {fasta_path}")
        return False
    print(os.getcwd())

    names_file = f"{Config.scratch}/names.{os.getpid()}.txt"
    tmp_fasta  = f"{Config.scratch}/tmp_fasta.{os.getpid()}.txt"
    with open(names_file, "w") as outnm:
        print(stable_transl_id, file=outnm)
    cmd = f"{Config.seqtk} subseq {canonical_fasta} {names_file}"
    run_subprocess(cmd, stdoutfnm=tmp_fasta)
    os.remove(names_file)

    os.chdir(home)

    if os.path.exists(tmp_fasta) and os.path.getsize(tmp_fasta) > 0:
        # name = "_".join([sp[:3] for sp in species.split("_")]) # +[gene_name])
        name = species
        with open(tmp_fasta, "r") as inf, open(out_fasta, "a") as outnm:
            for line in inf:
                if line[0] == ">":
                    print(f">{name}", file=outnm)
                else:
                    print(line.strip(), file=outnm)
        os.remove(tmp_fasta)
        return True
    else:
        print(f"in write_to_fasta(): {tmp_fasta} does not exist or is empty")
        return False


def reorder_seqs(in_afa, species_sorted, out_afa, trivial=None, name_prefix=None):
    seq = {}
    name = 'anon'
    with open(in_afa) as inf:
        for line in inf:
            line = line.strip()
            if len(line) == 0: continue
            if line[0] == '>':
                name = line[1:]
                seq[name] = ""
            else:
                seq[name] += line

    with open(out_afa, "w") as outf:
        for name in species_sorted:
            outname = trivial[name] if trivial else name
            if name_prefix: outname = f"{name_prefix}_{outname}"
            outf.write(f">{outname}\n")
            for i in range(0, len(seq[name]), 100):
                outf.write(f"{seq[name][i:i + 100]}\n")


#########################################
def main():

    if len(sys.argv) < 2:
        print("usage: %s <gene symbol> [trivial] [prepend]" % sys.argv[0])
        print("trivial = use trivial species name; prepend = prepend gene name")
        exit()

    gene_name = sys.argv[1]
    trivial = "trivial" in sys.argv
    prepend = "prepend" in sys.argv  # prepends geen synbol to gene name
    ref_species = 'homo_sapiens'  # the orthologue table is filled only here, for the moment

    out_fasta = f"{gene_name}.orthos.fasta"
    out_afa = f"{gene_name}.orthos.afa"
    for fnm in [out_fasta, out_afa]:
        if os.path.exists(fnm): os.remove(fnm)

    home = os.getcwd()

    cursor = mysql_using_env_creds()

    qry = f"select ensembl_gene_id  from identifier_maps.hgnc where approved_symbol='{gene_name}'"
    ensembl_stable_gene_id = hard_landing_search(cursor, qry)[0][0]

    [all_species, ensembl_db_name] = get_species(cursor)
    # ? species_sort(cursor, all_species, 'homo_sapiens')

    switch_to_db(cursor, ensembl_db_name[ref_species])
    qry = "select gene_id from gene where stable_id='%s'" % ensembl_stable_gene_id
    gene_id = hard_landing_search(cursor, qry)[0][0]

    ref_stable_transl_id = gene2stable_canon_transl_id(cursor, gene_id, ensembl_db_name[ref_species])

    write_to_fasta(home, ref_species, ref_stable_transl_id, out_fasta)

    species_in_the_almt = [ref_species]
    qry = "select  cognate_gene_id, cognate_genome_db_id from orthologues where gene_id=%d" % gene_id
    for line in list(error_intolerant_search(cursor, qry)):
        [cognate_gene_id, cognate_genome_db_id] = line
        qry = f"select db_name from exolocator_meta.db_names where genome_db_id={cognate_genome_db_id}"
        db_name = hard_landing_search(cursor, qry)[0][0]
        stable_transl_id = gene2stable_canon_transl_id(cursor, cognate_gene_id, db_name)
        species = db_name.split("core")[0].rstrip("_")
        if species not in all_species: continue
        print(db_name, species, cognate_gene_id, stable_transl_id)
        ok = write_to_fasta(home, species, stable_transl_id, out_fasta)
        print(ok)
        if ok: species_in_the_almt.append(species)

    cmd = f"{Config.muscle} -in {out_fasta} -out tmp.afa"
    subprocess.call(["bash", "-c", cmd])

    species_sorted = species_sort(cursor, species_in_the_almt, ref_species)
    trivial_names = get_trivial(cursor, species_sorted) if trivial else None
    name_prefix = gene_name if prepend else None
    reorder_seqs('tmp.afa', species_sorted, out_afa, trivial_names, name_prefix)
    if os.path.exists('tmp.afa'): os.remove('tmp.afa')

    mysql_server_conn_close(cursor)

    return True


#########################################
if __name__ == '__main__':
    main()
