#!/usr/bin/python3 -u
import os

from sys import  stderr
from el_utils.ensembl import *
from config import Config
from el_utils.utils import die_if_not_nonzero_file


# https://www.uniprot.org/help/canonical_and_isoforms
# "he various UniProtKB distribution formats (Flat Text, XML, RDF/XML) display only the canonical sequence. "


def read_primary_accession_numbers():
    acc_number_table = "/storage/databases/uniprot/sec2primary_acc_id_map.txt"
    die_if_not_nonzero_file(acc_number_table)
    primary = {}
    with open(acc_number_table) as inf:
        reading = False
        trigger = "Secondary AC"
        for line in inf:
            if not reading and line[:len(trigger)] == trigger:
                reading = True
            elif reading:
                if not line[0].isalnum(): continue
                [sec, prim] = line.strip().split()
                primary[sec] = prim
    return primary


def check_primary(primary, uniprot_ids):
    other_ids = {}
    uids = uniprot_ids.replace(' ', '').split(",")
    # some genes  can code for completely different proteins eg. CDKNA2 gene coding for CDKNA2 and ARF proteins
    # https://www.uniprot.org/uniprot/P42771
    # https://www.uniprot.org/uniprot/Q8N726
    primary_acc_ids = list(filter(lambda uid: uid in primary.values(), uids))
    if not primary_acc_ids:
        # the list that Uniprot provides is sec --> primary map
        # if an id is not among the primary keys it may mean (and seems to work)
        # that there are no secondary ids that map to this one
        if len(uids) == 1:
            primary_acc_ids = uids.copy()
        else:
            print(f"primary key not found uniprot_ids: {uniprot_ids}", file=stderr)
            exit(1)
    for paid in primary_acc_ids: other_ids[paid] = []
    secondary_acc_ids = list(filter(lambda uid: uid in primary.keys(), uids))
    for said in secondary_acc_ids:
        other_ids[primary[said]].append(said)
    return primary_acc_ids, other_ids


def find_pdb_paths_in_swissmodel_dir(sm_home, uniprot_id):
    uni_sub_path = "/".join(uniprot_id[i:i+2] for i in range(0, len(uniprot_id), 2))
    sm_path = f"{sm_home}/{uni_sub_path}/swissmodel"
    if not os.path.exists(sm_path): return []
    return [f"{uni_sub_path}/swissmodel/{pdb}" for pdb in os.listdir(sm_path)]

#########################################
def main():
    sm_home = "/storage/databases/swissmodel/human/"
    with open("top_genes.tsv") as inf:
        gene_names = [line.strip().split("\t")[0] for line in inf]

    db = connect_to_mysql(Config.mysql_conf_file)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species(cursor)
    switch_to_db(cursor, ensembl_db_name['homo_sapiens'])
    primary = read_primary_accession_numbers()

    print("\t".join(["gene", "popularity rank", "gene_id", "canonical_transcript",
                     "mulitple _primary_acc_ids", "uniprot_primary_acc_ids",
                     "some_paths_use_non-primary_acc_id", "structure_paths"]))
    for gene in gene_names:
        if gene == "non-coding": continue
        out = [gene, f"{gene_names.index(gene)+1}"]

        qry  = f"select  ensembl_id_by_ensembl, ensembl_gene_id, uniprot_id_by_uniprot "
        qry += f"from identifier_maps.hgnc where approved_symbol='{gene}'"
        ret = error_intolerant_search(cursor, qry)
        if not ret:
            out.append("gene entry not found in db")
        elif ret[0][0] == "" and ret[0][1] == "":
            out.append("ensembl transcript id not found")
        else:
            [ensembl_id_by_ensembl, ensembl_gene_id_by_hgnc, uniprot_ids] = ret[0]
            ensembl_gene_id = ensembl_id_by_ensembl if ensembl_id_by_ensembl != "" else ensembl_gene_id_by_hgnc
            [primary_acc_ids, other_ids] = check_primary(primary, uniprot_ids)
            out.append(ensembl_gene_id)

            qry  = "select transcript.stable_id  from transcript, gene "
            qry += "where gene.canonical_transcript_id = transcript.transcript_id "
            qry += f"and gene.stable_id = '{ensembl_gene_id}' "
            ensembl_transcript_id = hard_landing_search(cursor, qry)[0][0]
            out.extend([ensembl_transcript_id,  str(len(primary_acc_ids) > 1) , ", ".join(primary_acc_ids)])

            pdb_paths = []
            some_paths_non_primary = False
            for primary_acc_id in primary_acc_ids:  # see the comments in check_primary function
                pdb_paths.extend(find_pdb_paths_in_swissmodel_dir(sm_home, primary_acc_id))
                for uniprot_id in other_ids[primary_acc_id]:
                    alt_uniprot_id_paths = find_pdb_paths_in_swissmodel_dir(sm_home, uniprot_id)
                    if not alt_uniprot_id_paths: continue
                    some_paths_non_primary = True
                    pdb_paths.extend(alt_uniprot_id_paths)
            out.extend([str(some_paths_non_primary)] + pdb_paths)
        print("\t".join(out))

    cursor.close()

#########################################
if __name__ == '__main__':
    main()
