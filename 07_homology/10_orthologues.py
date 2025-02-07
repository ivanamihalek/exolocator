#!/usr/bin/python3 -u
# we are collecting orthologues for human only
# the possibility that "A is ortholog of B, and B is ortholog of C, but A is not ortholog of C"
# is something we don't want to know about right now


from time import time

from dotenv import load_dotenv

from config import Config
from el_utils.ensembl import *
from el_utils.processes import *

#########################################
# ./el_utils/kernprof.py -l <calling script>.py
# python3 -m line_profiler <calling script>.py.lprof
# @profile
############
def get_orthologues(cursor, compara_db, homology_ids, qry_stable_id, verbose=False):
    orthos = {}
    # multiple threads are accessing the same table here, this is just now qorking particularly well
    for homid in homology_ids:
        qry = f"select description from {compara_db}.homology_human where homology_id={homid}"
        homology_description = hard_landing_search(cursor, qry)[0][0]

        qry = f"select gene_member_id from {compara_db}.homology_member_human where homology_id={homid}"
        ret = error_intolerant_search(cursor, qry)
        if not ret: continue  # strange, but possible

        gmis = ",".join([str(r[0]) for r in ret])
        qry = f"select stable_id, genome_db_id from {compara_db}.gene_member where gene_member_id in ({gmis})"
        ret = error_intolerant_search(cursor, qry)
        if not ret: continue

        for stable_id, genome_db_id in ret:
            if stable_id == qry_stable_id: continue
            qry = "select db_name from exolocator_meta.db_names where genome_db_id=%d" % genome_db_id
            ret = error_intolerant_search(cursor, qry)
            if not ret: continue  # this is genome we are not using - one of many mice strains or some such

            db_name = ret[0][0]

            qry = f"select gene_id from {db_name}.gene where stable_id='{stable_id}'"
            ret = error_intolerant_search(cursor, qry)
            if not ret: continue  # investigate at some point why this happens
            ortho_gene_id = ret[0][0]

            if not homology_description in orthos:  orthos[homology_description] = []
            orthos[homology_description].append([ortho_gene_id, genome_db_id])

    return orthos


#########################################
# writing directly to db way too slow
def write_orthologues(ortho_file_handle, gene_id, orthos):
    for ortho in orthos:
        [ortho_gene_id, ortho_genome_db_id] = ortho
        ortho_file_handle.write("{}\t{}\t{}\n".format(gene_id, ortho_gene_id, ortho_genome_db_id))
        ortho_file_handle.flush()
    return


def open_files():
    pid = os.getpid()
    filehandle = {'orthologue': open(f"raw_tables/orthologues.{pid}.tsv", "w"),
                  'unresolved_ortho': open(f"raw_tables/unresolved_orthos.{pid}.tsv", "w")}
    return filehandle


def collect_orthologues(cursor, ensembl_compara_name, ensembl_db_name, stable_id, homology_ids, filehandle):
    switch_to_db(cursor, ensembl_db_name['homo_sapiens'])
    search_db(cursor, "set autocommit=1")

    ortho_table = {'ortholog_one2one': 'orthologue', 'apparent_ortholog_one2one': 'orthologue',
                   'possible_ortholog': 'unresolved_ortho', 'ortholog_one2many': 'unresolved_ortho',
                   'ortholog_many2many': 'unresolved_ortho'}

    gene_id = stable2gene(cursor, stable_id, ensembl_db_name['homo_sapiens'])
    # in compara table, get everything that homology has to say about
    # the possible orthologues
    # find all orthologous pairs suggested for this gene
    # f get_orthologues(cursor, compara_db, homology_ids, gene_stable_id, verbose=False)
    orthos = get_orthologues(cursor, ensembl_compara_name, homology_ids, stable_id, verbose=False)

    for ortho_type in ['ortholog_one2one', 'possible_ortholog', 'apparent_ortholog_one2one',
                       'ortholog_one2many', 'ortholog_many2many']:
        if ortho_type not in orthos or not orthos[ortho_type]: continue
        # the triple returned for each ortho type is [ortho_stable, species,  genome_db_id]
        table = ortho_table[ortho_type]
        write_orthologues(filehandle[table], gene_id, orthos[ortho_type])


def core_loop(genes, other_args):
    [ensembl_compara_name, ensembl_db_name] = other_args

    cursor =  mysql_using_env_creds()

    filehandle = open_files()

    time1 = time0 = time()
    ct = 0
    for gene_member_id, stable_id in genes:
        qry = f"select homology_id from {ensembl_compara_name}.homology_member_human  where gene_member_id={gene_member_id}"
        ret = time_qry(cursor, qry, verbose=False)
        if not ret: continue
        homology_ids = [r[0] for r in ret]
        collect_orthologues(cursor, ensembl_compara_name, ensembl_db_name, stable_id, homology_ids, filehandle)
        ct += 1
        if not ct % 100:
            print(os.getpid(), ct, "  last 100: %.1f mins" % ((time() - time1) / 60))
            time1 = time()

    minutes = (time() - time0) / 60
    print(os.getpid(), "done in %.1f mins" % minutes)

    for fh in filehandle.values(): fh.close()

    mysql_server_conn_close(cursor)


def main():
    cursor = mysql_using_env_creds()

    ensembl_compara_name = get_compara_name(cursor)
    [all_species, ensembl_db_name] = get_species(cursor)

    qry = f"select gene_member_id, stable_id from {ensembl_compara_name}.gene_member "
    qry += "where taxon_id=9606 and biotype_group='coding'"
    genes = hard_landing_search(cursor, qry)
    mysql_server_conn_close(cursor)

    core_loop(genes, [ensembl_compara_name, ensembl_db_name])

    return True


#########################################
if __name__ == '__main__':
    main()

'''
none of this worked:
create  table  scratch_gene_mem_subtable  engine=MYISAM  as select gene_member_id, 
stable_id from gene_member where taxon_id=9606  and biotype_group='coding' ;

Query OK, 23471 rows affected (0.66 sec)
Records: 23471  Duplicates: 0  Warnings: 0

'''
