#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  get_species, get_gene_ids, gene2stable, stable2gene
from   el_utils.ensembl import  get_compara_name, get_description, gene2exon_list
from   el_utils.ensembl import  get_exon_pepseq, genome_db_id2species, species2taxid
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  Map
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader

#########################################
def get_maps(cursor, ensembl_db_name, exon_id, is_known):
    
    maps = []

    switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
    qry  = "select * from exon_map where exon_id = %d " % exon_id
    qry += " and exon_known = %d " % is_known
    rows = search_db (cursor, qry)
    if not rows or "ERROR" in rows[0]:
        return []

    for row in rows:
        map = Map()
        map.load_from_db(row, cursor)
        maps.append(map)

    return maps



#########################################
def main():

    db     = connect_to_mysql()
    cursor = db.cursor()
    cfg    = ConfigurationReader()
    acg    = AlignmentCommandGenerator()
    # find db ids adn common names for each species db
    [all_species, ensembl_db_name] = get_species (cursor)
    

    species  = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    # for each human gene
    for gene_id in gene_ids:
        print gene_id, gene2stable(cursor, gene_id), get_description (cursor, gene_id)

        # find all exons we are tracking in the database
        human_exons = gene2exon_list(cursor, gene_id)
        for human_exon in human_exons:
            print "\t", human_exon.exon_id
            # find all orthologous exons the human exon  maps to
            maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)

            # output to fasta:
            seqname   = "{0}:{1}:{2}".format('homo_sapiens', human_exon.exon_id, human_exon.is_known)
            switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
            pepseq    = get_exon_pepseq (cursor, human_exon.exon_id, human_exon.is_known)
            headers   = [seqname]
            sequences = [pepseq]
            for map in maps:
                switch_to_db (cursor, ensembl_db_name[map.species_2])
                pepseq  = get_exon_pepseq (cursor, map.exon_id_2, map.exon_known_2)
                if (not pepseq):
                    continue
                seqname = "{0}:{1}:{2}".format(map.species_2, map.exon_id_2, map.exon_known_2)
                headers.append(seqname)
                sequences.append(pepseq)
                print map
                print "pepseq: ", pepseq
            fasta_fnm = "{0}/{1}.fa".format( cfg.dir_path['scratch'], human_exon.exon_id)
            output_fasta (fasta_fnm, headers, sequences)

            # align
            afa_fnm = "{0}/{1}.afa".format( cfg.dir_path['scratch'], human_exon.exon_id)
            mafftcmd = acg.generate_mafft_command (fasta_fnm, afa_fnm)
            ret      = commands.getoutput(mafftcmd)
            
            print fasta_fnm, afa_fnm
            # read in the alignment
            inf = erropen(afa_fnm, "r")
            for line in inf:
                if ('>' in line):
                    print line,
            # store the alignment as bitstring <<<<<<<<<<<<<<<<<
            # note the store_or_update function
            exit (1)



#########################################
if __name__ == '__main__':
    main()
