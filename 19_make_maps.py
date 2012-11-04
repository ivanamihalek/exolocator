#!/usr/bin/python


import MySQLdb, commands
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  get_species, get_gene_ids, gene2stable, stable2gene
from   el_utils.ensembl import  get_compara_name, get_description, gene2exon_list
from   el_utils.ensembl import  get_exon_pepseq, genome_db_id2species
from   el_utils.utils   import  erropen
from   el_utils.map     import  Map
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator

#########################################
def  get_orthos (cursor, gene_id):

    orthos = []
    qry  = "select cognate_gene_id, cognate_genome_db_id from orthologue "
    qry += "where gene_id = %d"% gene_id
    
    rows = search_db (cursor, qry)
    if (not rows):
        return []

    for row in rows:
        ortho_gene_id   = row[0]
        ortho_genome_db = row[1]
        # note: the cursor will be pointing to compara db after this
        species = genome_db_id2species (cursor, ortho_genome_db)
        orthos.append([ortho_gene_id,species] )
    
    return orthos

#########################################
def  get_unresolved_orthos (cursor, gene_id):

    orthos = []
    qry  = "select cognate_gene_id, cognate_genome_db_id from unresolved_ortho "
    qry += "where gene_id = %d"% gene_id
    
    rows = search_db (cursor, qry)
    if (not rows):
        return []

    for row in rows:
        ortho_gene_id   = row[0]
        ortho_genome_db = row[1]
        # note: the cursor will be pointing to compara db after this
        species = genome_db_id2species (cursor, ortho_genome_db)
        orthos.append([ortho_gene_id,species] )
    
    return orthos

#########################################
def decorate_and_concatenate (exon_seqs_protein, hook):

    decorated_seq = hook
    for  pepseq in exon_seqs_protein:

        # does the peptide already contain the hook:
        if (hook in pepseq):
            print ">>>>>>>>>>>!!!!!!!!"
            print " hook  already  present "
            print pepseq
            exit(1)

        decorated_seq += pepseq+hook
    
    return decorated_seq

#########################################
def make_map (cursor, ensembl_db_name, acg, ortho_species, human_exons, ortho_exons):

    map = Map()

    print 'making map for', ortho_species

    hook = "WWWPWWW"

    # 1) choose exons that I need
    # 2) sort them by start position in the gene

    # switch database to human
    switch_to_db(cursor, ensembl_db_name['homo_sapiens'])
    protein_seq = []
    for exon in human_exons:
        pepseq = get_exon_pepseq (cursor, exon.exon_id)
        if not pepseq:
            continue
        protein_seq.append(pepseq)
    human_seq = decorate_and_concatenate (protein_seq, hook)

    # switch database to the ortho species
    switch_to_db(cursor, ensembl_db_name[ortho_species])
    protein_seq = []
    for exon in ortho_exons:
        pepseq = get_exon_pepseq (cursor, exon.exon_id)
        if not pepseq:
            continue
        protein_seq.append(pepseq)
        
    ortho_seq  = decorate_and_concatenate (protein_seq, hook)
 
    fastafile  = "{0}.fa".format (exon.exon_id)
    afafile    = "{0}.afa".format (exon.exon_id)
    
    outf = erropen (fastafile, "w")
    print >> outf, ">"+'homo_sapiens'
    print >> outf, human_seq
    print >> outf, ">"+ortho_species
    print >> outf, ortho_seq
    outf.close()
   
    mafftcmd = acg.generate_mafft_command (fastafile, afafile)
    ret      = commands.getoutput(mafftcmd)
 
    print afafile

    exit (1)
    return map
    
#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()
    acg    = AlignmentCommandGenerator()

    [all_species, ensembl_db_name] = get_species (cursor)

    species='homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])

    if 0:
        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')
    else:
        #gene_ids = [stable2gene(cursor,'ENSG00000116729')] # wls
        gene_ids = [stable2gene(cursor,'ENSG00000156970')] # BUB1B

    for gene_id in gene_ids:

        print get_description (cursor, gene_id)
        
        # get _all_ exons
        switch_to_db (cursor, ensembl_db_name[species])
        human_exons = gene2exon_list(cursor, gene_id)
        if (not human_exons):
            print 'no exons for ', gene_id
            sys.exit(1)

        switch_to_db (cursor, ensembl_db_name[species])
        known_orthologues = get_orthos (cursor, gene_id)
        for [ortho_gene_id, ortho_species] in known_orthologues:
            print ortho_gene_id, ortho_species
            ortho_exons = gene2exon_list(cursor, ortho_gene_id, db_name= ensembl_db_name[ortho_species] )
            if not ortho_exons:
                print 'no exons for ', species, ortho_gene_id
            map =  make_map (cursor, ensembl_db_name, acg, ortho_species, human_exons, ortho_exons)

        print "==============================="

        switch_to_db (cursor, ensembl_db_name[species])
        known_orthologues = get_unresolved_orthos (cursor, gene_id)
        for ortho in known_orthologues:
            print ortho

        exit(1)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
