#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re, sys
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  get_maps, map2exon, Map
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial

from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.special_gene_sets  import *
from   el_utils.processes import  parallelize, get_process_id


from time      import  time
from Bio       import  SeqIO
from bitstring import  Bits

verbose = False

#########################################
def multiple_exon_alnmt(gene_list, db_info):


    print "process pid: %d, length of gene list:%d" % ( get_process_id(), len(gene_list))

    [local_db, ensembl_db_name] = db_info


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
    

    species  = 'homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)

    # for each human gene
    gene_ct = 0
    tot  = 0
    ok   = 0
    no_maps        = 0
    no_pepseq      = 0
    no_orthologues = 0
    min_similarity = cfg.get_value('min_accptbl_exon_sim')

    for gene_id in gene_list:

        start = time()
        gene_ct += 1
        if  not gene_ct%10: print gene_ct, "genes out of", len(gene_list)

        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        print gene_ct, len(gene_ids),  gene_id,  gene2stable(cursor, gene_id), get_description (cursor, gene_id)

        human_exons = filter (lambda e: e.is_known==1 and e.is_coding and e.covering_exon<0, gene2exon_list(cursor, gene_id))
        human_exons.sort(key=lambda exon: exon.start_in_gene)

        headers = []
        sequences = {}
        for human_exon in human_exons:
            
            tot += 1
            # find all orthologous exons the human exon  maps to
            maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
            if verbose: 
                print "\texon no.", tot, " id", human_exon.exon_id,
                if not maps: 
                    print " no maps"
                    print human_exon
                print 
            if not maps: 
                no_maps += 1
                continue

            # output to fasta:
            seqname   = "{0}:{1}:{2}".format('homo_sapiens', human_exon.exon_id, human_exon.is_known)
            switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
            [exon_seq_id, pepseq, pepseq_transl_start, pepseq_transl_end, 
             left_flank, right_flank, dna_seq] = get_exon_seqs (cursor, human_exon.exon_id, human_exon.is_known)
            
            # human seqeunce
            if (not pepseq):
                if verbose and  human_exon.is_coding and  human_exon.covering_exon <0: # this should be a master exon
                    print "no pep seq for",  human_exon.exon_id, "coding ", human_exon.is_coding,
                    print "canonical: ",  human_exon.is_canonical
                    print "length of dna ", len(dna_seq)
                no_pepseq += 1
                continue

            hassw = False
            for map in maps:

                switch_to_db (cursor, ensembl_db_name[map.species_2])
                if map.similarity < min_similarity: continue
                exon    = map2exon(cursor, ensembl_db_name, map)
                pepseq  = get_exon_pepseq (cursor,exon)
                if (not pepseq):
                    continue
                if  map.source == 'sw_sharp':
                    exon_known_code = 2
                    hassw = True
                elif  map.source == 'usearch':
                    exon_known_code = 3
                    hassw = True
                else:
                    exon_known_code = map.exon_known_2
                seqname = "{0}:{1}:{2}".format(map.species_2, map.exon_id_2, exon_known_code)
                headers.append(seqname)
                sequences[seqname] = pepseq

 
            if (len(headers) <=1 ):
                if verbose: print "single species in the alignment"
                no_orthologues += 1
                continue
            
            fasta_fnm = "{0}/{1}.fa".format( cfg.dir_path['scratch'], human_exon.exon_id)
            output_fasta (fasta_fnm, headers, sequences)

            # align
            afa_fnm  = "{0}/{1}.afa".format( cfg.dir_path['scratch'], human_exon.exon_id)
            mafftcmd = acg.generate_mafft_command (fasta_fnm, afa_fnm)

            #mafftcmd = "muscle -in" + fasta_fnm + " -out" + afa_fnm
            ret      = commands.getoutput(mafftcmd)


            # read in the alignment
            inf = erropen(afa_fnm, "r")
            human_seq_seen = False
            for record in SeqIO.parse(inf, "fasta"):
                ### store the alignment as bitstring
                # Generate the bitmap
                bs         = Bits(bin='0b' + re.sub("[^0]","1", str(record.seq).replace('-','0')))
                msa_bitmap = bs.tobytes()
                
                # Retrieve information on the cognate
                cognate_species, cognate_exon_id, cognate_exon_known = record.id.split(':')
                if cognate_exon_known == '2':
                    source = 'sw_sharp'
                elif cognate_exon_known == '3':
                    source = 'usearch'
                else:
                    source = 'ensembl'
                if (cognate_species == 'homo_sapiens'):
                    human_seq_seen = True
                cognate_genome_db_id = species2genome_db_id(cursor, cognate_species) # moves the cursor
                switch_to_db(cursor, ensembl_db_name['homo_sapiens']) # so move it back to homo sapiens
                # Write the bitmap to the database
                #if (cognate_species == 'homo_sapiens'):
                if verbose and (source=='sw_sharp' or source=='usearch'):
                    print "storing"
                    print human_exon.exon_id, human_exon.is_known
                    print cognate_species, cognate_genome_db_id, cognate_exon_id, cognate_exon_known, source
                    #print MySQLdb.escape_string(msa_bitmap)
                    if not msa_bitmap:
                        print "no msa_bitmap"
                        continue
                    print
                store_or_update(cursor, "exon_map",    {"cognate_genome_db_id":cognate_genome_db_id,
                   "cognate_exon_id":cognate_exon_id   ,"cognate_exon_known"  :cognate_exon_known,
                   "source": source, "exon_id" :human_exon.exon_id, "exon_known":human_exon.is_known},
                  {"msa_bitstring":MySQLdb.escape_string(msa_bitmap)})

            ok += 1
            commands.getoutput("rm "+afa_fnm+" "+fasta_fnm)

        if verbose: print " time: %8.3f\n" % (time()-start);

    print "tot: ", tot, "ok: ", ok
    print "no maps ",   no_pepseq
    print "no pepseq ", no_pepseq
    print "no orthologues  ", no_orthologues
    print


#########################################
def main():
    
    no_threads = 12
    special    = None


    if len(sys.argv) > 1 and  len(sys.argv)<3:
        print "usage: %s <set name> <from> <to>" % sys.argv[0]
        exit(1) # after usage statement
    elif len(sys.argv)>=3:

        special = sys.argv[1]
        special = special.lower()
        if special == 'none': special = None

        genes_from = int(sys.argv[2])
        genes_to   = ""
        if len(sys.argv)>=4:
            genes_to   = int(sys.argv[3])

    local_db = False

    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
    else:
        db  = connect_to_mysql    (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)

    print "running ", sys.argv[0]

    if special:
        print "using", special, "set"
        if special == 'complement':
            gene_list = get_complement_ids(cursor, ensembl_db_name, cfg)
        else:
            gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special )
        multiple_exon_alnmt, gene_list, [local_db, ensembl_db_name])

    else:
        print "using all protein coding genes"
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        gene_list = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        if genes_to:
            multiple_exon_alnmt ( gene_list[genes_from, genes_to], [local_db, ensembl_db_name]) 
        else:
            multiple_exon_alnmt ( gene_list[genes_from:], [local_db, ensembl_db_name]) 
            
 

    cursor.close()
    db.close()

    #parallelize (no_threads, multiple_exon_alnmt, gene_list, [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
    main()

'''
    #for gene_id in [412667]: #  wls
    #for gene_id in [378768]: #  p53
    #for gene_id in [378766]: #  dynein
'''

