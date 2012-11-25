#!/usr/bin/python
# make the best alignment we can using the maps
# we currently have at hand

import MySQLdb, commands, re
from   el_utils.mysql   import  connect_to_mysql, connect_to_db
from   el_utils.mysql   import  switch_to_db,  search_db, store_or_update
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  get_maps, Map
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader


from Bio       import SeqIO
from bitstring import Bits


#########################################
def main():

    local_db = False

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
    no_pepseq      = 0
    no_orthologues = 0
    #for gene_id in gene_ids:
    # for gene_id in [412667]: #  wls
    #for gene_id in [378768]: #  p53
    for gene_id in [378766]: #  dynein
        gene_ct += 1
        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        #print gene_ct, len(gene_ids),  gene2stable(cursor, gene_id), get_description (cursor, gene_id)

        # find all exons we are tracking in the database
        human_exons = gene2exon_list(cursor, gene_id)
        for human_exon in human_exons:
           
            tot += 1
            # find all orthologous exons the human exon  maps to
            maps = get_maps(cursor, ensembl_db_name, human_exon.exon_id, human_exon.is_known)
            if not maps: continue
            # output to fasta:
            seqname   = "{0}:{1}:{2}".format('homo_sapiens', human_exon.exon_id, human_exon.is_known)
            switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
            [exon_seq_id, pepseq, left_flank, right_flank, dna_seq]  \
                = get_exon_seqs (cursor, human_exon.exon_id, human_exon.is_known)
            if (not pepseq):
                if ( human_exon.is_coding and  human_exon.covering_exon <0): # this should be a master exon
                    print "no pep seq for",  human_exon.exon_id, "coding ", human_exon.is_coding,
                    print "canonical: ",  human_exon.is_canonical
                    print "length of dna ", len(dna_seq)
                no_pepseq += 1
                continue
            headers   = [seqname]
            sequences = {seqname:pepseq}
            for map in maps:
                switch_to_db (cursor, ensembl_db_name[map.species_2])
                pepseq  = get_exon_pepseq (cursor, map.exon_id_2, map.exon_known_2)
                if (not pepseq):
                    continue
                seqname = "{0}:{1}:{2}".format(map.species_2, map.exon_id_2, map.exon_known_2)
                headers.append(seqname)
                sequences[seqname] = pepseq

            fasta_fnm = "{0}/{1}.fa".format( cfg.dir_path['scratch'], human_exon.exon_id)
            output_fasta (fasta_fnm, headers, sequences)
            
            if (len(headers) <=1 ):
                #print "single species in the alignment"
                no_orthologues += 1
                continue

            # align
            afa_fnm  = "{0}/{1}.afa".format( cfg.dir_path['scratch'], human_exon.exon_id)
            mafftcmd = acg.generate_mafft_command (fasta_fnm, afa_fnm)
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
                if (cognate_species == 'homo_sapiens'):
                    human_seq_seen = True
                cognate_genome_db_id = species2genome_db_id(cursor, cognate_species) # moves the cursor
                switch_to_db(cursor, ensembl_db_name['homo_sapiens']) # so move it back to homo sapiens
                # Write the bitmap to the database
                #if (cognate_species == 'homo_sapiens'):
                store_or_update(cursor, "exon_map",    {"cognate_genome_db_id":cognate_genome_db_id,
                   "cognate_exon_id":cognate_exon_id   ,"cognate_exon_known"  :cognate_exon_known,
                   "exon_id"        :human_exon.exon_id,"exon_known"          :human_exon.is_known},
                  {"msa_bitstring":MySQLdb.escape_string(msa_bitmap)})

            ok += 1
            commands.getoutput("rm "+afa_fnm+" "+fasta_fnm)
    
    print "tot: ", tot, "ok: ", ok
    print "no pepseq ", no_pepseq
    print "no_orthologues  ", no_orthologues

#########################################
if __name__ == '__main__':
    main()
