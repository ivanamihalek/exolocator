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
from   el_utils.threads import  parallelize
from   el_utils.custom  import  get_theme_ids
from   random           import  choice

from time      import  time
from Bio       import  SeqIO
from bitstring import  Bits



#########################################
def multiple_exon_alnmt(species_list, db_info):


    [local_db, ensembl_db_name] = db_info

    verbose  = False

    if local_db:
        db     = connect_to_mysql()
        cfg    = ConfigurationReader()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()


    for species in species_list:
        if species == 'homo_sapiens': continue
        print
        print "############################"
        print  species

        switch_to_db (cursor,  ensembl_db_name[species])
        gene_ids = get_gene_ids (cursor, biotype='protein_coding')
        #gene_ids = get_theme_ids(cursor, cfg, 'wnt_pathway')
        if not gene_ids:
            print "no gene_ids"
            continue


        gene_ct       = 0
        tot           = 0
        ok            = 0
        no_maps       = 0
        no_pepseq     = 0
        no_paralogues = 0
        #for gene_id in gene_ids:
        for gene_id in [378128]: #   
        #for sample_ct in range(10):
            #gene_id = choice(gene_ids)

            if verbose: start = time()
            gene_ct += 1
            if not gene_ct%100: print species, gene_ct, "genes out of", len(gene_ids)
            if verbose: 
                print
                print gene_id, gene2stable(cursor, gene_id), get_description (cursor, gene_id)

            # get the paralogues - only the representative for  the family will have this 
            paralogues = get_paras (cursor, gene_id)  
            if not paralogues:
                if verbose:  print "\t not a template or no paralogues"
                continue

            if verbose:  print "paralogues: ", paralogues

            # get _all_ exons
            template_exons = gene2exon_list(cursor, gene_id)
            if (not template_exons):
                if verbose: print 'no exons for ', gene_id
                continue

            # find all template  exons we are tracking in the database
            for template_exon in template_exons:

                if verbose: print template_exon.exon_id
                maps = get_maps(cursor, ensembl_db_name, template_exon.exon_id,
                                template_exon.is_known, species=species, table='para_exon_map')

                if not maps:
                    no_maps += 1
                    continue

                # output to fasta:
                seqname   = "{0}:{1}:{2}".format('template', template_exon.exon_id, template_exon.is_known)
                [exon_seq_id, pepseq, pepseq_transl_start, pepseq_transl_end, 
                 left_flank, right_flank, dna_seq] = get_exon_seqs (cursor, template_exon.exon_id, template_exon.is_known)
                if (not pepseq):
                    if ( template_exon.is_coding and  template_exon.covering_exon <0): # this should be a master exon
                        print "no pep seq for",  template_exon.exon_id, "coding ", template_exon.is_coding,
                        print "canonical: ",  template_exon.is_canonical
                        print "length of dna ", len(dna_seq)
                        no_pepseq += 1
                    continue
                
                tot += 1

                sequences = {seqname:pepseq}
                headers   = [seqname]
                for map in maps:
                    pepseq  = get_exon_pepseq (cursor, map.exon_id_2, map.exon_known_2)
                    if (not pepseq):
                        continue
                    seqname = "{0}:{1}:{2}".format('para', map.exon_id_2, map.exon_known_2)
                    headers.append(seqname)
                    sequences[seqname] = pepseq

                fasta_fnm = "{0}/{1}.fa".format( cfg.dir_path['scratch'], template_exon.exon_id)
                output_fasta (fasta_fnm, headers, sequences)

 
                if (len(headers) <=1 ):
                    print "single species in the alignment (?)"
                    no_paralogues += 1
                    continue

                # align
                afa_fnm  = "{0}/{1}.afa".format( cfg.dir_path['scratch'], template_exon.exon_id)
                mafftcmd = acg.generate_mafft_command (fasta_fnm, afa_fnm)
                ret      = commands.getoutput(mafftcmd)

                # read in the alignment
                inf = erropen(afa_fnm, "r")
                if not inf:
                    print gene_id
                    continue
                template_seq_seen = False
                for record in SeqIO.parse(inf, "fasta"):
                    ### store the alignment as bitstring
                    # Generate the bitmap
                    bs         = Bits(bin='0b' + re.sub("[^0]","1", str(record.seq).replace('-','0')))
                    msa_bitmap = bs.tobytes()
                    # Retrieve information on the cognate
                    label, cognate_exon_id, cognate_exon_known = record.id.split(':')
                    if (label == 'template'):
                        template_seq_seen = True
                    # Write the bitmap to the database
                    #print "updating: ", template_exon.exon_id
                    store_or_update(cursor, "para_exon_map", {"cognate_exon_id"    :cognate_exon_id,
                                                         "cognate_exon_known" :cognate_exon_known,
                                                         "exon_id"            :template_exon.exon_id,
                                                         "exon_known"         :template_exon.is_known},
                                    {"msa_bitstring":MySQLdb.escape_string(msa_bitmap)})
                inf.close()
                ok += 1
                commands.getoutput("rm "+afa_fnm+" "+fasta_fnm)
            if verbose: print " time: %8.3f\n" % (time()-start);
 
        print species, "tot: ", tot,  "ok: ", ok
        print "no maps       ", no_pepseq
        print "no pepseq     ", no_pepseq
        print "no paralogues ", no_paralogues
        print


#########################################
def main():
    no_threads = 1

    local_db = False

    if local_db:
        db = connect_to_mysql()
    else:
        db = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    cursor.close()
    db    .close()

    parallelize (no_threads, multiple_exon_alnmt, all_species, [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
    main()

'''
    #for gene_id in [412667]: #  wls
    #for gene_id in [378768]: #  p53
    #for gene_id in [378766]: #  dynein
'''

