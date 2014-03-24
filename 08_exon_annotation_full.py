#!/usr/bin/python

import MySQLdb
import commands, sys
from   random           import choice
from   el_utils.mysql   import connect_to_mysql, search_db, switch_to_db, check_null
from   el_utils.mysql   import store_or_update
from   el_utils.ensembl import  *
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.translation import crop_dna, translation_bounds, translate
from   el_utils.threads     import parallelize
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.special_gene_sets  import get_theme_ids
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

#########################################
def get_phase(cursor, exon_id):

    qry  = "select is_coding, phase, gene_id from gene2exon where exon_id = %d" % exon_id
    rows = search_db(cursor, qry)
    if (rows):
        [is_coding, phase, gene_id] = rows[0]
    else:
        [is_coding, phase, gene_id] = [0,0,0]

    return [is_coding, phase, gene_id]


#########################################
def phase2offset(phase):
    if phase > 2:
        phase = phase%3
    if phase==0:
        offset = 0
    else:
        offset = 3-phase
    return offset

########################################
def pep_exon_seqs(species_list, db_info):

    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql          (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    #####################################
    for species in species_list:

        print
        print "############################"
        print  species
        sys.stdout.flush()

        if not switch_to_db(cursor, ensembl_db_name[species]):
            return False

        if (species=='homo_sapiens'):
            gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
        else:
            gene_ids = get_gene_ids (cursor, biotype='protein_coding')

        tot         = 0
        gene_ct     = 0
        no_pepseq   = 0
        pepseq_ok   = 0
        short_dna   = 0
        no_exon_seq = 0
        translation_fail = 0

        #for all protein coding genes in a species
        for gene_id in gene_ids:

            # for all exons in the gene
            exons = gene2exon_list(cursor, gene_id)
            if (not exons):
                print 'no exons for gene', gene_id
                continue

            for exon in exons:

                #####################################                
                if (not exon.is_coding or  exon.covering_exon > 0):
                    continue 
                tot += 1
                exon_seqs = get_exon_seqs(cursor, exon.exon_id, exon.is_known)
                if (not exon_seqs):
                    no_exon_seq += 1
                    continue                   
                [exon_seq_id, pepseq, pepseq_transl_start, 
                 pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs
                if len(dna_seq)<3:
                    short_dna += 1
                    continue

                #####################################                
                mitochondrial        = is_mitochondrial(cursor, gene_id)
                [seq_start, seq_end] = translation_bounds (cursor, exon.exon_id)
                dna_cropped          = crop_dna (seq_start, seq_end, dna_seq)
                [offset, pepseq]     = translate (dna_cropped, exon.phase, mitochondrial)

                if ( not pepseq): # usually some short pieces (end in pos 4 and such)
                    translation_fail += 1
                    continue
 
                if seq_start is None: seq_start = 1
                if seq_start == 0: seq_start = 1
                start = seq_start+offset-1
                end   = start + 3*len(pepseq)
                
                dnaseq  = Seq (dna_seq[start:end], generic_dna)
                if (mitochondrial):
                    pepseq2 = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
                else:
                    pepseq2 = dnaseq.translate().tostring()

                if (not pepseq == pepseq2):
                    start = -10
                    end   = -10
                    translation_fail += 1
                else:
                    pepseq_ok += 1

                qry  = "update exon_seq "
                qry += "set protein_seq   = '%s',  "  %  pepseq
                qry += " pepseq_transl_start =  %d, " %  start
                qry += " pepseq_transl_end   =  %d  " %  end
                qry += " where exon_seq_id =  %d  "   %  exon_seq_id
                rows = search_db (cursor, qry)
                if (rows):
                    rows = search_db (cursor, qry, verbose = True)
                    continue
            ####################################
            if not gene_ids.index(gene_id)%1000:
                print "%50s:  %5.1f%% " %  (species, 100*(float( gene_ids.index(gene_id) +1 )/len(gene_ids))  )
                sys.stdout.flush()
       
                 
        print species, "done"
        print "total coding exons ", tot
        print "no exon seq info   ", no_exon_seq
        print "short dna          ", short_dna
        print "transl failure     ", translation_fail
        print "pepseq ok          ", pepseq_ok
        sys.stdout.flush()


########################################
def pep_exon_seqs_special (gene_list, db_info):

    [local_db, ensembl_db_name] = db_info
    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql          (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    tot         = 0
    gene_ct     = 0
    no_pepseq   = 0
    pepseq_ok   = 0
    short_dna   = 0
    no_exon_seq = 0
    translation_fail = 0

    #####################################
    for gene_id in gene_list:

        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        orthologues  = get_orthos (cursor, gene_id, 'orthologue') # get_orthos changes the db pointer

        switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
        orthologues += get_orthos (cursor, gene_id, 'unresolved_ortho')

        for [ortho_gene_id, ortho_species] in [[gene_id,'homo_sapiens']] + orthologues:
 
            print ">>> ", ortho_species, ortho_gene_id
            switch_to_db (cursor, ensembl_db_name[ortho_species])

            # for all exons in the gene
            exons = gene2exon_list(cursor, ortho_gene_id)
            if (not exons):
                print 'no exons for gene', ortho_gene_id
                continue

            for exon in exons:

                #####################################                
                if (not exon.is_coding or  exon.covering_exon > 0):
                    continue 
                tot += 1
                exon_seqs = get_exon_seqs(cursor, exon.exon_id, exon.is_known)
                if (not exon_seqs):
                    no_exon_seq += 1
                    continue                   
                [exon_seq_id, pepseq, pepseq_transl_start, 
                 pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs
                if len(dna_seq)<3:
                    short_dna += 1
                    continue

                #####################################                
                mitochondrial        = is_mitochondrial(cursor, ortho_gene_id)
                [seq_start, seq_end] = translation_bounds (cursor, exon.exon_id)
                dna_cropped          = crop_dna (seq_start, seq_end, dna_seq)
                [offset, pepseq]     = translate (dna_cropped, exon.phase, mitochondrial)

                if ( not pepseq): # usually some short pieces (end in pos 4 and such)
                    translation_fail += 1
                    continue
 
                if seq_start is None: seq_start = 1
                if seq_start == 0: seq_start = 1
                start = seq_start+offset-1
                end   = start + 3*len(pepseq)
                
                dnaseq  = Seq (dna_seq[start:end], generic_dna)
                if (mitochondrial):
                    pepseq2 = dnaseq.translate(table="Vertebrate Mitochondrial").tostring()
                else:
                    pepseq2 = dnaseq.translate().tostring()

                if (not pepseq == pepseq2):
                    start = -10
                    end   = -10
                    translation_fail += 1
                else:
                    pepseq_ok += 1

                qry  = "update exon_seq "
                qry += "set protein_seq   = '%s',  "  %  pepseq
                qry += " pepseq_transl_start =  %d, " %  start
                qry += " pepseq_transl_end   =  %d  " %  end
                qry += " where exon_seq_id =  %d  "   %  exon_seq_id
                rows = search_db (cursor, qry)
                if (rows):
                    rows = search_db (cursor, qry, verbose = True)
                    continue

        ####################################
        if not gene_list.index(gene_id)%1000:
            print "%5.1f%% " %  (100*(float( gene_list.index(gene_id) +1 )/len(gene_list))  )
            sys.stdout.flush()
       
                 

#########################################
def main():

    no_threads = 1
    special    = ''
    local_db = False

    if len(sys.argv) > 1 and  len(sys.argv)<3  or len(sys.argv) >= 2 and sys.argv[1]=="-h":
        print "usage: %s <set name> <number of threads>" % sys.argv[0]
        exit(1) # after usage statment
    elif len(sys.argv)==3:
        special = sys.argv[1].lower()
        if special == 'none': special = None
        no_threads = int(sys.argv[2])


    if local_db:
        db     = connect_to_mysql()
        cfg    = ConfigurationReader()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg    = ConfigurationReader       (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)
    print '======================================='
    print sys.argv[0]
    if special:
        print "using", special, "set"
        gene_list = get_theme_ids (cursor,  ensembl_db_name, cfg, special )
 
    cursor.close()
    db    .close()

    if not special:
        parallelize (no_threads, pep_exon_seqs, all_species, [local_db, ensembl_db_name] )
    else:
        parallelize (no_threads, pep_exon_seqs_special, gene_list,  [local_db, ensembl_db_name])



#########################################
if __name__ == '__main__':
    main()
