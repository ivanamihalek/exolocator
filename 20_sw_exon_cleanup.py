#!/usr/bin/python


import MySQLdb
from   el_utils.mysql   import  connect_to_mysql, search_db, switch_to_db
from   el_utils.ensembl import  *
# BioPython
from Bio.Seq      import Seq

#########################################
def main():
    
    no_threads = 1
    local_db   = False

    if local_db:
        db     = connect_to_mysql()
    else:
        db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)

    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)

    tot    = 0
    tot_ok = 0
    for species in all_species:

        if species=='homo_sapiens':
            continue

        switch_to_db(cursor, ensembl_db_name[species])
        qry      = "select * from sw_exon"
        sw_exons = search_db(cursor, qry)

        if not sw_exons:
            print "no sw exons for ", species
            continue

        ct = 0
        ok = 0
        for sw_exon in sw_exons:
            ct += 1
            
            has_stop = False
            has_NNN  = False
            [sw_exon_id, gene_id, start_in_gene, end_in_gene, human_exon_id,
             exon_seq_id, strand, phase, has_NNN, has_stop, has_3p_ss, has_5p_ss] = sw_exon
            [exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, 
             left_flank, right_flank, dna_seq] = get_exon_seq_by_db_id (cursor, exon_seq_id, ensembl_db_name[species])


            ############################
            # check for NNN stuff
            if 0:
                human_coding = is_coding (cursor, human_exon_id, ensembl_db_name['homo_sapiens'])
                if human_coding and protein_seq and '*' in protein_seq:
                    has_stop = True
                if 'NNN' in dna_seq:
                    has_NNN  = True
                switch_to_db(cursor, ensembl_db_name[species])
                qry  = "update sw_exon set has_stop = %d,  " % (1 if has_stop else 0)
                qry += "has_NNN = %d " % (1 if has_NNN else 0)
                qry += " where exon_id = %d " %  int(sw_exon_id) 
                rows = search_db (cursor, qry)
                if  rows:
                    print species
                    print qry
                    print rows
                    exit (1)
                if not has_stop and not has_NNN:
                    ok  += 1

            ########################
            # see if the  protein seq matches the quoted boundaries
            # coding dna sequence:
            cds = dna_seq[pepseq_transl_start:pepseq_transl_end]
            translated_cds = Seq(cds).translate().tostring()
            if not  translated_cds == protein_seq:
                print Seq(cds).translate().tostring()
                print protein_seq
                exit(1)

            #print exon_seq_id, human_coding, protein_seq

        tot_ok += ok
        tot += ct

        print species, ct, ok
        #exit(1)
    
    print tot, tot_ok

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
