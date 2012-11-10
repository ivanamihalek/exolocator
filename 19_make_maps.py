#!/usr/bin/python


import MySQLdb, commands, re
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
def decorate_and_concatenate (exon_seqs_protein):
    decorated_seq = ""
    count = 1
    for  pepseq in exon_seqs_protein:
        padded_count = "{0:03d}".format(count)
        decorated_seq += 'B'+padded_count+pepseq+'Z'
        count += 1

    return decorated_seq

#########################################
def moveB (seq):
    seqlist = list(seq)
    begin_pattern = re.compile("B\d{3}\-+")

    for begin_label in begin_pattern.finditer(seq):
        start =  begin_label.start()
        end   =  begin_label.end()
        label =  seq[start:end]
        for i in range(4):
            seqlist[start+i] = '-'
        for i in range(4):
            seqlist[end-4+i] = label[i]

    return "".join(seqlist)

#########################################
def find_exon_positions(seq):

    exon_position = {}
    
    exon_pattern = re.compile("B.*?Z")
    for match in exon_pattern.finditer(seq):
        start       = match.start()
        end         = match.end()
        exon_seq_no = seq[start+1:start+4]
        exon_position[exon_seq_no] = [start+4, end-2]  # python's non-inclusive convention, minus Z

    return  exon_position


#########################################
def  pad_the_alnmt (exon_seq_human, human_start, exon_seq_other, other_start):
    
    seq_human = ""
    seq_other = ""

    padding = ""
    if ( human_start > other_start):
        for i in range (human_start-other_start):
            padding += "-"
    seq_human = padding + exon_seq_human


    padding = ""
    if ( other_start > human_start):
        for i in range (other_start-human_start):
            padding += "-"
    seq_other = padding + exon_seq_other

    if ( len(seq_human) >  len(seq_other)):
        padding = ""
        for i in range  (len(seq_human)-len(seq_other)):
            padding += "-"
        seq_other += padding

    if ( len(seq_other) >  len(seq_human)):
        padding = ""
        for i in range  (len(seq_other)-len(seq_human)):
            padding += "-"
        seq_human += padding


    return [seq_human, seq_other] 

    
#########################################
def alignment_line (seq_human, seq_other):

    alignment_line = []

    if ( not len(seq_human) ==  len(seq_other) ):
        print "alignment_line:  the seqeunces must be aligned"
        exit(1)
    else:
        length = len(seq_human)

    for i in range(length):
        if not seq_human[i] == "-" and  not seq_other[i] == "-":
            alignment_line.append ("|")

        elif seq_human[i] == "-" and  seq_other[i] == "-":
            #pass
            alignment_line.append ("-")

        elif (seq_human[i]  == "-" ):
            alignment_line.append ("A")

        elif (seq_other[i]  == "-" ):
            alignment_line.append ("B")
    return  "".join(alignment_line)


#########################################
def cigar_line (seq_human, seq_other):

    cigar_line     = []

    alignment_line = []

    if ( not len(seq_human) ==  len(seq_other) ):
        print "alignment_line:  the seqeunces must be aligned"
        exit(1)
    else:
        length = len(seq_human)

    for i in range(length):
        if not seq_human[i] == "-" and  not seq_other[i] == "-":
            alignment_line.append ("M")

        elif seq_human[i] == "-" and  seq_other[i] == "-":
            pass
            #alignment_line.append ("-")

        elif (seq_human[i]  == "-" ):
            alignment_line.append ("A")

        elif (seq_other[i]  == "-" ):
            alignment_line.append ("B")

            
    prev_char = alignment_line[0]
    count     = 1
    for i in range(1,len(alignment_line)):
        if ( alignment_line[i] == prev_char):
            count += 1
        else:
            cigar_line.append( "{0}{1}".format(count, prev_char))
            prev_char = alignment_line[i]
            count     = 1
                               
    cigar_line.append("{0}{1}".format(count, prev_char))

    return  "".join(cigar_line)


#########################################
def unfold_cigar_line (seq_A, seq_B, cigar_line):

    seq_A_aligned = ""
    seq_B_aligned = ""


    char_pattern = re.compile("\D")
    a_ct = 0
    b_ct = 0
    prev_end  = 0

    for match in char_pattern.finditer(cigar_line):
        this_start       = match.start()
        no_repeats = int(cigar_line[prev_end:this_start])
        prev_end         = match.end()
        alignment_instruction = cigar_line[this_start]

        if alignment_instruction == 'M':
            seq_A_aligned += seq_A[a_ct:a_ct+no_repeats]
            a_ct  += no_repeats
            seq_B_aligned += seq_B[b_ct:b_ct+no_repeats]
            b_ct  += no_repeats

        elif alignment_instruction == 'A':
            seq_A_aligned += '-'*no_repeats 
            seq_B_aligned += seq_B[b_ct:b_ct+no_repeats]
            b_ct  += no_repeats
            
        elif alignment_instruction == 'B':
            seq_A_aligned += seq_A[a_ct:a_ct+no_repeats]
            a_ct  += no_repeats
            seq_B_aligned +=  '-'*no_repeats 
           

    return [seq_A_aligned, seq_B_aligned]


#########################################
def maps_evaluate (human_exons, ortho_exons, aligned_seq, exon_positions):

    maps = []

    if len(aligned_seq.keys()) > 2:
        print "right now the mapping implemented for two species only"
        exit (1)

    for species in aligned_seq.keys():
        if species == 'homo_sapiens': continue
        other_species = species
        break

    print "other_species ", other_species
    print "human ", len(human_exons), len( exon_positions['homo_sapiens'])
    print "other ", len(ortho_exons), len( exon_positions[other_species])

 
    for human_exon_ct in range(len(human_exons)):

        padded_count_human = "{0:03d}".format(human_exon_ct+1)
        if ( not  exon_positions['homo_sapiens'].has_key(padded_count_human) ):
            continue
        [human_start, human_end] = exon_positions['homo_sapiens'][padded_count_human]


        for ortho_exon_ct in range(len(ortho_exons)):

            padded_count_ortho = "{0:03d}".format(ortho_exon_ct+1)
            if ( not  exon_positions[other_species].has_key(padded_count_ortho) ):
                continue
            [other_start, other_end] = exon_positions[other_species][padded_count_ortho]

            if ( human_start <= other_start <= human_end or
                 human_start <= other_end <= human_end):
                
                map = Map()
                map.species_1 = 'homo_sapiens'
                map.species_2 = other_species
                
                map.exon_id_1 = human_exons[human_exon_ct].exon_id
                map.exon_id_2 = ortho_exons[ortho_exon_ct].exon_id

                print map.exon_id_1, " --> ", map.exon_id_2

                exon_seq_human = aligned_seq['homo_sapiens'][human_start:human_end]
                exon_seq_other = aligned_seq[other_species][other_start:other_end]
                [seq_human, seq_other] = pad_the_alnmt (exon_seq_human,human_start,
                                                        exon_seq_other, other_start)
                ciggy = cigar_line (seq_human, seq_other)
                print seq_human
                print alignment_line (seq_human, seq_other)
                print seq_other
                print ciggy
                 
                print 

                [seq1, seq2] = unfold_cigar_line (seq_human.replace('-',''), seq_other.replace('-',''), ciggy)
                print seq1
                print seq2
                print 
                print 

                #exit (1)

                map.cigar_line = ciggy
                #map.similarity = similarity (seq[], seq[])
                
                maps.append(map)
    #exit (1)

    return maps

#########################################
def find_relevant_exons (cursor, all_exons):

    relevant_exons = []
    protein_seq    = []

    # 1) choose exons that I need
    for exon in all_exons:
        if (not exon.is_coding or  exon.covering_exon > 0):
            continue
        relevant_exons.append(exon)

    # 2) sort them by their start position in the gene
    to_remove = []
    relevant_exons.sort(key=lambda exon: exon.start_in_gene)
    for i in range(len(relevant_exons)):
        exon   = relevant_exons[i]
        pepseq = get_exon_pepseq (cursor, exon.exon_id)
        if not pepseq:
            print "no seq found for ", exon.exon_id, exon.is_known
            to_remove.append(i)
            continue
        pepseq = pepseq.replace ('X', '')
        if  not pepseq:
            to_remove.append(i)
            continue
        #print exon.start_in_gene,  exon.end_in_gene, exon.covering_exon , pepseq
        print "*"+pepseq+"*"
        protein_seq.append(pepseq)

    # remove all exons for which the pepseq  is empty
    for i in range (0, len(to_remove), -1):
        del relevant_exons[to_remove[i]]
 
    print " ** ", to_remove
    print protein_seq
    print

    return [relevant_exons, protein_seq]

#########################################
def output_fasta (filename, header, sequence):
    outf = erropen (filename, "w")
    print >> outf, ">"+header
    print >> outf, sequence
    outf.close()
    return

#########################################
def make_maps (cursor, ensembl_db_name, acg, ortho_species, human_exons, ortho_exons):

    maps = []

    print 'making map for', ortho_species
    switch_to_db(cursor,  ensembl_db_name['homo_sapiens'])
    print '    human seq'
    [relevant_human_exons, exon_protein_seq] = find_relevant_exons (cursor, human_exons)
    # concatenate exons' seqs
    human_seq = decorate_and_concatenate (exon_protein_seq)

    #print "##############################"
    switch_to_db(cursor,  ensembl_db_name[ortho_species])
    print '    other seq'
    [relevant_ortho_exons, exon_protein_seq] = find_relevant_exons (cursor, ortho_exons)
    # concatenate exons' seqs
    ortho_seq  = decorate_and_concatenate (exon_protein_seq)
 
    # output
    fastafile1  = "{0}.1.fa".format (relevant_human_exons[0].exon_id)
    output_fasta (fastafile1, 'homo_sapiens', human_seq)
    # 
    fastafile2  = "{0}.2.fa".format (relevant_ortho_exons[0].exon_id)
    output_fasta (fastafile2, ortho_species, ortho_seq)

    # align
    swsharpcmd = acg.generate_SW_peptide (fastafile1, fastafile2)
    ret        = commands.getoutput(swsharpcmd)

    # read in the almt
    header_pattern = re.compile(">\s*(\S+?)\s")
    aligned_seq = {}
    for line in ret.split('\n'):
        if ('error' in line):
            print " >>>>> ", line
            print fastafile1+" "+fastafile2
            print
            exit (1)
        if '>' in line:
            match   = re.search(header_pattern, line)
            if not match:
                print " >>>>> ", line
                print fastafile1+" "+fastafile2
                print
                exit (1)   
            species = match.groups()[0]
            aligned_seq[species]  = ""
        else:
            aligned_seq[species] += line

    commands.getoutput("rm "+fastafile1+" "+fastafile2)
       
    # find the positions of the exons in the alignment
    exon_positions = {}
    for species, seq in aligned_seq.iteritems():
        # move B to beginning of each exon sequence
        seq = moveB(seq)
        #beginning and end of each exon in the alignment
        exon_positions[species] = find_exon_positions(seq)

    # fill in the actual map values
    maps = maps_evaluate (relevant_human_exons, relevant_ortho_exons, aligned_seq, exon_positions)

 
    return maps
    
#########################################
def store (cursor, map):
    fixed_fields  = {}
    update_fields = {}
    fixed_fields ['exon_id']              = map.exon_id_2
    fixed_fields ['cognate_exon_id']      = map.exon_id_2
    fixed_fields ['cognate_genome_db_id'] = species2genome_db(map.species_2)
    update_fields['cigar_line']           = map.cigar_line
    update_fields['similarity']           = map.similarity
    update_fields['source']               = 'ensembl'
    #####
    store_or_update (cursor, 'exon_map', fixed_fields, update_fields)

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()
    acg    = AlignmentCommandGenerator()

    [all_species, ensembl_db_name] = get_species (cursor)

    species='homo_sapiens'
    switch_to_db (cursor,  ensembl_db_name[species])
    gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1)


    for gene_id in [310359]:

        print gene_id, gene2stable(cursor, gene_id), get_description (cursor, gene_id)
        
        # get _all_ exons
        switch_to_db (cursor, ensembl_db_name[species])
        human_exons = gene2exon_list(cursor, gene_id)
        if (not human_exons):
            print 'no exons for ', gene_id
            sys.exit(1)

        print "==============================="
        print "    one2one orthologues        "
        switch_to_db (cursor, ensembl_db_name[species])
        known_orthologues = get_orthos (cursor, gene_id)
        for [ortho_gene_id, ortho_species] in known_orthologues:
            print " ortho gene id ", ortho_gene_id, ortho_species
            ortho_exons = gene2exon_list(cursor, ortho_gene_id, db_name=ensembl_db_name[ortho_species] )
            if not ortho_exons:
                print 'no exons for ', species, ortho_gene_id
            maps = make_maps (cursor, ensembl_db_name, acg, ortho_species, human_exons, ortho_exons)            
            # store the maps into the database
            # store (cursor, maps, alignment_id)
        print "==============================="
            

        print "==============================="
        print "    unresolved orthologues     "
        switch_to_db (cursor, ensembl_db_name[species])
        unresolved_orthologues = get_unresolved_orthos (cursor, gene_id)
        for [ortho_gene_id, ortho_species] in unresolved_orthologues:
            print " ortho gene id ", ortho_gene_id, ortho_species
            ortho_exons = gene2exon_list(cursor, ortho_gene_id, db_name=ensembl_db_name[ortho_species] )
            if not ortho_exons:
                print 'no exons for ', species, ortho_gene_id
            maps = make_maps (cursor, ensembl_db_name, acg, ortho_species, human_exons, ortho_exons)
            # store the maps into the database
            # store (cursor, maps, alignment_id)
        print "==============================="
 

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
