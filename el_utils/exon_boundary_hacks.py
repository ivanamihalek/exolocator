#########################################
def move_Z (peptide_alnmt, global_bdry_position, template_name):

    # see which boundaries can be re-aligned with human
    for name, seq in peptide_alnmt.iteritems():
        if name==template_name: continue
        pep_exons = peptide_alnmt[name].split ('-Z-')
        
        bdry_position  = 0
        for exon_ct in range(1,len(pep_exons)):
            bdry_position += len(pep_exons[exon_ct-1])
            if not bdry_position in global_bdry_position:
                for pos in global_bdry_position:
                    if bdry_position < pos:
                        gaps ='-'*(pos-bdry_position)
                        if peptide_alnmt[name][bdry_position+3:pos+3] ==gaps:
                            pepseq  = peptide_alnmt[name][:bdry_position]
                            pepseq += gaps+'-Z-'
                            pepseq += peptide_alnmt[name][pos+3:]
                            peptide_alnmt[name] = pepseq
                            bdry_new = False
                            break
                    else:
                       gaps ='-'*(bdry_position - pos) 
                       if peptide_alnmt[name][pos:bdry_position] ==gaps:
                            pepseq  = peptide_alnmt[name][:pos]
                            pepseq += '-Z-'+gaps
                            pepseq += peptide_alnmt[name][bdry_position+3:]
                            peptide_alnmt[name] = pepseq
                            bdry_new = False
                            break
            bdry_position += 3
    

#########################################
def correct_overlap_Z(peptide_alnmt, name, pos, global_bdry_position):
    #pdb.set_trace()
    #print name
    
    length       = len(peptide_alnmt[name])
    if pos+3 >= length: return
    triple       = not '-' in peptide_alnmt[name][pos:pos+3]
    insert_left  = peptide_alnmt[name][pos] == '-' 
    insert_left  = insert_left or  triple and (not pos or  peptide_alnmt[name][pos-1]=='-')
    insert_right = peptide_alnmt[name][pos+2] == '-'
    insert_right = insert_right or  triple and (not pos>=length-3 or peptide_alnmt[name][pos+3]=='-')
    if insert_left:
        # insert gaps on the left
        if triple:
            insert = '---'
            insert_start = pos

        elif peptide_alnmt[name][pos:pos+2] == '--':
            insert = '-'
            insert_start = pos+2
        else:
            insert = '--'
            insert_start = pos+1

        for name2, seq2 in peptide_alnmt.iteritems():
            if peptide_alnmt[name2][pos:pos+3] == '-Z-':
                pepseq  = seq2[:pos+3]
                pepseq += insert
                pepseq += seq2[pos+3:]
                peptide_alnmt[name2] = pepseq
            else:
                pepseq  = seq2[:insert_start]
                pepseq += insert
                pepseq += seq2[insert_start:]
                peptide_alnmt[name2] = pepseq
        # also: shift the boundary positions
        for i in range(len(global_bdry_position) ):
            if global_bdry_position[i] <= pos: continue
            global_bdry_position[i] += len(insert)

    elif insert_right:

        # insert gaps on the right
        if triple:
            insert = '---'
            insert_start = pos+3
        elif peptide_alnmt[name][pos+1:pos+3] == '--':
            insert = '-'
            insert_start = pos+1
        else:
            insert = '--'
            insert_start = pos+2


        for name2, seq2 in peptide_alnmt.iteritems():
            if peptide_alnmt[name2][pos:pos+3] == '-Z-':
                pepseq  = seq2[:pos]
                pepseq += insert
                pepseq += seq2[pos:]
                peptide_alnmt[name2] = pepseq
            else:
                pepseq  = seq2[:insert_start]
                pepseq += insert
                pepseq += seq2[insert_start:]
                peptide_alnmt[name2] = pepseq

        # also: shift the boundary positions
        for i in range(len(global_bdry_position) ):
            if global_bdry_position[i] < pos: continue
            global_bdry_position[i] += len(insert)


#########################################
def splice_gaps(alnmt, to_splice):

    number_of_previous_splices = 0
    for pos in to_splice:
        pos_corrected = pos + number_of_previous_splices

        number_of_previous_splices += 1
        for name in alnmt.keys():
            temp = alnmt[name]
            alnmt[name] = temp[:pos_corrected]+'-'+temp[pos_corrected:]


#########################################
# what I really need is the alignment program that enforces boundaries ...
def boundary_cleanup (peptide_alnmt, sorted_seq_names): 

    global_bdry_position = []

    # check that no boundary marker is missing one '-'
    delimiter = re.compile("Z")
    
    for name, seq in peptide_alnmt.iteritems():
        to_splice = []
        for match in delimiter.finditer(seq):
            if not match: continue
            if match.start() and not seq[match.start()-1] == '-':
                to_splice.append(match.start())
            if  match.start()+1 < len(seq) and not seq[match.start()+1] == '-':
                to_splice.append(match.start()+1)
        splice_gaps(peptide_alnmt, to_splice)

    # use sorted species as anchors
    for sorted_seq_name in sorted_seq_names:
        pep_exons = peptide_alnmt[sorted_seq_name].split ('-Z-')
        bdry_position  = 0
        for exon_ct in range(1,len(pep_exons)):
            bdry_position += len(pep_exons[exon_ct-1])
            global_bdry_position.append(bdry_position)
            bdry_position += 3
 
        global_bdry_position = sorted(global_bdry_position)

        # see which boundaries can be re-aligned with anchor's Z
        move_Z (peptide_alnmt, global_bdry_position, sorted_seq_name)

        # find sequences that have amino acid overlapping with the 
        # exon boundary marker for our "anchor" sequence
        for name, seq in peptide_alnmt.iteritems():
            if name==sorted_seq_name: continue  # this is the anchor

            for pos in global_bdry_position:
                # this one has the bounday marker itself
                if 'Z' in peptide_alnmt[name][pos:pos+3]:
                    continue
                # this one is not a problem
                if peptide_alnmt[name][pos:pos+3] == '---':
                    continue
                # this needs correction
                correct_overlap_Z(peptide_alnmt, name, pos, global_bdry_position)

        # see which boundaries can be re-aligned with anchor's Z
        move_Z (peptide_alnmt, global_bdry_position, sorted_seq_name)
            
#########################################
def find_exon_boundaries (peptide_alignment):

    global_bdry_position = []
    local_bdry_position  = {}

    delimiter = re.compile("-Z-")
    for name, seq in peptide_alignment.iteritems():
        local_bdry_position[name] = []
        for match in delimiter.finditer(seq):
            if not match: continue
            pos = match.start()
            if not pos in global_bdry_position:
                global_bdry_position.append(pos)
            if not pos in local_bdry_position[name]:
                local_bdry_position[name].append(pos)

    global_bdry_position = sorted(global_bdry_position)

    return global_bdry_position, local_bdry_position

