import sys, os
import string

#########################################
def erropen (file,mode):
    of = None
    try:
        of = open (file,mode)
    except:
        print "error opening ", file
        sys.exit(1)

    return of

#########################################
def mkdir_p (path):
    try:
        os.makedirs(path)
    except: 
        sys.exit(1) 

#########################################
def output_fasta (filename, headers, sequence):
    outf = erropen (filename, "w")
    for header  in  headers:
        if not sequence.has_key(header): continue
        print >> outf, ">"+header
        print >> outf, sequence[header]
    outf.close()

    return
#########################################
def input_fasta (filename):
    sequence = {}
    header   = ""
    inf = erropen (filename, "r")
    for line  in  inf:
        if '>' in line:
            header   = line.rstrip().replace('>', "").replace(' ', "")
            sequence[header] = ""
        elif header: 
            sequence[header] += line.rstrip()
       
    inf.close()

    return sequence

#########################################
def parse_aln_name (name):
    fields     = name.split("_")
    exon_id    = int(fields[-2])
    exon_known = int(fields[-1])
    species    =  "_".join(fields[:-2])
    return [species, exon_id, exon_known]

#########################################
def  fract_identity (cigar_line):

    fraction = 0

    char_pattern = re.compile("\D")
    total_length = 0
    common       = 0
    prev_end     = 0
    lengthA      = 0
    lengthB      = 0
    for match in char_pattern.finditer(cigar_line):
        this_start       = match.start()
        no_repeats = int(cigar_line[prev_end:this_start])
        alignment_instruction = cigar_line[this_start]
        prev_end = match.end()

        total_length += no_repeats
        if alignment_instruction == 'M':
            common  += no_repeats
            lengthA += no_repeats
            lengthB += no_repeats
        elif alignment_instruction == 'A':
            lengthB += no_repeats
        elif alignment_instruction == 'B':
            lengthA += no_repeats
            
    shorter = lengthA if lengthA<=lengthB else lengthB

    if shorter == 0: return fraction # fraction is still set to 0

    if total_length:
        fraction = common/float(shorter)
        
    return  fraction

#########################################
def  pairwise_fract_identity (seq1, seq2):
    
    fract_identity = 0.0
    if ( not len(seq1)):
        return fract_identity

    for i in range(len(seq1)):
        if (seq1[i] == '-'): continue
        if seq1[i] == seq2[i]: fract_identity += 1.0
    
    fract_identity /= float(len(seq1))
    return fract_identity

#########################################
def  pairwise_fract_similarity (seq1, seq2):
    
    is_similar_to = {}

    # this is rather crude ...
    for  char in string.printable:
        is_similar_to[char] = char

    is_similar_to['I'] = 'V';
    is_similar_to['L'] = 'V';
    is_similar_to['S'] = 'T';
    is_similar_to['D'] = 'E';
    is_similar_to['K'] = 'R';
    is_similar_to['Q'] = 'N';
    is_similar_to['.'] = '.';
    is_similar_to['-'] = '.';

    is_similar_to['A'] = 'V';
    is_similar_to['M'] = 'V';
    is_similar_to['G'] = 'V';
    is_similar_to['F'] = 'Y';
    is_similar_to['H'] = 'R';

    fract_similarity = 0.0
    if ( not len(seq1)):
        return fract_similarity

    for i in range(len(seq1)):
        if (seq1[i] == '-'): continue
        if is_similar_to[seq1[i]] == is_similar_to[seq2[i]]: fract_similarity += 1.0
    
    fract_similarity /= float(len(seq1))
    return fract_similarity

