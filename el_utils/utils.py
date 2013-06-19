import sys, os,  re, commands
import string, random
from subprocess import Popen, PIPE, STDOUT
from tempfile   import NamedTemporaryFile
from math       import sqrt


#########################################
def unfold_cigar_line (seq_A, seq_B, cigar_line):

    seq_A_aligned = ""
    seq_B_aligned = ""

    char_pattern = re.compile("\D")
    a_ct     = 0
    b_ct     = 0
    prev_end = 0

    for match in char_pattern.finditer(cigar_line):
        this_start       = match.start()
        no_repeats = int(cigar_line[prev_end:this_start])
        alignment_instruction = cigar_line[this_start]
        prev_end = match.end()

        if alignment_instruction == 'M':
            seq_A_aligned += seq_A[a_ct:a_ct+no_repeats]
            a_ct         += no_repeats
            seq_B_aligned += seq_B[b_ct:b_ct+no_repeats]
            b_ct  += no_repeats

        elif alignment_instruction == 'A':
            seq_A_aligned += '-'*no_repeats 
            seq_B_aligned += seq_B[b_ct:b_ct+no_repeats]
            b_ct  += no_repeats
            
        elif alignment_instruction == 'B':
            seq_A_aligned += seq_A[a_ct:a_ct+no_repeats]
            a_ct          += no_repeats
            seq_B_aligned += '-'*no_repeats 

    return [seq_A_aligned, seq_B_aligned]


#########################################
def sw_search (cfg, acg, query_seq, target_seq, delete= True):

    resultstr  = ""

    # save fasta (temporarily)
    random_str   = ''.join(random.choice(string.ascii_letters + string.digits) for x in range(6))

    qry_filename = "{0}/qry_{1}.fa".format  (cfg.dir_path['scratch'], random_str)
    qry_file = erropen( qry_filename, "w")
    qry_file.write(">query\n"+query_seq+"\n")
    qry_file.close()


    tgt_filename = "{0}/tgt_{1}.fa".format  (cfg.dir_path['scratch'], random_str)
    tgt_file = erropen( tgt_filename, "w")
    tgt_file.write(">tartget\n"+target_seq+"\n")
    tgt_file.close()


    # do  SW# search
    swsharpcmd = acg.generate_SW_nt (qry_filename, tgt_filename)
    print swsharpcmd
    resultstr  = commands.getoutput (swsharpcmd)


    if 'Segmentation' in  resultstr:
        print swsharpcmd
        print  " ** ", resultstr
        resultstr = ""

    if delete:
        cmd = "rm  {0}  {1} ".format (qry_filename, tgt_filename)
        stdout  = commands.getoutput (cmd)


    return resultstr


#########################################
def usearch (cfg, acg, query_seq, target_seq, delete= True):

    resultstr  = ""

    # NmedTMpFile was giving me some of the wirdest errors with usearch
    # something about amino acid files, but when I ran it from cmd line it was ok
    
    random_str   = ''.join(random.choice(string.ascii_letters + string.digits) for x in range(6))

    qry_filename = "{0}/qry_{1}.fa".format  (cfg.dir_path['scratch'], random_str)
    qry_file = erropen( qry_filename, "w")
    qry_file.write(">query\n"+query_seq+"\n")
    qry_file.close()


    total_length = 0
    outstr    = ""
    chunksize = 15000

    while total_length < len(target_seq):

        if (total_length+chunksize >= len(target_seq)):
            upper_bound   = len(target_seq)
        else:
            upper_bound   = total_length+chunksize

        outstr += ">piece_%d_%d"% (total_length, upper_bound)
        outstr += "\n"
        outstr += target_seq[total_length:upper_bound]
        if not outstr[-1] == '\n': outstr += "\n"
   
        total_length += chunksize


    total_length = 10000
    while total_length < len(target_seq):

        if (total_length+chunksize >= len(target_seq)):
            upper_bound   = len(target_seq)
        else:
            upper_bound   = total_length+chunksize

        outstr += ">piece_%d_%d"% (total_length, upper_bound)
        outstr += "\n"
        outstr += target_seq[total_length:upper_bound]
        if not outstr[-1] == '\n': outstr += "\n"
   
        total_length += chunksize

    name = "tgt_blah"
    tgt_filename = "{0}/tgt_{1}.fa".format  (cfg.dir_path['scratch'], random_str)
    tgt_file = erropen( tgt_filename, "w")
    tgt_file.write(outstr)
    tgt_file.close()
 
    # 
    outname = "{0}/{1}.out".format  (cfg.dir_path['scratch'], random_str)
   

    # do  usearch
    cmd = acg.generate_usearch_nt (qry_filename, tgt_filename, outname)
    stdout  = commands.getoutput (cmd)
    if 'Segmentation' in  stdout:
        print  cmd
        print  " ** ", stdout
        resultstr = ""

    # what happens if the search fails?
    resultstr =  commands.getoutput ("cat "+ outname)
    
    if delete:
        cmd = "rm  {0}  {1}  {2}".format (qry_filename, tgt_filename, outname)
        stdout  = commands.getoutput (cmd)


    return resultstr


#########################################
def parse_sw_output (resultstr):
		    
    best_match = None
    longest    = -1

    for r in (f.splitlines() for f in resultstr.split("#"*80+"\n")):
        
        if len(r) < 14: continue # Skip blank or malformed results

        # Parse result
        seqlen = min(int(re.split('\D+',r[1])[1]),int(re.split('\D+',r[3])[1]))
        identity, matchlen = map(int, re.split('\D+', r[7])[1:3])
        #similarity = int(re.split('\D+',r[8])[1])
        #gaps       = int(re.split('\D+',r[9])[1])
        #score      = float(r[10].split()[1])

        
        # Reject if identity too low or too short -- seqlen is the length of the query
        if identity < 0.4*seqlen or identity < 10: continue

        # FOUND AN EXON!
        # ... but lets keep what might be the best match
        if matchlen > longest:
            longest = matchlen
            [search_start, search_end, template_start, template_end] = map(int,re.split('\D+',r[6])[1:5])

            aligned_qry_seq = ""
            for row in r[13::3]: # every third row
                seq = re.split('\s+',row)[2]
                aligned_qry_seq += seq
            aligned_target_seq = ""
            for row in r[12::3]: # every third row
                seq = re.split('\s+',row)[2]
                aligned_target_seq += seq
            best_match = [search_start-1, search_end-1, template_start-1, template_end-1, 
                          aligned_target_seq, aligned_qry_seq]

    return best_match 


#########################################
def parse_usearch_output (resultstr):
	

    best_match = None
    longest    = -1

    read_start = -1
    lines  =  resultstr.split("\n")
    table_entry = {}

    for lineno in range(len(lines)):
        if 'QueryLo-Hi' in lines[lineno]: 
            for ct in range  (lineno+1,len(lines)):
            
                fields = re.split('\s+', lines[ct])
                if len(fields) < 7: break
                target =  fields[-1]

                [search_start, search_end]      =  map (int, re.split('\D+', fields[5])[:2])
                [template_start, template_end]  =  map (int, re.split('\D+', fields[4])[:2])

                table_entry[target] = [search_start, search_end, template_start, template_end]

        elif ' Query ' in lines[lineno]: 
            read_start = lineno
            match =  re.search('\d+', lines[lineno])
            seqlen  = int(match.group())

        elif 'Evalue' in lines[lineno] and  read_start >0: # the result is in the previous couple of lines
            
            [matchlen, number_matching, identity] = map(int, re.split('\D+', lines[lineno])[:3])
            if  matchlen < 0.4*seqlen or identity < 10: continue
                                                        
            # FOUND AN EXON!
            # ... but lets keep what might be the best match
            if matchlen > longest:

                longest = matchlen

                target_line_fields = lines[read_start+1].split (" ")
                target = target_line_fields[-1][1:] # get rid of ">"
                
                [search_start, search_end, template_start, template_end] = table_entry[target]

                aligned_qry_seq = ""
                for row in lines[read_start+3:lineno-3:4]: # every fourth row
                    seq = re.split('\s+',row)[3]
                    aligned_qry_seq += seq
                aligned_target_seq = ""
                for row in lines[read_start+5:lineno-1:4]: # every fourth row
                    seq = re.split('\s+',row)[3]
                    aligned_target_seq += seq

                # just one more piece of hacking:
                # we had to chopup the query into smaller pieces to make usearch happy
                # where is the match, counted from the beginning of the whole thing?
                name_fields = target.split ("_")
                piece_from  = int(name_fields[-2])
                #piece_to   = int(name_fields[-1])
                best_match  = [piece_from+search_start-1, piece_from+search_end-1, template_start-1, 
                               template_end-1, aligned_target_seq, aligned_qry_seq]
                

    return best_match 



#########################################
def get_best_filename(names_string):
    names_list = names_string.split()
    chromo_names = [n for n in names_list if 'chromosom' in n]
    if chromo_names: return chromo_names[0]
    else: return names_list[0]

#########################################
def get_fasta (acg, species, searchname, searchfile, searchstrand, searchstart, searchend):

    fasta    = None
    fastacmd = acg.generate_fastacmd_gene_command(species, searchname, searchfile, 
                                                  searchstrand, searchstart, searchend)
    p      = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
    fasta, errmsg = p.communicate()
    if errmsg and 'From location cannot be greater than' in errmsg:
        match = re.search('\d+', errmsg )
        max_searchend = match.group(0)
        if max_searchend < searchstart: return fasta

        fastacmd = acg.generate_fastacmd_gene_command(species, searchname, 
                                                      searchfile, searchstrand, searchstart, max_searchend)
        p      = Popen(fastacmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
        fasta, errmsg = p.communicate()


    return fasta


#########################################
def check_seq_length(sequence, msg):

    if not sequence.values():
        print "no sequences (?)"
        return False
    aln_length = len(sequence.values()[0])
    if not aln_length:
        print "aln length zero (?)"
        return False
    for name, seq in sequence.iteritems():
        if not len(seq) == aln_length:
            print msg, 
            print "seq length check failure  for",  name, " length: ", len(seq),  "aln_length", aln_length
            afa_fnm = msg+'.afa'
            print "writing the offending alnmt to ", afa_fnm
            output_fasta (afa_fnm, sequence.keys(), sequence)
            return False
    return True

#########################################
def strip_gaps (sequence):

    seq_stripped = {}

    all_gaps = {}  

    if not check_seq_length(sequence, 'in_strip_gaps'): 
        return ""
    
    aln_length = len(sequence.itervalues().next())

    if aln_length is None or aln_length==0:
        print "aln length zero (?)"
        return sequence
    
    for name, seq in sequence.iteritems():
        if not len(seq): 
            continue
        sequence[name] = seq.replace("-Z-", "BZB")

    for pos in range(aln_length):
        all_gaps[pos] = True
        for name, seq in sequence.iteritems():
            if not len(seq): 
                continue
            if (not seq[pos]=='-'):
                all_gaps[pos] = False
                break


    for name, seq in sequence.iteritems():
        if not len(seq): 
            continue
        seq_stripped[name] = ""
        for pos in range(aln_length):
            if all_gaps[pos]: continue
            seq_stripped[name] += seq[pos]


    for name, seq in seq_stripped.iteritems():
        if not len(seq): 
            continue
        seq_stripped[name] = seq_stripped[name].replace("BZB", "-Z-")

    return seq_stripped

#########################################
def erropen (file,mode):
    of = None
    try:
        of = open (file,mode)
    except:
        print "error opening ", file
        return None

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
    if not outf: return False

    for header  in  headers:
        if not sequence.has_key(header): continue
        print >> outf, ">"+header
        chunk_size   = 50
        chunk_number =  1
        while chunk_number*chunk_size <= len (sequence[header]):
            start = (chunk_number-1)*chunk_size
            end = start+chunk_size
            print >> outf, sequence[header][start:end]
            chunk_number += 1
        if chunk_number*chunk_size > len (sequence[header]):
            start = (chunk_number-1)*chunk_size
            print >> outf, sequence[header][start:]
    outf.close()

    return True

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
            
    if lengthA<=lengthB:
        shorter = lengthA
    else:
        shorter = lengthB
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

    common_length = 0.0
    for i in range(len(seq1)):
        if (seq1[i] == '-' or seq2[i] == '-'): continue
        if is_similar_to[seq1[i]] == is_similar_to[seq2[i]]: fract_similarity += 1.0
        common_length += 1.0

    if not common_length: return 0.0

    fract_similarity /= common_length
    return fract_similarity

#########################################
def  pairwise_tanimoto (seq1, seq2):
    
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

    tanimoto = 0.0
    if ( not len(seq1) or  not len(seq1)):
        return tanimoto

    similar_length = 0.0
    common_length = 0.0
    for i in range(len(seq1)):
        if (seq1[i] == '-' or seq2[i] == '-'): continue
        if is_similar_to[seq1[i]] == is_similar_to[seq2[i]]: similar_length += 1.0
        common_length += 1.0

    if not common_length: return 0.0

    tanimoto = sqrt(float(similar_length*similar_length)/(len(seq1)*len(seq2)))
    return tanimoto

