import sys, os

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
