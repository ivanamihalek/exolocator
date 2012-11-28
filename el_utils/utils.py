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
        print >> outf, ">"+header
        print >> outf, sequence[header]
    outf.close()

    return
