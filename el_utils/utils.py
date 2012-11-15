import sys

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
def output_fasta (filename, headers, sequences):
    outf = erropen (filename, "w")
    for i in range (len(headers)):
        print >> outf, ">"+headers[i]
        print >> outf, sequences[i]
    outf.close()

    return
