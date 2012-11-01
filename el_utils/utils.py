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
