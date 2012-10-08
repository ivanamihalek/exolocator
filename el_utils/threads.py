
import threading

def parallelize (no_threads, embarassingly_pllbl_fn, list, other_args):


    if (no_threads < 1):
        print "number of threads is expected to be >= 1"
        return False

    if (no_threads = 1):
        ret = embarassingly_pllbl_fn (list, other_args)
        return ret


    #########################################
    # nontrivial
    for chunk in range (1,len(list)+1):
        if ((no_threads+1)*chunk >= len(list)):
            break
    
    # check my math
    tot = 0
    for thr in range (no_threads):
        thr_from = thr*chunk
        thr_to   = thr_from+chunk
           
        if (thr_from >= len(list)):
            break
        if (thr == no_threads-1):
            thr_to = len(list)
            tot  += len(species_list[thr_from:thr_to])

        if ( len(species_list) != tot ):
            print " something's off with the balancing ... "
            return False

    # run
    for thr in range (no_threads):
        thr_from = thr*chunk
        thr_to   = thr_from+chunk
        if (thr_from >= len(list)):
            break
        if (thr == no_threads-1):
            thr_to = len(list)
                
        thread = threading.Thread (target=embarassingly_pllbl_fn, 
                                   args=(list[thr_from:thr_to], other_args))
        try:
            thread.start()
        except:
            print "Error: unable to start thread"
            return False
    
        
