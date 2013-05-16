#!/usr/bin/python

import MySQLdb, subprocess, re, commands


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
                [template_start, template_end]  =  map (int,re.split('\D+', fields[4])[:2])

                table_entry[target] = [search_start, search_end, template_start, template_end]

            continue

        elif ' Query ' in lines[lineno]: 
            read_start = lineno
            match =  re.search('\d+', lines[lineno])
            seqlen  = int(match.group())

        elif 'Evalue' in lines[lineno] and  read_start >0: # the result is in the previous couple of lines
            
            [matchlen, number_matching, identity] = map(int, re.split('\D+', lines[lineno])[:3])
            if  matchlen < 0.4*seqlen or identity < 10: continue
            print matchlen, number_matching, identity
                                                        
            # FOUND AN EXON!
            # ... but lets keep what might be the best match
            if matchlen > longest:
                #[search_start, search_end, template_start, template_end] = table_entry[target]

                aligned_qry_seq = ""
                for row in lines[read_start+3:lineno-3:4]: # every fourth row
                    seq = re.split('\s+',row)[3]
                    print seq
                    aligned_qry_seq += seq
                aligned_target_seq = ""
                for row in  lines[read_start+5:lineno-1:4]: # every fourth row
                    seq = re.split('\s+',row)[3]
                    aligned_target_seq += seq

                # just one more piece of hacking:
                # we had to chopup the query into smaller pieces to make usearch happy
                # where is the match, counted from the beginning of the whole thing?
                best_match = [search_start-1, search_end-1, template_start-1, 
                              template_end-1, aligned_target_seq, aligned_qry_seq]
                

    print best_match

    return best_match




#############################################################

resultstr = commands.getoutput ("cat in.test")
parse_usearch_output(resultstr)
