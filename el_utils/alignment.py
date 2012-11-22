
from Bio.SubsMat import MatrixInfo as matlist

#########################################
def exon_aware_smith_waterman (seq1, seq2):


    matrix   = matlist.blosum62
    alphabet = map(chr, range(65, 91))
    digits   = map(chr, range(48, 58))

    for character in alphabet:

        if (character == 'B'):
            matrix[(character, 'B')] =  30
            matrix[('B', character)] =  30
            matrix[(character, 'Z')] = -30
            matrix[('Z', character)] = -30

        elif (character == 'Z'):
            matrix[(character, 'Z')] =  30
            matrix[('Z', character)] =  30
            matrix[(character, 'B')] = -30
            matrix[('B', character)] = -30
        else:
            matrix[(character, 'Z')] =  0
            matrix[('Z', character)] =  0
            matrix[(character, 'B')] =  0
            matrix[('B', character)] =  0

    for digit1 in digits:
        
        for digit2 in digits:
            matrix[(digit1, digit2)] = 10

        for character in alphabet:
            matrix[(digit1, character)] = 0
            matrix[(character, digit1)] = 0

        matrix[(digit1, 'Z')] =  -30
        matrix[('Z', digit1)] =  -30
        matrix[(digit1, 'B')] =  -30
        matrix[('B', digit1)] =  -30

        
    for any1 in alphabet+digits:
        for any2 in alphabet+digits:
            if (not matrix.has_key( (any1, any2))):
                if (not matrix.has_key( (any2, any1))):
                    matrix[(any1, any2)] = 0
                else:
                    matrix[(any1, any2)] = matrix[(any2, any1)]

    ######################################################3
    custom_gap_penalty_x = False
    custom_gap_penalty_y = False
    gap_opening   =  -5.0
    gap_extension =  -3.0
    endgap        = 0.0
    use_endgap    = False

    far_away = -1

    max_i = len(seq1)
    max_j = len(seq2)

    # allocation, initialization
    F         = [[0.0]*(max_j+1) for x in xrange(max_i+1)]
    direction = [[''] *(max_j+1) for x in xrange(max_i+1)]
    map_i2j   = [far_away]*max_i
    map_j2i   = [far_away]*max_j


    F_max   = far_away
    F_max_i = 0.0
    F_max_j = 0.0
  
    i_between_exons = True
    j_between_exons = True
    #
    for i in range(max_i+1):

        if i > 0:
            if seq1[i-1] == 'B':
                i_between_exons = False
            elif seq1[i-1] == 'Z':
                i_between_exons = True
       
        for  j in range(max_j+1):
         
            if j > 0:
                if seq2[j-1] == 'B':
                    j_between_exons = False
                elif seq2[j-1] == 'Z':
                    j_between_exons = True
               
	    if ( not i and not j ):
		F[0][0] = 0
		direction[i][j] = 'd'
		continue
	    
	    
	    if ( i and j ):
		if ( direction[i-1][j] == 'i' ) :
		    #  gap extension  */
		    if  j_between_exons:
			penalty =  0
                    else :
			if (use_endgap and j==max_j):
                            penalty = endgap 
                        else:
                            penalty = gap_extension
		    
                else :
		    #  gap opening  */
		    if  j_between_exons:
			penalty =  0
                    else:
			if (use_endgap and j==max_j):  
                            penalty = endgap 
                        else: 
                            penalty = gap_opening
		    
		
                i_sim =  F[i-1][j] + penalty
		
		if ( direction[i][j-1] == 'j' ) :
		    if i_between_exons:
			penalty = 0
                    else :
                        if (use_endgap and i==max_i):
                            penalty = endgap 
                        else:
                            penalty = gap_extension
		    
                else:
		    if i_between_exons:
			penalty = 0
                    else :
                        if (use_endgap and i==max_i):
                            penalty = endgap
                        else:
                            penalty = gap_opening
		    
		
                j_sim = F[i][j-1] + penalty
       	
		
		diag_sim =  F[i-1][j-1] + matrix[(seq1[i-1], seq2[j-1])]
		

		max_sim = diag_sim
		direction[i][j] = 'd'
		if ( i_sim > max_sim ):
		    max_sim = i_sim
		    direction[i][j] = 'i'
		
		if ( j_sim > max_sim ) :
		    max_sim = j_sim
		    direction[i][j] = 'j'
		
		
            elif j:
		
		if ( custom_gap_penalty_y and custom_gap_penalty_y[j-1] < 0 ) :
		    penalty = custom_gap_penalty_y[j-1]
                else:
		    if ( use_endgap) :
			penalty = endgap
                    else:
			if ( direction[i][j-1] =='j' ) :
			    penalty = gap_extension
                        else :
			    penalty = gap_opening
			
		    
		
		j_sim = F[i][j-1] + penalty
		max_sim = j_sim
		direction[i][j] = 'j'

		
            elif i:
		
		if ( custom_gap_penalty_x and custom_gap_penalty_x[i-1] < 0 ) :
		    penalty = custom_gap_penalty_x[i-1]
                else :
		    if ( use_endgap) :
			penalty = endgap
                    else :
			if ( direction[i-1][j] == 'i' ) :
			    penalty =  gap_extension
                        else :
			    penalty =  gap_opening
			
		    
		
		i_sim = F[i-1][j] + penalty
		max_sim = i_sim
		direction[i][j] = 'i'
		

	     

	    if (max_sim < 0.0 ): max_sim = 0.0
	    
	    F[i][j] = max_sim
	    if ( F_max < max_sim ) :
		# TODO: tie break here */
		F_max = max_sim
		F_max_i = i
		F_max_j = j
		
	 
    # retrace
    for i in range (F_max_i, max_i+1):
        map_i2j[i-1] = far_away
    for j in range (F_max_j, max_j+1):
        map_j2i[j-1] = far_away

    i = F_max_i
    j = F_max_j
    aln_score = F[i][j] 


    while ( i>0 or  j >0 ):

	if ( i<0 or j<0 ):
	    print "Retracing error"
            exit(1)
	
        if direction[i][j] == 'd':
	    map_i2j [i-1] = j-1
	    map_j2i [j-1] = i-1
	    i-= 1
	    j-= 1 
	    
	elif direction[i][j] == 'i':
	    map_i2j [i-1] = far_away
	    i-= 1 
	   
	elif direction[i][j] == 'j':
	    map_j2i [j-1] = far_away
	    j-= 1 
	   
        else: 
            print "Retracing error"
            exit(1)	
	
    
    print aln_score

    i = 0
    j = 0
    done = False
    aligned_seq_1 = ""
    aligned_seq_2 = ""
    while not done:
        if (map_i2j[i] == j):
            aligned_seq_1 += seq1[i]
            aligned_seq_2 += seq2[j]
            i += 1
            j += 1
        elif (map_i2j[i] < 0 ):
            aligned_seq_1 += seq1[i]
            aligned_seq_2 += '-'
            i += 1
        elif (map_j2i[j] < 0 ):
            aligned_seq_1 += '-'
            aligned_seq_2 += seq2[j]
            j += 1
           
        done = ( j>=max_j and i >= max_i)



    print ">seq1"
    print aligned_seq_1
    print ">seq2"
    print aligned_seq_2

    exit(1)

 
    return [map_i2j, map_j2i]


   

#########################################
#########################################
#########################################
#########################################
#########################################
def smith_waterman (similarity):

    
    custom_gap_penalty_x = False
    custom_gap_penalty_y = False
    gap_opening   = 0.0
    gap_extension = 0.0
    endgap        = 0.0
    use_endgap    = False

    far_away = -1

    max_i = len(similarity)
    max_j = len(similarity[0])

    # allocation, initialization
    F         = [[0.0]*(max_j+1) for x in xrange(max_i+1)]
    direction = [[''] *(max_j+1) for x in xrange(max_i+1)]
    map_i2j   = [far_away]*max_i
    map_j2i   = [far_away]*max_j


    F_max   = far_away
    F_max_i = 0.0
    F_max_j = 0.0
  
    i_between_exons = True
    j_between_exons = True
    for i in range(max_i+1):
        for  j in range(max_j+1):


	    if ( not i and not j ):
		F[0][0] = 0
		direction[i][j] = 'd'
		continue
	    
	    
	    if ( i and j ):
		if ( direction[i-1][j] == 'i' ) :
		    #  gap extension  */
		    if ( custom_gap_penalty_x and custom_gap_penalty_x[i-1] < 0 ) :
			penalty = custom_gap_penalty_x[i-1]
                    else :
			if (use_endgap and j==max_j):
                            penalty = endgap 
                        else:
                            penalty = gap_extension
		    
                else :
		    #  gap opening  */
		    if ( custom_gap_penalty_x and custom_gap_penalty_x[i-1] < 0 ) :
			penalty = custom_gap_penalty_x[i-1]
                    else :
			if (use_endgap and j==max_j):  
                            penalty = endgap 
                        else: 
                            penalty = gap_opening
		    
		
                i_sim =  F[i-1][j] + penalty
		
		if ( direction[i][j-1] == 'j' ) :
 		    if ( custom_gap_penalty_y and custom_gap_penalty_y[j-1] < 0 ) :
			penalty = custom_gap_penalty_y[j-1]
                    else :
                        if (use_endgap and i==max_i):
                            penalty = endgap 
                        else:
                            penalty = gap_extension
		    
                else:
		    if ( custom_gap_penalty_y and custom_gap_penalty_y[j-1] < 0 ) :
			penalty = custom_gap_penalty_y[j-1]
                    else :
                        if (use_endgap and i==max_i):
                            penalty = endgap
                        else:
                            penalty = gap_opening
		    
		
                j_sim = F[i][j-1] + penalty
       	
		
		diag_sim =  F[i-1][j-1] + similarity [i-1][j-1] 
		

		max_sim = diag_sim
		direction[i][j] = 'd'
		if ( i_sim > max_sim ):
		    max_sim = i_sim
		    direction[i][j] = 'i'
		
		if ( j_sim > max_sim ) :
		    max_sim = j_sim
		    direction[i][j] = 'j'
		
		
            elif j:
		
		if ( custom_gap_penalty_y and custom_gap_penalty_y[j-1] < 0 ) :
		    penalty = custom_gap_penalty_y[j-1]
                else:
		    if ( use_endgap) :
			penalty = endgap
                    else:
			if ( direction[i][j-1] =='j' ) :
			    penalty = gap_extension
                        else :
			    penalty = gap_opening
			
		    
		
		j_sim = F[i][j-1] + penalty
		max_sim = j_sim
		direction[i][j] = 'j'

		
            elif i:
		
		if ( custom_gap_penalty_x and custom_gap_penalty_x[i-1] < 0 ) :
		    penalty = custom_gap_penalty_x[i-1]
                else :
		    if ( use_endgap) :
			penalty = endgap
                    else :
			if ( direction[i-1][j] == 'i' ) :
			    penalty =  gap_extension
                        else :
			    penalty =  gap_opening
			
		    
		
		i_sim = F[i-1][j] + penalty
		max_sim = i_sim
		direction[i][j] = 'i'
		

	     

	    if (max_sim < 0.0 ): max_sim = 0.0
	    
	    F[i][j] = max_sim
	    if ( F_max < max_sim ) :
		# TODO: tie break here */
		F_max = max_sim
		F_max_i = i
		F_max_j = j
		
	 
    # retrace
    for i in range (F_max_i, max_i+1):
        map_i2j[i-1] = far_away
    for j in range (F_max_j, max_j+1):
        map_j2i[j-1] = far_away

    i = F_max_i
    j = F_max_j
    aln_score = F[i][j] 


    while ( i>0 or  j >0 ):

	if ( i<0 or j<0 ):
	    print "Retracing error"
            exit(1)
	
        if direction[i][j] == 'd':
	    map_i2j [i-1] = j-1
	    map_j2i [j-1] = i-1
	    i-= 1
	    j-= 1 
	    
	elif direction[i][j] == 'i':
	    map_i2j [i-1] = far_away
	    i-= 1 
	   
	elif direction[i][j] == 'j':
	    map_j2i [j-1] = far_away
	    j-= 1 
	   
        else: 
            print "Retracing error"
            exit(1)	
	
    

 
    return [map_i2j, map_j2i]
