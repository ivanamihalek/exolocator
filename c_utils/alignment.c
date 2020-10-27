/* see https://docs.python.org/3.8/extending/extending.html*/
/* compilation: */
/* gcc -shared -I/usr/include/python3.8/ -lpython3.8 -o alignment.so alignment.c -fPIC -O */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define METH_VARARGS  0x0001




# define ASCII_SIZE  128


/***************************************************/
void * emalloc(int	size)
{
	void * ptr;
	if ((ptr = calloc(size, 1)) == NULL) {
        fprintf (stderr,  "emalloc: no memory for %u bytes", size);
        return NULL;
	}

	return ptr;
}
/***************************************************/
char **cmatrix(int rows, int columns){
	char **m;
	int i;
	/* allocate pointers to rows */
	m=(char **) malloc(rows*sizeof(char*));
	if (!m)  {
        fprintf (stderr,"row allocation failure  in chmatrix().\n");
        return NULL;
	}
	/* allocate rows and set pointers to them */
	m[0]=(char *) calloc( rows*columns, sizeof(char));
	if (!m[0]) {
        fprintf (stderr,"column allocation failure in chmatrix().\n");
        return NULL;
	}
	for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
	/* return pointer to array of pointers to rows */
	return m;
}
/***************************************************/
int **imatrix(int rows, int columns){
	int **m;
	int i;
	/* allocate pointers to rows */
	m=(int **) malloc(rows*sizeof(int*));
	if (!m)  {
        fprintf (stderr,"row allocation failure  in chmatrix().\n");
        return NULL;
	}
	/* allocate rows and set pointers to them */
	m[0]=(int *) calloc( rows*columns, sizeof(int));
	if (!m[0]) {
        fprintf (stderr,"column allocation failure in chmatrix().\n");
        return NULL;
	}
	for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
	/* return pointer to array of pointers to rows */
	return m;
}
/***************************************************/
void free_imatrix(int **m) {
	free(m[0]);
	free(m);
}
void free_cmatrix(char **m) {
	free(m[0]);
	free(m);
}
/***************************************************/
int load_sim_matrix (int ** similarity) {

	char *amino_acid_order = "ABCDEFGHIKLMNPQRSTVWXYZ";
	char *digits           = "0123456789";
	int aao_strlen, ctr, i, j;
	int char_i, char_j;
	int DELIMITER_SCORE = 5;
	int blosum62[]={
	 4,
	-2,  4,
	 0, -3,  9,
	-2,  4, -3,  6,
	-1,  1, -4,  2,  5,
	-2, -3, -2, -3, -3,  6,
	 0, -1, -3, -1, -2, -3,  6,
	-2,  0, -3, -1,  0, -1, -2,  8,
	-1, -3, -1, -3, -3,  0, -4, -3,  4,
	-1,  0, -3, -1,  1, -3, -2, -1, -3,  5,
	-1, -4, -1, -4, -3,  0, -4, -3,  2, -2,  4,
	-1, -3, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5,
	-2,  3, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6,
	-1, -2, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7,
	-1,  0, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,
	-1, -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5,
	 1,  0, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,
	 0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,
	 0, -3, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4,
	-3, -4, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,
	 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1, -2, -1,
	-2, -3, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2, -1,  7,
	-1,  1, -3,  1,  4, -3, -2,  0, -3,  1, -3, -1,  0, -1,  3,  0,  0, -1, -2, -3, -1, -2,  4};

    /* first, the plain old blosum values for standard amino acid types */
	aao_strlen = strlen(amino_acid_order);
	ctr = 0;
	for(i=0;i<aao_strlen;i++){
		char_i = (int) amino_acid_order [i];
		for (j=0;j<=i;j++){
			char_j = (int) amino_acid_order [j];
			/* similarity matrix is dimensioned to hold all ascii characters */
			similarity[char_i][char_j] = similarity[char_j][char_i] = blosum62[ctr];
			ctr++;
		}
	}

    /* I am repurposing X, B, and Z characters */
	/* X is neutral - this is unknown amino acid*/
	char_i = (int) 'X';
	for (j=0;j<=i;j++){
		char_j = (int) amino_acid_order [j];
		similarity[char_i][char_j] = similarity[char_j][char_i] = 0;
	}

	/* special treatment for B, Z */
	/* the sequence passed to thei function consists of concatenated exons of the format
        padded_count = "{0:03d}".format(count) (the count id presumably the order in which the exon appears)
        decorated_seq += 'B' + padded_count + pepseq + 'Z'
        for exaple, for pepseq = "JKHHG", count=5, the exon is decorated into 'B005JKHHGZ' */
	char_i = (int) 'B'; /* B is the left flank of an exon*/
	for (j=0;j<=i;j++){
        char_j = (int) amino_acid_order [j];
        if ( char_j == 'B' ) {
             /* B is very similar to B */
            similarity[char_i][char_j] = similarity[char_j][char_i] =  DELIMITER_SCORE;
        } else if ( char_j == 'Z' ) {
           /* Z definet;y does not want to be aligned with B' if  DELIMITER_SCORE,
                the penalty of -10 is the largest anywhere*/
             similarity[char_i][char_j] = similarity[char_j][char_i] = -2*DELIMITER_SCORE;
        } else {
            similarity[char_i][char_j] = similarity[char_j][char_i] =  0;
        }
	}
	char_i = (int) 'Z'; /* Z is the right flank of an exon*/
	for (j=0;j<=i;j++){
        char_j = (int) amino_acid_order [j];
        if ( char_j == 'Z') {
            /* Z is very similar to Z */
            similarity[char_i][char_j] = similarity[char_j][char_i] =  DELIMITER_SCORE;
        } else if ( char_j == 'B') {
            /* Z definet;y does not want to be aligned with B' if  DELIMITER_SCORE,
                the penalty of -10 is the largest anywhere - didn't we atke care of that above?*/
            similarity[char_i][char_j] = similarity[char_j][char_i] = -2*DELIMITER_SCORE;
        } else {
            similarity[char_i][char_j] = similarity[char_j][char_i] =  0;
        }
	}

	int digit_strlen = strlen(digits);

	for(i=0;i<digit_strlen;i++){
        char_i = digits[i];
        /* award to align digit to digit */
        for(j=0;j<digit_strlen;j++){
            char_j = digits[j];
            similarity[char_i][char_j] = similarity[char_j][char_i] = DELIMITER_SCORE;
        }
        /* digit to amino acid */
        /* I cannot penalize that too extremely, because one exon might be longer than the other */
        for(j=0;j<aao_strlen;j++){
            char_j = amino_acid_order[j];
            similarity[char_i][char_j] = similarity[char_j][char_i] = -DELIMITER_SCORE/2;
        }
        /* Z and B to digit*/
        char_j = 'Z';
        similarity[char_i][char_j] = similarity[char_j][char_i] = -3*DELIMITER_SCORE;
        char_j = 'B';
        similarity[char_i][char_j] = similarity[char_j][char_i] = -3*DELIMITER_SCORE;
    }

	return 0;
}

/***************************************************/
/*
 * Function to be called from Python
 */
static PyObject* smith_waterman_context(PyObject* self, PyObject* args)
{
	char *seq1 = NULL;
	char *seq2 = NULL;
	char retstr[100]   = {'\0'};
	int  len1, len2;
	int i, j;
	int gap_opening, gap_extension;
	static int ** similarity = NULL;


	PyArg_ParseTuple(args, "s#s#ii", &seq1, &len1, &seq2, &len2, &gap_opening, &gap_extension);

	if (!seq1 || !seq2) {
		sprintf (retstr, "no seq in py_smith_waterman_context");
		return Py_BuildValue("s", retstr);
	}

	/* passing a matrix this way is all to painful, so we'll elegantyly hardcode it: */
	if ( !similarity) {
		similarity = imatrix(ASCII_SIZE, ASCII_SIZE);
		if (!similarity) {
			sprintf (retstr, "error alloc matrix space");
			return Py_BuildValue("s", retstr);
		}
		load_sim_matrix (similarity);
   }


	/**********************************************************************************/
	//int gap_opening   =  -5; // used in 15_make_maps
	//int gap_extension =  -3;
	//char gap_character = '-'
	//int gap_opening    =  -3;  // used in 25_db_migration/06_make_alignments
	//int gap_extension  =   0;
	char gap_character = '#'; // why am I using this as a gap character?
	int endgap         =   0;
	int use_endgap     =   0;

	int far_away = -1;

	int max_i    = len1;
	int max_j    = len2;

	// allocation, initialization
	int  **F         = NULL;
	char **direction = NULL;
	int *map_i2j     = NULL;
	int *map_j2i     = NULL;

	if ( ! (F= imatrix (max_i+1, max_j+1)) ) {
		sprintf (retstr, "error alloc matrix space");
		return Py_BuildValue("s", retstr);
	}
	if ( ! (direction = cmatrix (max_i+1, max_j+1)) ) {
		sprintf (retstr, "error alloc matrix space");
		return Py_BuildValue("s", retstr);
	}
	if (! (map_i2j = emalloc( (max_i+1)*sizeof(int))) ) {
		sprintf (retstr, "error alloc matrix space");
		return Py_BuildValue("s", retstr);
	}
	if (! (map_j2i = emalloc( (max_j+1)*sizeof(int))) ) {
		sprintf (retstr, "error alloc matrix space");
		return Py_BuildValue("s", retstr);
	}
	for (i=0; i<=max_i; i++) map_i2j[i]=far_away;
	for (j=0; j<=max_j; j++) map_j2i[j]=far_away;


	int F_max   = far_away;
	int F_max_i = 0;
	int F_max_j = 0;
	int penalty = 0;
	int i_sim, j_sim, diag_sim, max_sim;

	int i_between_exons = 1;
	int j_between_exons = 1;
	//
	for (i=0; i<=max_i; i++) {

		if (i > 0) {
			if (seq1[i-1] == 'B') {
				i_between_exons = 0;
			} else if ( seq1[i-1] == 'Z'){
				i_between_exons = 1;
			}
		}
		for (j=0; j<=max_j; j++) {

			if (j > 0) {
				if (seq2[j-1] == 'B') {
					j_between_exons = 0;
				} else if (seq2[j-1] == 'Z') {
					j_between_exons = 1;
				}
			}

			if ( !i && !j ){
				F[0][0] = 0;
				direction[i][j] = 'd';
				continue;
			}

			/* neither i nor j are 0 */
			if ( i && j ){

				/**********************************/
				penalty =  0;
				if ( direction[i-1][j] == 'i' ) {
					//  gap extension
					if  (j_between_exons) {
						penalty =  0;
					} else {
						if (use_endgap && j==max_j){
							penalty = endgap;
						} else {
							penalty = gap_extension;
						}
					}
				} else {
					//  gap opening  */
					if  (j_between_exons) {
						penalty =  0;
					} else {
						if (use_endgap && j==max_j){
							penalty = endgap;
						} else{
							penalty = gap_opening;
						}
					}
				}
				i_sim =  F[i-1][j] + penalty;

				/**********************************/
				penalty =  0;
				if ( direction[i][j-1] == 'j' ) {
					//  gap extension
					if (i_between_exons) {
						penalty = 0;
					} else {
						if (use_endgap && i==max_i){
							penalty = endgap;
						} else{
							penalty = gap_extension;
						}
					}
				} else {
					//  gap opening  */
					if  (i_between_exons) {
					   penalty =  0;
					} else {
						if (use_endgap && i==max_i){
							penalty = endgap;
						} else {
							penalty = gap_opening;
						}
					}

				}
				j_sim = F[i][j-1] + penalty;

				/**********************************/
				diag_sim =  F[i-1][j-1] + similarity [seq1[i-1]][seq2[j-1]];

				/**********************************/
				max_sim         = diag_sim;
				direction[i][j] = 'd';
				if ( i_sim > max_sim ){
					max_sim = i_sim;
					direction[i][j] = 'i';
				}
				if ( j_sim > max_sim ) {
					max_sim = j_sim;
					direction[i][j] = 'j';
				}



			/* i is 0, j is not */
			} else if (j) {

				penalty =  0;
				if (j_between_exons) {
					penalty = 0;
						} else {
					if (use_endgap) {
					penalty = endgap;
							} else {
					if ( direction[i][j-1] =='j' ) {
						penalty = gap_extension;
								} else {
						penalty = gap_opening;
					}
					}
				}
				j_sim   = F[i][j-1] + penalty;
				max_sim = j_sim;
				direction[i][j] = 'j';

			/* j is 0, i is not */
			} else if (i) {

				penalty =  0;
				if (i_between_exons) {
					penalty = 0;
				} else {
					if ( use_endgap) {
						penalty = endgap;
					} else {
						if ( direction[i-1][j] == 'i' ) {
							penalty =  gap_extension;
						} else {
							penalty =  gap_opening;
						}
					}
				}
				i_sim   = F[i-1][j] + penalty;
				max_sim = i_sim;
				direction[i][j] = 'i';
			}

			if (max_sim < 0.0 ) max_sim = 0.0;

			F[i][j] = max_sim;
			if ( F_max < max_sim ) {
				// TODO{ tie break here */
				F_max = max_sim;
				F_max_i = i;
				F_max_j = j;
			}

		} /*end loop over j*/
	} /*end loop over i*/
		
		
	 
	i = F_max_i;
	j = F_max_j;
	// aln_score = F[i][j] ;


	while ( i>0 || j >0 ){

		if ( i<0 || j<0 ){
			sprintf (retstr, "Retracing error");
			return Py_BuildValue("s", retstr);
		}

		if (direction[i][j] == 'd'){
			map_i2j [i-1] = j-1;
			map_j2i [j-1] = i-1;
			i-= 1;
			j-= 1;
		} else if (direction[i][j] == 'i') {
			map_i2j [i-1] = far_away;
			i-= 1 ;

		} else if (direction[i][j] == 'j') {
			map_j2i [j-1] = far_away;
			j-= 1 ;

		} else{
			sprintf (retstr, "Retracing error");
			return Py_BuildValue("s", retstr);
		}
	}
	
	char * aligned_seq_1 = NULL;
	char * aligned_seq_2 = NULL;

	/* (lets hope it gets properly freed in the main program */
	if (! (aligned_seq_1 = emalloc( (len1+len2)*sizeof(char))) ) {
		sprintf (retstr, "error alloc array space");
		return Py_BuildValue("s", retstr);
	}
	if (! (aligned_seq_2 = emalloc( (len1+len2)*sizeof(char))) ) {
		sprintf (retstr, "error alloc array space");
		return Py_BuildValue("s", retstr);
	}

	i = 0;
	j = 0;
	int done = 0;
	int pos  = 0;
	while (!done) {

		if (j>=max_j && i>=max_i){
			done = 1;
		} else if (j<max_j && i<max_i){

			if (map_i2j[i] == j){
				aligned_seq_1[pos] = seq1[i];
				aligned_seq_2[pos] = seq2[j];
				i += 1;
				j += 1;
			} else if (map_i2j[i] < 0){
				aligned_seq_1[pos] = seq1[i];
				aligned_seq_2[pos] = gap_character;
				i += 1;
			} else if (map_j2i[j] < 0){
				aligned_seq_1[pos] = gap_character;
				aligned_seq_2[pos] = seq2[j];
				j += 1;
			}

		} else if (j<max_j){
			aligned_seq_1[pos] = gap_character;
			aligned_seq_2[pos] = seq2[j];
			j += 1;
		} else {
			aligned_seq_1[pos] = seq1[i];
			aligned_seq_2[pos] = gap_character;
			i += 1;
		}
		pos ++;
	}

	free_imatrix(F);
	free_cmatrix(direction);
	free(map_i2j);
	free(map_j2i);

	return Py_BuildValue("ss", aligned_seq_1, aligned_seq_2 );

}


/*
 * Bind Python function names to our C functions
 * this will be called as alignment.smith_waterman_context() from python
 */
static PyMethodDef alignment_methods[] = {
	{"smith_waterman_context", smith_waterman_context, METH_VARARGS, "exon-aware alignment"},
	{NULL, NULL, 0, NULL}        /* Sentinel, whatever that is */
};


/*
 * The method table must be referenced in the module definition structure: (this is python3)
 */
static struct PyModuleDef alnmtmodule =
{
    PyModuleDef_HEAD_INIT,
    "alignment", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    alignment_methods /* module methods come here */
};


/*
 * Python3 calls this to let us initialize our module
 * The initialization function must be named PyInit_name(), where name is the name of the module
 */
PyMODINIT_FUNC
PyInit_alignment(void)
{
    return PyModule_Create(&alnmtmodule);
}





