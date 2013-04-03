#!/usr/bin/python

import MySQLdb
from   mysql import search_db, switch_to_db
from   exon  import Exon

#########################################
def  get_analysis_dict(cursor):
    source = {}
    qry  = "select analysis_id, logic_name  from analysis"
    rows = search_db (cursor, qry)
    if (not rows):
        print "blah?"
        return False
    for row in rows:
        source[row[0]] = row[1]
    return source

#########################################
def  exon_id2gene_id (cursor, ensembl_db_name, exon_id, is_known):

    switch_to_db(cursor, ensembl_db_name)
    if is_known==2: # sw_sharp exon
        qry  = "select gene_id from sw_exon where "
        qry += "exon_id = %d " % exon_id
    else:
        qry  = "select gene_id from gene2exon where "
        qry += "exon_id = %d and is_known = %d " % (exon_id, is_known)
    
    rows = search_db (cursor, qry)
    if (not rows or 'ERROR' in rows[0]):
        rows = search_db (cursor, qry, verbose = True)
        return 0

    return rows[0][0]
    
#########################################
def  get_paras (cursor, gene_id):

    paras = []
    qry  = "select cognate_gene_id from paralogue "
    qry += " where gene_id = %d " % gene_id
    
    rows = search_db (cursor, qry)
    if (not rows):
        return []

    for row in rows:
        para_gene_id   = row[0]
        paras.append(para_gene_id)
    
    return paras

#########################################
def  get_orthos (cursor, gene_id, table):

    orthos = []
    qry  = "select cognate_gene_id, cognate_genome_db_id from "+table
    qry += " where gene_id = %d"% gene_id
    
    rows = search_db (cursor, qry)
    if (not rows):
        return []

    for row in rows:
        ortho_gene_id   = row[0]
        ortho_genome_db = row[1]
        # note: the cursor will be pointing to compara db after this
        species = genome_db_id2species (cursor, ortho_genome_db)
        orthos.append([ortho_gene_id,species] )
    
    return orthos

#########################################
def get_description (cursor, gene_id, db_name = None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False
    qry  = "select description from gene where gene_id = %d " % gene_id
    rows = search_db(cursor, qry)
    if rows:
        return rows[0][0]

    return ""

#########################################
def get_logic_name(cursor, analysis_id, db_name = None):
    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False
    qry = "SELECT logic_name FROM analysis WHERE analysis_id = %d" % analysis_id
    rows    = search_db (cursor, qry)
    if (not rows):
        logic_name = ''
    else:
        logic_name = rows[0][0]
    return logic_name 

#########################################
def is_coding (cursor, exon_id, db_name=None):
    
    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    qry = "select is_coding from gene where gene_id = %d " % int(exon_id)
    rows = search_db (cursor, qry)
    if ( not rows):
        return False

    return rows[0][0]>0
    
#########################################
def get_biotype (cursor, exon_id, db_name=None):
    
    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    qry = "select biotype from gene where gene_id = %d " % int(exon_id)
    rows = search_db (cursor, qry)
    if ( not rows):
        return False

    return rows[0][0]
    
#########################################
def get_status (cursor, exon_id, db_name=None):
    
    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    qry = "select status from gene where gene_id = %d " % int(exon_id)
    rows = search_db (cursor, qry)
    if ( not rows):
        return False

    return rows[0][0]
    
#########################################
def is_coding_exon (cursor, exon_id, is_known, db_name=None):
    
    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False
    qry = "select is_coding from gene2exon where exon_id = %d and is_known = %d" % (exon_id, is_known)
    rows = search_db (cursor, qry)
    if ( not rows):
        return False

    return rows[0][0]>0

#########################################
def get_gene_coordinates (cursor, gene_id, db_name=None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return None

    qry  = "select seq_region_id, seq_region_start, seq_region_end, seq_region_strand  "
    qry += " from gene "
    qry += " where gene_id = %d" %  gene_id
    rows = search_db (cursor, qry)


    if ( not rows or  isinstance(rows[0], str) and 'error' in rows[0].lower()):
         search_db (cursor, qry, verbose = True)
         return None

    return rows[0]

#########################################
def is_mitochondrial (cursor, gene_id, db_name=None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    qry  = "select seq_region.name from seq_region, gene  "
    qry += " where seq_region.seq_region_id =  gene.seq_region_id "
    qry += " and gene.gene_id = %d" %  gene_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         return []
    seq_name = rows[0][0]

    is_mitochondrial = (seq_name == 'MT')

    return is_mitochondrial

#########################################
def get_exon_pepseq (cursor, exon, db_name=None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    if  exon.analysis_id > 0:
        exon_id  = exon.exon_id
        is_known = exon.is_known
        qry  = "select protein_seq  "
        qry += " from exon_seq where exon_id = %d and is_known = %d" % (exon_id, is_known)
    else:
        exon_seq_id = exon.exon_seq_id
        qry  = "select protein_seq "
        qry += " from  exon_seq where exon_seq_id = %d" % exon_seq_id
        
    rows = search_db(cursor, qry)
    if (not rows):
        #rows = search_db(cursor, qry, verbose = True)
        return []

    protein_seq = rows[0][0]
    if (protein_seq is None):
        protein_seq = ""
  
    return protein_seq

#########################################
def get_exon_dnaseq (cursor, exon, db_name=None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return None

    if  exon.analysis_id > 0:
        exon_id  = exon.exon_id
        is_known = exon.is_known
        qry  = "select dna_seq  "
        qry += " from exon_seq where exon_id = %d and is_known = %d" % (exon_id, is_known)
    else:
        exon_seq_id = exon.exon_seq_id
        qry  = "select dna_seq "
        qry += " from  exon_seq where exon_seq_id = %d" % exon_seq_id
        
    rows = search_db(cursor, qry)
    if (not rows):
        #rows = search_db(cursor, qry, verbose = True)
        return None

    dna_seq = rows[0][0]
  
    return dna_seq

#########################################
def  get_sw_seq_id (cursor, exon_id, db_name=None):
    
    if (db_name):
        if not switch_to_db(cursor, db_name):
            return -1
    qry  = "select exon_seq_id "
    qry += " from sw_exon where exon_id = %d" % exon_id
    rows = search_db(cursor, qry)
    if not rows or not rows[0][0]:
        return -1
    
    return int(rows[0][0])

#########################################
def get_exon_seq_by_db_id (cursor, exon_seq_id, db_name=None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    qry  = "select exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
    qry += " left_flank, right_flank, dna_seq  "
    qry += " from exon_seq where exon_seq_id = %d" % exon_seq_id
    rows = search_db(cursor, qry)
    if (not rows):
        #rows = search_db(cursor, qry, verbose = True)
        return []

    [exon_seq_id, protein_seq, pepseq_transl_start, 
     pepseq_transl_end, left_flank, right_flank, dna_seq] = rows[0]
    if (protein_seq is None):
        protein_seq = ""
    if (left_flank  is None):
        left_flank = ""
    if (right_flank is None):
        right_flank = ""
    if (dna_seq is None):
        dna_seq = ""
    
    return [exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, left_flank, right_flank, dna_seq]

#########################################
def get_exon_phase (cursor, exon_id, is_known, db_name=None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    if is_known==2: # sw exon
        qry  = "select phase from sw_exon"
    elif is_known==1:
        qry  = "select phase from exon"
    elif is_known==0:
        qry  = "select start_phase from prediction_exon"
    else:
        return None

    qry += " where exon_id = %d " %  exon_id
    rows = search_db(cursor, qry)

    if not rows:
        rows = search_db(cursor, qry, verbose = True)
        return None

    
    return  rows[0][0]

#########################################
def get_exon_seqs (cursor, exon_id, is_known, db_name=None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    if is_known==2: # sw exon
        qry  = "select exon_seq.exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
        qry += " left_flank, right_flank, dna_seq  from  exon_seq "
        qry += " join sw_exon on exon_seq.exon_seq_id = sw_exon.exon_seq_id  "
        qry += " where sw_exon.exon_id = %d " %  exon_id
    else:
        qry  = "select exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
        qry += " left_flank, right_flank, dna_seq  "
        qry += " from  exon_seq where exon_id = %d and is_known = %d" % (exon_id, is_known)

    rows = search_db(cursor, qry)

    if not rows or len(rows[0]) < 7:
        print "using db ", db_name
        rows = search_db(cursor, qry, verbose = True)
        return []

    [exon_seq_id, protein_seq, pepseq_transl_start, 
     pepseq_transl_end, left_flank, right_flank, dna_seq] = rows[0]
    if (protein_seq is None):
        protein_seq = ""
    if (left_flank is None):
        left_flank = ""
    if (right_flank is None):
        right_flank = ""
    if (dna_seq is None):
        dna_seq = ""
    
    return [exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, left_flank, right_flank, dna_seq]


#########################################
def get_exon (cursor, exon_id, is_known=None, db_name=None):

    exon = Exon ()

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return exon

    if is_known and is_known==2:
        # sw# exon
        qry  = "select * from sw_exon where exon_id = %d"   % exon_id
        rows = search_db(cursor, qry, verbose=False)
        if (not rows):
            return exon
        exon.load_from_sw_exon (rows[0])
    else:
        qry  = "select * from gene2exon where exon_id = %d" % exon_id
        if is_known: qry += " and is_known = %s " % is_known
        rows = search_db(cursor, qry, verbose=False)
        if (not rows):
            return exon
        exon.load_from_gene2exon (rows[0])

    return exon


#########################################
def gene2exon_list (cursor, gene_id, db_name=None):

    exons = []

    if (db_name): 
        if not switch_to_db(cursor, db_name):
            return False

    qry  = "select * from gene2exon where gene_id = %d " % gene_id
    rows = search_db(cursor, qry)
    if (not rows):
        rows = search_db(cursor, 'select database()')
        print "database ", rows[0][0]
        rows = search_db(cursor, qry, verbose = True)
        return []

    for row in rows:
        exon = Exon()
        if (not exon.load_from_gene2exon(row)):
            continue
        exons.append(exon)

    return exons

########################################
def get_canonical_exons (cursor, gene_id):

    exons = gene2exon_list (cursor, gene_id)
    if (not exons):
        print gene2stable (cursor, gene_id = gene_id), " no exons found ", ct, tot
        exit(1) # shouldn't happen at this point

    # sorting exons in place by their start in gene:
    exons.sort(key=lambda exon: exon.start_in_gene)

    canonical_coding_exons = []
    reading = False
    for exon in exons:
        if (not exon.is_canonical): 
             continue
        if (not exon.canon_transl_start is None):
            reading = True
        if (reading):
            canonical_coding_exons.append(exon)
        if (not exon.canon_transl_end is None):
            break  

    return canonical_coding_exons

#########################################
def get_selenocysteines (cursor, gene_id):

    selenoC_pos = []

    canonical_transl_id = gene2canon_transl(cursor, gene_id)

    qry  = "select value from translation_attrib "
    qry += " where attrib_type_id = 12 and  translation_id = %d " % canonical_transl_id

    rows = search_db (cursor, qry)
    if (not rows):
       return []

    for row in rows:
        blah  = row[0].split (" ")
        start = int(blah[0])
        end   = int(blah[1])
        for pos in range(start, end+1):
            selenoC_pos.append(pos-1)

    return selenoC_pos


#########################################
def gene2stable_canon_transl(cursor, gene_id, db_name=None):

    if (db_name and not switch_to_db(cursor, db_name)):
            return False
 
    qry  = "select translation.stable_id  from translation, gene "
    qry += " where gene.canonical_transcript_id = translation.transcript_id "
    qry += " and gene.gene_id = %d " % gene_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return  rows[0][0]

#########################################
def gene2canon_transl(cursor, gene_id, db_name=None,):

    if  (db_name and not switch_to_db(cursor, db_name)):
            return False
 
    qry  = "select translation.translation_id  from translation, gene "
    qry += " where gene.canonical_transcript_id = translation.transcript_id "
    qry += " and gene.gene_id = %d " % gene_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return  rows[0][0]




########
def stable2gene (cursor, stable_id=None, db_name=None):

    if (not stable_id):
        return 0

    if (db_name and not switch_to_db(cursor, db_name)):
            return False

    qry = "select gene_id from gene where stable_id='%s'" % stable_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return 0

    return int(rows[0][0])
    
########
def gene2stable (cursor, gene_id=None, db_name=None):

    if (not gene_id):
        return ""

    if  (db_name):
        qry  = "use %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            rows = search_db (cursor, qry, verbose = True)
            print rows
            exit (1)


    qry = "select stable_id from gene where gene_id=%d" % gene_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return rows[0][0]
    
########
def exon2stable (cursor, exon_id=None, db_name=None):

    if (not exon_id):
        return ""

    if  (db_name):
        qry  = "use %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            rows = search_db (cursor, qry, verbose = True)
            print rows
            exit (1)


    qry = "select stable_id from exon where exon_id=%d" % exon_id
    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return rows[0][0]
    

########
def get_gene_ids (cursor, db_name=None, biotype = None, is_known = None):

    gene_ids = []
    
    if  (db_name):
        qry  = "use %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            rows = search_db (cursor, qry, verbose = True)
            print rows
            exit (1)

    qry = "select gene_id from gene"

    if ( biotype or is_known):
        qry +=  " where "
        if ( biotype):
           qry += "biotype='%s'" % biotype
        if (biotype and is_known):
            qry += " and "
        if (is_known):
           qry += "status='known'"

    rows = search_db (cursor, qry, verbose = False)
    
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return []
    else:
        if ('Error' in rows[0]):
            print rows[0]
            return []

        for row in rows:
            if ( not type(row[0]) is long ):
                print row
                exit(1)
            gene_ids.append(int(row[0]))
    
    return gene_ids




########
def get_species (cursor):

    ensembl_db_name = {}
    all_species     = []

    qry = "show databases like '%core%'"
    cursor.execute(qry)

    rows = cursor.fetchall()
    if (not rows):
        print "No databases with 'core' in the name found"
        return 1

    for row in rows:
        db_name    = row[0]
        name_token = db_name.split ('_')
        species = name_token[0]
        i = 1
        while not name_token[i] == 'core':
            species += "_"+ name_token[i]
            i       += 1
        ensembl_db_name[species] = db_name
        all_species.append(species)

    return all_species, ensembl_db_name

########
def get_compara_name (cursor):

    qry = "show databases like '%compara%'"
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return rows[0][0]

########
def species2taxid (cursor, species):

    switch_to_db (cursor, get_compara_name (cursor))
    qry  = "select taxon_id from genome_db where name = '%s'" % species
    rows = search_db (cursor, qry)
    if (not rows):
        search_db (cursor, qry, verbose = True)
        return ""
    
    return rows[0][0]

########
def species2genome_db_id (cursor, species):


    switch_to_db (cursor, get_compara_name (cursor))

    qry  = "select genome_db_id from genome_db where name = '%s'" % species

    rows = search_db (cursor, qry)
    if (not rows):
        search_db (cursor, qry, verbose = True)
        return 0
    
    return int(rows[0][0])

########
def genome_db_id2species (cursor, genome_db_id):


    switch_to_db (cursor, get_compara_name (cursor))

    qry  = "select name from genome_db where genome_db_id = %d" % int(genome_db_id)

    rows = search_db (cursor, qry)
    if (not rows):
        search_db (cursor, qry, verbose = True)
        return ""
    
    return rows[0][0]
