#!/usr/bin/python

import MySQLdb
from   mysql import search_db, switch_to_db
from   exon  import Exon
import commands

#########################################
def get_canonical_transcript_id (cursor, gene_id, db_name=None):

    if db_name and not switch_to_db(cursor, db_name):
        return []

    qry     = "select canonical_transcript_id"
    qry    += " from  gene where gene_id=%d"  %  gene_id
    rows    = search_db (cursor, qry, verbose=False)

    if (not rows):
        rows    = search_db (cursor, qry, verbose=True)
        return ""
    elif ( 'Error' in rows[0]):
        print  rows[0]
        return ""

    return rows[0][0]

#########################################
def get_canonical_coding_exons (cursor, gene_id, db_name=None):

    if db_name and not switch_to_db(cursor, db_name):
        return []

    exons_all =  gene2exon_list (cursor, gene_id)
    if not exons_all:  return []

    exons = filter (lambda x: x.is_coding and x.is_canonical, all_exons)
    # now, the problem is that an exon can be coding, 
    # but not in the canonical version of the transcript
    exons = exons.sort(key=lambda exon: exon.start_in_gene)
    # is this gene on the forward or on the reverse strand?
    ret  = get_gene_region (cursor, gene_id)
    if not ret:  return []
    [seq_region_id, seq_region_start, seq_region_end, strand] = ret
    if strand < 0: exons.reverse()
    # is there info about the beginning and the end of canonical translation?
    canonical_transcript_id  = get_canonical_transcript_id (cursor, gene_id, db_name=None)
    if not canonical_transcript_id: return []
    ret = get_canonical_coordinates (cursor, canonical_transcript_id)
    if not ret or not len(ret) == 4: return []

    [canonical_start_in_exon, canonical_start_exon_id,
     canonical_end_in_exon, canonical_end_exon_id] = ret
    if canonical_start_exon_id is None or  canonical_end_exon_id is None: return []
    
    # filter the exons that are within the start and end bracket
    canonical_exons = []
    reading = 0
    for exon in exons:
        if exon.exon_id == canonical_start_exon_id:  reading = 1
        if reading: canonical_exons.append(exon)
        if exon.exon_id == canonical_end_exon_id:   break

    if strand < 0: canonical_exons.reverse()

    return canonical_exons



#########################################
def get_gene_region (cursor, gene_id, is_known=None):

    qry     = "select seq_region_id, seq_region_start, seq_region_end, "
    qry    += " seq_region_strand "
    qry    += " from  gene  where  gene_id=%d"  %  gene_id
    if (not is_known is None and is_known):
        qry  += " and  status='known' "
    rows    = search_db (cursor, qry, verbose=False)

    if (not rows):
        rows = search_db (cursor, qry, verbose=True)
        return []
    elif ( 'Error' in rows[0]):
        print  rows[0]
        return []

    return rows[0]

#########################################
def get_known_exons (cursor, gene_id, species):

    exons = []

    qry  = "select distinct exon_transcript.exon_id from  exon_transcript, transcript "
    qry += " where exon_transcript.transcript_id = transcript.transcript_id "
    qry += " and transcript.gene_id = %d " % gene_id

    rows = search_db (cursor, qry)
    
    if (not rows ):
        return []
    if ('Error' in rows[0]):
        search_db (cursor, qry, verbose = True)
        return []

    # get the region on the gene
    ret = get_gene_region (cursor, gene_id)
    if  ret:
        [gene_seq_id, gene_region_start, gene_region_end, 
         gene_region_strand] = ret
    else:
        print "region not retrived for ", species, gene_id
        return []

    exon_ids = []
    for row in rows:
        exon_ids.append(row[0])

    for exon_id in exon_ids:
        qry = "select * from exon where exon_id=%d" % exon_id
        rows = search_db (cursor, qry)
        if (not rows or 'Error' in rows[0]):
            search_db (cursor, qry, verbose = True)
            continue
        exon         = Exon()
        exon.gene_id = gene_id
        exon.load_from_ensembl_exon (gene_region_start, gene_region_end, rows[0])
        exons.append(exon)

    return exons


#########################################
def get_predicted_exons (cursor, gene_id, species):

    exons = []

    # get the region on the gene
    ret = get_gene_region (cursor, gene_id)
    if  ret:
        [gene_seq_id, gene_region_start, gene_region_end, 
         gene_region_strand] = ret
    else:
        print "region not retrived for ", species, gene_id
        return []

    qry    = "SELECT  * FROM  prediction_exon  WHERE seq_region_id = %d "  %  gene_seq_id
    qry   += " AND  seq_region_start >= %d AND seq_region_start <= %d " %  \
        (gene_region_start, gene_region_end)
    qry   += " AND  seq_region_end   >= %d AND seq_region_end   <= %d " %  \
        (gene_region_start, gene_region_end)
    rows   = search_db (cursor, qry)

    if (not rows):
        return []
    for row in rows:
        exon         = Exon()
        exon.gene_id = gene_id
        exon.load_from_ensembl_prediction (gene_region_start, gene_region_end, row)
        exons.append(exon)
 
    return exons

#########################################
def get_novel_exons (cursor, gene_id, table):

    exons = []

    qry  = "select * from %s " % table
    qry += " where gene_id = %d " % int(gene_id)
    rows = search_db (cursor, qry)
    if not rows: return exons

    for row in rows:
        exon         = Exon()
        exon.load_from_novel_exon (row, table)
        exons.append(exon)
    return exons



#########################################
def  get_primary_seq_info (cursor, gene_id, species):

    # seq identifier from gene table
    qry  = "select seq_region_id, seq_region_start, seq_region_end, seq_region_strand"
    qry += " from gene where gene_id = %d" % gene_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         exit(1)
    [seq_region_id, seq_region_start, seq_region_end, seq_region_strand] = rows[0]
    

    qry  = "select name, file_name from seq_region "
    qry += " where seq_region_id= %d" %  seq_region_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         return []
    [seq_name, file_names] = rows[0]
    # Ivana:
    # What is the proper way to find out whether the seq_region is mitochondrial?
    # EMily:
    # It's in the seq_region table, under name. Name will either be a chromosome
    # number, MT for mitochrondria or the name of a contig.
    mitochondrial = is_mitochondrial(cursor, gene_id)

    return [seq_name, file_names, seq_region_start, 
            seq_region_end, seq_region_strand, mitochondrial]

#########################################
def  get_alt_seq_info (cursor, gene_id, species):

    # seq identifier from gene table
    qry  = "select seq_region_id, seq_region_strand from gene where gene_id = %d" % gene_id
    rows = search_db (cursor, qry)
    if ( not rows):
         search_db (cursor, qry, verbose = True)
         exit(1)
    [seq_region_id, seq_region_strand] = rows[0]
    
    # check whether we have "assembly exception"
    # we do not want 'PAR' regions, though:
    '''
    The pseudo-autosomal regions are homologous DNA sequences on the (human) X and Y chromosomes. 
    They allow the pairing and crossing-over of these sex chromosomes the same way the autosomal 
    chromosomes do during meiosis. 
    As these genomic regions are identical between X and Y, they are oftentimes only stored once.
    '''

    qry  = "select seq_region.name,  assembly_exception.exc_seq_region_start, assembly_exception.exc_seq_region_end "
    qry += "from seq_region, assembly_exception "
    qry += "where seq_region.seq_region_id = assembly_exception.exc_seq_region_id "
    qry += "and assembly_exception.seq_region_id = %d" % seq_region_id
    qry += " and not assembly_exception.exc_type = 'PAR'"
    rows = search_db (cursor, qry)
    if (rows):
        [seq_name, seq_region_start, seq_region_end] = rows[0]
        qry = " select distinct file_name from seq_region where seq_region.name = '%s' " % seq_name
        rows = search_db (cursor, qry)
        file_names = ""
        for row in rows:
            if file_names:
                file_names += " "
            file_names += row[0]

        mitochondrial = is_mitochondrial (cursor, gene_id)
        return [seq_name, file_names, seq_region_start, 
                seq_region_end, seq_region_strand, mitochondrial]
    else:
        return []

#########################################
def extract_gene_seq (acg, species, seq_name, file_names, seq_region_strand,  
                      seq_region_start, seq_region_end):

    # now the question is, which file do I use if there are several options?
    first_choice  = ""
    second_choice = ""
    for file_name in file_names.split(" "):
        if  '.chromosome.' in  file_name:
            first_choice = file_name
        elif '.toplevel.' in  file_name:
            second_choice = file_name
            
    fasta_db_file = ""
    if first_choice:
        fasta_db_file = first_choice
    elif second_choice:
        fasta_db_file = second_choice

    if not fasta_db_file:
        print "failed to decide on fasta_db_file:"
        print file_names
        exit(1)


    # extract gene sequence  from fasta db
    fastacmd = acg.generate_fastacmd_gene_command(species, seq_name, fasta_db_file,
                                                  seq_region_strand,  seq_region_start,    
                                                  seq_region_end)
    ret = commands.getoutput(fastacmd)
    if not ret:
        print "no refturn for fastacmd for", species, gene_id
        print "fastacmd: ", fastacmd
        exit (1)
    if ('ERROR' in ret):
        print "Error running fastacmd: ", fastacmd
        print ret
        if 'Ignoring sequence location' in ret:
            print 'will ignore'
            print
        else:
            exit (1)
    gene_seq = ""
    reading = 0
    for line in ret.split("\n"):
        if ('>' in line):
            reading = 1
            continue
        if (not reading):
            continue
        line.rstrip()
        gene_seq += line

    return gene_seq



#########################################
def get_transcript_ids (cursor, gene_id, db_name=None):

    transcript_ids = []

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return []

    qry  = "select transcript_id, stable_id from transcript where gene_id = %d " % gene_id
    rows = search_db(cursor, qry)

    if not rows:
        return transcript_ids
    
    for row in rows:
        transcript_ids.append(row)

    return transcript_ids


#########################################
def get_translation_coords (cursor, transcript_id, db_name=None):
    
    if (db_name):
        if not switch_to_db(cursor, db_name):
            return []

    qry  = "select  seq_start, start_exon_id, seq_end, end_exon_id  "
    qry += " from translation where transcript_id = %d " % transcript_id
   
    rows = search_db(cursor, qry)

    if not rows:
        return []

    if 'error' in str(rows[0][0]).lower():
        return []

    return rows[0]

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
    if is_known==3: # sw_sharp exon
        qry  = "select gene_id from usearch_exon where "
        qry += "exon_id = %d " % exon_id
    elif is_known==2: # sw_sharp exon
        qry  = "select gene_id from sw_exon where "
        qry += "exon_id = %d " % exon_id
    else:
        qry  = "select gene_id from gene2exon where "
        qry += "exon_id = %d and is_known = %d " % (exon_id, is_known)
    
    rows = search_db (cursor, qry)
    if (not rows or 'ERROR' in rows[0]):
        rows = search_db (cursor, qry, verbose = True)
        return ""

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
    qry  = "select description from gene where gene_id = %d " % int(gene_id)
    rows = search_db(cursor, qry)
    if rows:
        return rows[0][0]

    return ""

#########################################
def get_logic_name(cursor, analysis_id, db_name = None):

    if analysis_id < 0:
        return ''

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
def is_coding (cursor, gene_id, db_name=None):
    
    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    qry = "select is_coding from gene where gene_id = %d " % int(gene_id)
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
def get_exon_pepseq (cursor, exon, db_name=None, verbose=False):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return False

    if exon.exon_seq_id:
        exon_seq_id = exon.exon_seq_id
        qry  = "select protein_seq "
        qry += " from exon_seq where exon_seq_id = %d" % exon_seq_id
    else:
        exon_id  = exon.exon_id
        is_known = exon.is_known
        qry  = "select protein_seq  "
        qry += " from exon_seq where exon_id = %d and is_known = %d" % (exon_id, is_known)

        
    rows = search_db(cursor, qry)


    if (not rows):
        if verbose:
            rows = search_db(cursor, qry, verbose = True)
            print rows
        return ""


    protein_seq = rows[0][0]
    if (protein_seq is None):
        protein_seq = ""
  
    return protein_seq

#########################################
def exon_seq_id2exon_id (cursor, exon_seq_id, db_name=None):

    if (db_name):
        if not switch_to_db(cursor, db_name):
            return ""
    qry  = "select exon_id, is_known from exon_seq where exon_seq_id = %d  " % int(exon_seq_id)
    rows = search_db(cursor, qry)
    if (not rows):
        return ""


    return rows[0]

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

    elif is_known==3: # usearch exon
        qry  = "select exon_seq.exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
        qry += " left_flank, right_flank, dna_seq  from  exon_seq "
        qry += " join usearch_exon on exon_seq.exon_seq_id = usearch_exon.exon_seq_id  "
        qry += " where usearch_exon.exon_id = %d " %  exon_id

    else:
        qry  = "select exon_seq_id, protein_seq, pepseq_transl_start, pepseq_transl_end, "
        qry += " left_flank, right_flank, dna_seq  "
        qry += " from  exon_seq where exon_id = %d and is_known = %d" % (exon_id, is_known)

    rows = search_db(cursor, qry)

    if not rows or not len(rows[0]) == 7:
        #rows = search_db(cursor, "select database()")
        #print "using db ", rows[0]
        #rows = search_db(cursor, qry, verbose = True)
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

    if is_known==2:
        # sw# exon
        qry  = "select * from sw_exon where exon_id = %d"   % exon_id
        rows = search_db(cursor, qry, verbose=False)
        if (not rows):
            return exon
        exon.load_from_novel_exon (rows[0], "sw_exon")
    elif is_known==3:
        # sw# exon
        qry  = "select * from usearch_exon where exon_id = %d"   % exon_id
        rows = search_db(cursor, qry, verbose=False)
        if (not rows):
            return exon
        exon.load_from_novel_exon (rows[0], "usearch_exon")
    else:
        qry  = "select * from gene2exon where exon_id = %d" % exon_id
        if is_known: qry += " and is_known = %s " % is_known
        rows = search_db(cursor, qry, verbose=False)
        if (not rows):
            return exon
        exon.load_from_gene2exon (rows[0])

    return exon


#########################################
def gene2exon_list (cursor, gene_id, db_name=None, verbose=False):

    exons = []

    if (db_name): 
        if not switch_to_db(cursor, db_name):
            return False

    qry  = "select * from gene2exon where gene_id = %d " % gene_id
    rows = search_db(cursor, qry)

    if (not rows):
        rows = search_db(cursor, 'select database()')
        if verbose:
            print "database ", rows[0][0]
            rows = search_db(cursor, qry, verbose = True)
            print rows
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
def stable2member (cursor, stable_id):
    
    # member_id refers to compara db
    # of which we need to have one
    qry = "select  member_id from member where stable_id = '%s'" % stable_id
    rows = search_db (cursor, qry)
    if (not rows or 'ERROR' in rows[0]):
        rows = search_db (cursor, qry, verbose = True)
        exit(1)
        return ""
    
    return int(rows[0][0])

########
def member2stable (cursor, member_id):
    
    # member_id refers to compara db
    # of which we need to have one
    qry = "select  stable_id from member where member_id = %d" % member_id
    rows = search_db (cursor, qry)
    if (not rows):
        rows = search_db (cursor, qry, verbose = True)
        return ""

    return rows[0][0]



########
def get_orthologues(cursor, ortho_type, member_id):

    # cursor must be pointing to compara database
    # the ortho pytpe i sone of the following: 'ortholog_one2one', 
    # 'ortholog_one2many', 'ortholog_many2many'

    orthos = []

    qry  = "select homology.homology_id from homology_member, homology "
    qry += " where homology_member.member_id =%d " % member_id
    qry += " and homology.homology_id = homology_member.homology_id "
    qry += " and  homology.description = '%s' " % ortho_type
    rows = search_db (cursor, qry)

    if (not rows):
        return [] # no orthologs here

    # for each homology id find the other member id
    for row in rows:
        homology_id = row[0]

        qry  = "select member_id from homology_member "
        qry += " where homology_id = %d"  % int(homology_id)
        qry += " and not  member_id = %d" % member_id
        rows2 = search_db (cursor, qry, verbose = True)
        if (not rows2):
            rows2 = search_db (cursor, qry, verbose = True)
            return []
        ortho_id     = rows2[0][0]
        
        qry  = "select  member.stable_id, genome_db.name, genome_db.genome_db_id "
        qry += " from member, genome_db "
        qry += " where member.member_id = %d " % ortho_id
        qry += " and genome_db.genome_db_id = member.genome_db_id"
        rows3 = search_db (cursor, qry, verbose = True)
        if (not rows3):
            rows3 = search_db (cursor, qry, verbose = True)
            return []
        [ortho_stable, species, genome_db_id]     = rows3[0]
        orthos.append([ortho_stable, species,  int(genome_db_id)])
        
    return orthos
        

########################
def get_orthologues_from_species(cursor, ensembl_db_name, ortho_type, member_id, species):

    # the ortho_type is one of the following: 'ortholog_one2one', 
    # 'ortholog_one2many', 'ortholog_many2many', 'possible_ortholog', 'apparent_ortholog_one2one'
    orthos = []

    # find genome db_id
    genome_db_id = species2genome_db_id (cursor, species)

    # make the cursor point to compara database - should be the responsibility of each function
    switch_to_db (cursor, get_compara_name (cursor))

    qry  = "select homology.homology_id from homology_member, homology "
    qry += " where homology_member.member_id =%d " % member_id
    qry += " and homology.homology_id = homology_member.homology_id "
    qry += " and  homology.description = '%s' "    % ortho_type
    rows = search_db (cursor, qry)

    if (not rows):
        return [] # no orthologs here

    # for each homology id find the other member id
    #print qry
    #print member_id, ortho_type, species, genome_db_id
    #print rows
    for row in rows:
        homology_id = row[0]
        #print "\t homology id:", homology_id
        switch_to_db (cursor, get_compara_name (cursor))
        qry  = "select member_id from homology_member "
        qry += " where homology_id = %d"  % int(homology_id)
        qry += " and not  member_id = %d" % member_id

        rows2  = search_db (cursor, qry, verbose = False)
        if (not rows2):
            #print "\t ",
            #rows2 = search_db (cursor, qry, verbose = True)
            continue
        for row2 in rows2:
            ortho_id     = row2[0]
            #print "\t\t ortho id:", ortho_id
            qry  = "select  stable_id  from member  "
            qry += " where member_id = %d "  % ortho_id
            qry += " and genome_db_id = %d " % genome_db_id
            rows3 = search_db (cursor, qry, verbose = False)
            if (not rows3):
                #print "\t\t ",
                #rows3 = search_db (cursor, qry, verbose = True)
                continue
            ortho_stable  = rows3[0][0]
            #print "\t\t ortho stable:", ortho_stable
            orthos.append(ortho_stable)
    if orthos:    
        switch_to_db (cursor, ensembl_db_name [species])
        orthos = map  (lambda gene_id:  stable2gene(cursor, gene_id), orthos)
    #print 'orthos:', orthos
    return orthos
        



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
def stable2exon (cursor, stable_id, db_name=None):


    if  (db_name):
        qry  = "use %s " % db_name
        rows = search_db (cursor, qry)
        if (rows):
            rows = search_db (cursor, qry, verbose = True)
            print rows
            exit (1)

    qry  = "select exon_id from exon where stable_id='%s'" % stable_id
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

    if ( biotype or not is_known is None):
        qry +=  " where "
        if ( biotype):
           qry += "biotype='%s'" % biotype
        if (biotype and not is_known is None):
            qry += " and "
        if (not is_known is None):
            if is_known:
                qry += "status='known'"
            else:
                qry += "status!='known'"
                

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
