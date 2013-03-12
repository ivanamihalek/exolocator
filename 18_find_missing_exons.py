#!/usr/bin/python

import MySQLdb, subprocess, re, commands
from tempfile     import NamedTemporaryFile
from bisect       import bisect
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna
from difflib      import SequenceMatcher
from   el_utils.mysql   import  *
from   el_utils.ensembl import  *
from   el_utils.utils   import  erropen, output_fasta
from   el_utils.map     import  Map, get_maps
from   el_utils.tree    import  species_sort
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader
from   el_utils.threads import  parallelize


#########################################
MAX_SEARCH_LENGTH = 100000

#########################################
class Seq_Region:
    def __init__(self, name, filename, strand, start, end):
        self.species = filename.split('.')[0].lower()
	self.name = name
	self.filename = filename
	self.strand = 1 if strand == 1 else 2
	self.start = start
	self.end = end
    def __str__(self):
        return "(" + ",".join((self.name,self.filename,str(self.strand),str(self.start),str(self.end))) + ")"

#########################################
def get_best_filename(names_string):
    names_list = names_string.split()
    chromo_names = [n for n in names_list if 'chromosom' in n]
    if chromo_names: return chromo_names[0]
    else: return names_list[0]

#########################################
def get_gene_start_end(cursor, ensembl_db_name, gene_id):
    switch_to_db (cursor, ensembl_db_name)
    qry  = "select seq_region_start, seq_region_end from gene where gene_id=%d" % gene_id
    rows = search_db(cursor, qry)
    if not rows: 
        return None
    return  map(int, rows[0])

#########################################
def get_seq_region(cursor, ensembl_db_name, exon_id, exon_known):

    switch_to_db (cursor, ensembl_db_name)
    qry  = "select name, file_name, seq_region_strand, seq_region_start, seq_region_end "
    if exon_known == 1:
        qry += "from exon "
    elif exon_known == 0:
        qry += "from prediction_exon "
    elif exon_known == 2:
        print "get_seq_region() from sw exon not implemented - need to finish it up from the gene region "
        return None
    else:
        print "error in get_seq_region() "
        return None
        
    qry += "join seq_region on exon.seq_region_id = seq_region.seq_region_id "
    qry += "where exon_id=%d" % exon_id
    rows = search_db(cursor, qry)
    if not rows:   return False
    if not len(rows[0]) ==5:
        print rows[0]
        return False
    name, file_names, strand, start, end = rows[0]
    return Seq_Region(name, get_best_filename(file_names), int(strand), int(start), int(end))

#########################################
def get_alt_seq_region(cursor,ensembl_db_name, exon_id, exon_known):

    switch_to_db (cursor, ensembl_db_name)
    qry  = "select name, file_name, seq_region_strand, exc_seq_region_start, exc_seq_region_end "
    if exon_known == 1:
        qry += "from exon "
    elif exon_known == 0:
        qry += "from prediction_exon "
    elif exon_known == 2:
        print "get_alt_seq_region() from sw exon not implemented - need to fishi it up from the gene region "
        return None
    else:
        print "error in get_alt_seq_region() "
        return None
    
    qry += "join assembly_exception as ae on ae.seq_region_id=exon.seq_region_id "
    qry += "join seq_region on seq_region.seq_region_id=ae.exc_seq_region_id "
    qry += "where exon_id=%d" % exon_id
    rows = search_db(cursor, qry)
    if not rows: return False
    name, file_names, strand, start, end = rows[0]
    return Seq_Region(name, get_best_filename(file_names), int(strand), int(start), int(end))

#########################################
def find_sw_exon_id (cursor, human_exon_id):

    query  = "select exon_id  from sw_exon "
    query += " where maps_to_human_exon_id = '%s' " % human_exon_id
    rows   = search_db(cursor, query)
    if not rows or  'Error' in rows[0]:
        return ''
    else:
        return int(rows[0][0])

#########################################
def find_exon_seq_id (cursor, sw_exon_id):

    query  = "select exon_seq_id  from exon_seq "
    query += " where exon_id = %d  and is_sw = 1 " % sw_exon_id
    rows   = search_db(cursor, query)
    if not rows or  'Error' in rows[0]:
        return ''
    else:
        return int (rows[0][0])

#########################################
def store_sw_exon (cursor, db_name, human_exon_id, gene_id, start_in_gene, 
		   end_in_gene,  strand,  dnaseq, left_flank, right_flank, pepseq,
                   has_NNN, has_stop, template_exon_id, template_species):

    switch_to_db(cursor, db_name)

    #################################
    fixed_fields  = {}
    update_fields = {}

    fixed_fields['maps_to_human_exon_id'] = human_exon_id

    update_fields['gene_id']          = gene_id
    update_fields['start_in_gene']    = start_in_gene
    update_fields['end_in_gene']      = end_in_gene
    update_fields['strand']           = strand
    update_fields['has_NNN']          = has_NNN
    update_fields['has_stop']         = has_stop

    update_fields['template_exon_id'] = int(template_exon_id)
    update_fields['template_species'] = template_species

    store_or_update (cursor, 'sw_exon', fixed_fields, update_fields)

    #################################
    sw_exon_id = find_sw_exon_id (cursor, human_exon_id)
    if not sw_exon_id: return False

    #################################
    fixed_fields  = {}
    update_fields = {}

    fixed_fields ['exon_id']             = sw_exon_id
    fixed_fields ['is_sw']               = 1

    update_fields['is_known']            = 0
    update_fields['dna_seq']             = dnaseq
    update_fields['left_flank']          = left_flank
    update_fields['right_flank']         = right_flank
    update_fields['protein_seq']         = pepseq
    update_fields['pepseq_transl_start'] = 0
    update_fields['pepseq_transl_end']   = len(dnaseq)

    store_or_update (cursor, 'exon_seq', fixed_fields, update_fields)

    #################################
    exon_seq_id = find_exon_seq_id (cursor, sw_exon_id)
    if not exon_seq_id: return False


    #################################
    qry  = "update sw_exon set exon_seq_id = {0} where exon_id = {1}".format(exon_seq_id, sw_exon_id)
    cursor.execute(qry)

    return sw_exon_id

#########################################
def translate(searchseq, search_start, search_end, mitochondrial):
    dnaseq = searchseq[search_start:search_end]
    if mitochondrial:
        try:
            pepseq = str(Seq(dnaseq,generic_dna).translate(table="Vertebrate Mitochondrial"))
        except:
            pepseq = ""
    else:
        try:
            pepseq = str(Seq(dnaseq,generic_dna).translate())
        except:
            pepseq = ""
  
    return pepseq

#########################################
def find_the_most_similar_frame (cursor, mitochondrial, template_pepseq, searchseq,  search_start, search_end):

    pepseq = translate(searchseq, search_start, search_end, mitochondrial)
    if pepseq==template_pepseq:
        return  [pepseq, search_start, search_end]

    s = SequenceMatcher(None, pepseq, template_pepseq)
    if s.ratio() > 0.5 and not '*' in pepseq:
        return  [pepseq, search_start, search_end]

    ratio = {}
    for b in range(-2, 3): # begin offset
        if search_start+b < 0: continue
        if search_start+b > len(searchseq): continue
        for e in range(-2, 3): # end offset
            if search_end + e < 0: continue
            if search_end + e > len(searchseq): continue
            
            pepseq = translate(searchseq, search_start+b, search_end+e,  mitochondrial)
            s = SequenceMatcher(None, pepseq, template_pepseq)
            key = (b,e)
            ratio[key] =  s.ratio()
    
    sorted_keys = sorted(ratio, key=lambda key: -ratio[key])
    
    key = sorted_keys[0]
    (b,e) = key
    if ratio[key] < 0.5:
        return ["", search_start, search_end]

    pepseq = translate(searchseq, search_start+b, search_end+e,  mitochondrial)

    if "*" in pepseq and len(pepseq) > len(template_pepseq): # if the stop is outside of the matching region with template, 
                                                      # return only the matching region
        print 'resolving a stop codon'
        s      = SequenceMatcher(None, pepseq, template_pepseq)
        blocks = s.get_matching_blocks()
        blocks = filter( lambda (i,j,n): n>0, blocks)
        if len(blocks)==1:
            (i, j, n) = blocks[0]
            if n == len(template_pepseq):
                b += 3*i
                e = search_start + b + 3*n - search_end  
                pepseq = pepseq[i:i+n]
        else:
            (i_b,j_b,n_b) = blocks[0]
            (i_e,j_e,n_e) = blocks[-1]
            matching_region_length = i_e+n_e-i_b
            if matching_region_length == len(template_pepseq):
                b += 3*i_b
                e = search_start + b + 3*(i_e+n_e-i_b)  - search_end
                pepseq = pepseq[i_b:i_e+n_e]
            

    return  [pepseq, search_start+b, search_end+e]
        

#########################################
def  check_sw_exon_exists(cursor, ensembl_db_name, species, maps_to_human_exon_id):
    switch_to_db (cursor, ensembl_db_name[species])
    qry  = "select count(1) from sw_exon where maps_to_human_exon_id=%d" % maps_to_human_exon_id
    rows = search_db(cursor, qry)
    if not rows or 'ERROR' in rows[0] or rows[0][0] == 0:
        return False
    else:
        return True
    

#########################################
def sw_search(acg, species,  prev_seq_region, next_seq_region, template_seq):

    resulststr = ""

    ###########################################################
    # determine the search region, adn extract it using fastascmd

    if not next_seq_region and not prev_seq_region:
        print "no regions specified in sw_searc()"
        return resulststr
    
    if not prev_seq_region:
        searchname   = next_seq_region.name
        searchfile   = next_seq_region.filename
        searchstrand = next_seq_region.strand
        if searchstrand==1:
            searchstart = next_seq_region.start - MAX_SEARCH_LENGTH
            searchend   = next_seq_region.start
        else:
            searchstart = next_seq_region.end
            searchend   = next_seq_region.end   + MAX_SEARCH_LENGTH

    elif not next_seq_region:
        searchname   = prev_seq_region.name
        searchfile   = prev_seq_region.filename
        searchstrand = prev_seq_region.strand

        if searchstrand==1:
            searchstart = prev_seq_region.end
            searchend   = prev_seq_region.end   + MAX_SEARCH_LENGTH
        else:
            searchstart = prev_seq_region.start - MAX_SEARCH_LENGTH
            searchend   = prev_seq_region.start

    else:
        if prev_seq_region.name == next_seq_region.name:
            searchname =  prev_seq_region.name
        else:
            print "prev_seq_region.name != next_seq_region.name ",  
            print prev_seq_region.name,  next_seq_region.name
            print "(cannot handle such cases yet)"
            return resulststr

        if prev_seq_region.filename == next_seq_region.filename:
            searchfile =  prev_seq_region.filename
        else:
            print "prev_seq_region.filename != next_seq_region.filename ",  
            print prev_seq_region.filename,  next_seq_region.filename
            print "(cannot handle such cases yet)"
            return resulststr

        if prev_seq_region.strand == next_seq_region.strand:
            searchstrand = prev_seq_region.strand
        else:
            print "prev_seq_region.strand != next_seq_region.strand ",  
            print prev_seq_region.strand,  next_seq_region.strand
            print "(cannot handle such cases yet)"
            return resulststr

        searchstart = prev_seq_region.end
        searchend   = next_seq_region.start

        if searchstart > searchend: # we are on the other strand
            searchstart = next_seq_region.end
            searchend   = prev_seq_region.start
        

    print searchname, searchfile, searchstrand, searchstart, searchend

    searchtmp = NamedTemporaryFile(delete=True)
    querytmp  = NamedTemporaryFile(delete=True)

    ###########################################################
    # extract search region  using fastacmd
    fastacmd = acg.generate_fastacmd_gene_command(species, searchname, searchfile, searchstrand, searchstart, searchend)
    #fasta = subprocess.check_output(fastacmd, shell=True)
    fasta  = commands.getoutput(fastacmd)
    searchtmp.write(fasta)
    searchtmp.flush()

    
    ###########################################################
    # Write out query sequences(? how many of them?)
    querytmp.write(">{0}\n{1}\n".format("query", template_seq)) # dna_seq
    querytmp.flush()

    ###########################################################
    # do  SW# search
    swsharpcmd = acg.generate_SW_nt(querytmp.name, searchtmp.name)
    resultstr  = commands.getoutput (swsharpcmd)
    searchtmp.close()
    querytmp.close()
    
    print resultstr
    exit(1)

    return resulststr
    
#########################################
def get_template (cursor, ensembl_db_name,  map_table, species, he):

    template_species = None
    template_seq     = None

    nearest_species = species_sort(cursor, map_table.keys(), species)[1:]
    exon = Exon()

    for nearest in nearest_species:
        if not map_table[nearest][he]: continue

        m = map_table[nearest][he]
        exon.exon_id     = m.exon_id_2
        exon.is_known    = m.exon_known_2
        exon.analysis_id = 1 # to be able to recycle get_exon_dna_seq()

        template_seq = get_exon_dnaseq (cursor, exon,  ensembl_db_name[species])
        if not template_seq:
            template_species = None
        else:
            template_species = nearest
            break


    return [template_species, template_seq]



#########################################
# he stands for "human exon"
def get_neighboring_region  (cursor, ensembl_db_name,  map_table, species, he, nbr_he):
    
    nbr_region = None

    if not nbr_he:
        # the bound will be my present left boud 
        #  + default range of search on the left
        pass
    else:
        nbr_map        = map_table[species][nbr_he]
        nbr_exon_id    = nbr_map.exon_id_2
        nbr_exon_known = nbr_map.exon_known_2

        nbr_region = get_seq_region(cursor, ensembl_db_name[species], nbr_exon_id, nbr_exon_known)

    return nbr_region


#########################################
def search_for_exons(human_gene_list, db_info):

    verbose = True

    # 
    [local_db, ensembl_db_name] = db_info
    if local_db:
        db  = connect_to_mysql()
        cfg = ConfigurationReader()
        acg = AlignmentCommandGenerator()
    else:
        db  = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        cfg = ConfigurationReader      (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    # find db ids and common names for each species db
    all_species, ensembl_db_name = get_species (cursor)

    # for each human gene
    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])

    gene_ct = 0
    found   = 0
    sought  = 0
    human_gene_list.reverse()
    #for human_gene_id in human_gene_list:
    #for human_gene_id in [370495]: # Known hit
    #for human_gene_id in [378768]: #  p53
    #for human_gene_id in [412667]: #  wls
    #for human_gene_id in [374433]: # Known hit
    #for human_gene_id in [397321]: # nice example, finds 23 out of 61

    for human_gene_id in [397176]:

	switch_to_db (cursor, ensembl_db_name['homo_sapiens'])

	# Get stable id and description of this gene -- DEBUG
	human_stable      = gene2stable    (cursor, human_gene_id)
        human_description = get_description(cursor, human_gene_id)
	if verbose:  print human_gene_id, human_stable, human_description

	# progress counter 
	gene_ct += 1
	if (not gene_ct%10): 
            print "processed ",   gene_ct, " out of ", len(human_gene_list), "genes"
            print "exons found: ", found, "out of ", sought, "sought"

	# find all human exons for this gene we are tracking in the database in order
	human_exons = [e for e in gene2exon_list(cursor, human_gene_id) if e.covering_exon < 0 and e.is_canonical and e.is_known]
	human_exons.sort(key=lambda exon: exon.start_in_gene)


	# make 'table' of maps, which is either pointer to the map if it exists, or None
	map_table  = {}
        for species in all_species:
            map_table[species] = {}
            for he in human_exons:
                map_table[species][he] = None

        for he in human_exons:
            maps_for_exon =  get_maps(cursor, ensembl_db_name, he.exon_id, he.is_known) # exon data
            for m in maps_for_exon:
               if m.similarity < 0.7: continue
               map_table[m.species_2][he] = m

        # get rid of species that do not have the gene at all
        for species in all_species:
            one_exon_found = False
            for he in human_exons:
                if map_table[species][he]:
                    one_exon_found = True
                    break
            if not one_exon_found:
                del map_table[species]
               
        # get rid of exons  that appear in no other species but human?
        bad_he = []
        for he in human_exons:
            one_species_found = False
            for species in  map_table.keys(): 
                if species =='homo_sapiens': continue
                if map_table[species][he]:
                    one_species_found = True
                    break
            if not one_species_found:
                bad_he.append(he)
        human_exons = filter (lambda he: not he in bad_he, human_exons)

        # keep track of nesrest neighbors for each human exon
        previous = {}
        next     = {}
        prev     = None
        for he in human_exons:
            previous[he]        = prev
            if prev: next[prev] = he
            prev = he
        next[he] = None

        # fill,  starting from the species that are nearest to the human
        species_sorted_from_human = species_sort(cursor,map_table.keys(),species)[1:]
        for species in species_sorted_from_human:
            # see which exons have which neighbors
            no_left  = []
            no_right = []
            has_both_neighbors = []
            for he in human_exons:
                m =  map_table[species][he]
                if m: continue
                prev = previous[he]
                nxt  = next[he]
                if prev and nxt and map_table[species][prev] and map_table[species][nxt]:
                    has_both_neighbors.append(he)
                elif not prev or not map_table[species][prev]:
                    no_left.append(he)
                elif not nxt or not map_table[species][nxt]:
                    no_right.append(he)
            
            # fill in exons that have both heigbors:
            if has_both_neighbors: print "\n", species, "\n both nbrs:"
            for he in has_both_neighbors:

                m = map_table[species][he]
                # get template (known exon from the nearest species)
                [template_species, template_seq] = get_template (cursor, ensembl_db_name,  map_table, species, he)
                if not template_species: continue

                print he.exon_id, "template species ", template_species
                # get previous region
                prev_seq_region = get_neighboring_region (cursor, ensembl_db_name, 
                                                          map_table, species, he, previous[he])
                next_seq_region = get_neighboring_region (cursor, ensembl_db_name, map_table, 
                                                          species, he, next[he])
                print "prev region: ", prev_seq_region
                print "next region: ", next_seq_region
                print 
                # do sw search
                resulststr = sw_search(acg, species, prev_seq_region, next_seq_region, template_seq)
                # parse the output
                
                # store if something found
                exit(1)
            print
            
            continue

            print "no left:",
            no_left.reverse()
            for he in no_left:
                m =  map_table[species][he]
                print he.exon_id,
            print
            
            print "no right:",
            for he in no_right:
                m =  map_table[species][he]
                print he.exon_id,
            print
            
                

        exit(1)
                

 
	#for species in all_species:
	for species in ['macropus_eugenii']:

            
            if species == "homo_sapiens" or species not in regions: continue # Skip if homo_sapiens or no exons for this species at all

            if verbose: print species

	    # Get sorted list of other species to search for exons for this species
	    switch_to_db (cursor, ensembl_db_name[species])
	    sorted_species_list = species_sort(cursor,all_species,species)[1:]

	    #for i in xrange(1,len(human_exons)-1): # we cannot find first or last exon, so skip them -- DEBUG
	    for he in  human_exons:

		

                sought += 1

                maps_to_human_exon_id = he.exon_id
                human_exon_pepseq     = get_exon_pepseq (cursor, he, ensembl_db_name['homo_sapiens'])

                if 0 and check_sw_exon_exists(cursor, ensembl_db_name, species, gene_id, maps_to_human_exon_id):
                    #print maps_to_human_exon_id, species, ">>> exon found"
                    found += 1
                    continue


		for gene_id  in regions[species].keys():


                    resultstr = sw_search(acg, exon_seq_id, protein_seq, species, searchname, searchfile, 
                                          searchstrand, searchstart, searchend)


		    for r in (f.splitlines() for f in resultstr.split("#"*80+"\n")):
			if len(r) < 14: continue # Skip blank or malformed results

			# Parse result
			seqlen = min(int(re.split('\D+',r[1])[1]),int(re.split('\D+',r[3])[1]))
			identity, matchlen = map(int,re.split('\D+',r[7])[1:3])
			#similarity = int(re.split('\D+',r[8])[1])
			#gaps = int(re.split('\D+',r[9])[1])
			#score = float(r[10].split()[1])

			# Reject if identity too low or too short
			if identity < 0.4*seqlen or identity < 10: continue

			# FOUND AN EXON!
			search_start, search_end, template_start, template_end = map(int,re.split('\D+',r[6])[1:5])

			# Lengthen find to align to coding frame in template
                        # what _is_ this?
			start_adjust    = (template_start-1) % 3
			search_start   -= start_adjust + 1
			template_start -= start_adjust + 1

			end_adjust      = -template_end % 3
			search_end     -= end_adjust
			template_end   -= end_adjust

			template_name      = int(r[2].split('>')[1])
			template_seq       = [s for s in nearestseqs if s[0] == template_name][0] 
			template_gene_id   = template_seq[7]
			template_searchseq = template_seq[6]
			template_dnaseq    = template_searchseq[template_start:template_end]

                        # why are we doing this?
			#if is_mitochondrial(cursor, template_gene_id, ensembl_db_name[nearestspecies]):
			#    template_pepseq    = str(Seq(template_dnaseq,generic_dna).translate(table="Vertebrate Mitochondrial"))
			#else:
			#    template_pepseq    = str(Seq(template_dnaseq,generic_dna).translate())
                        template_pepseq   = template_seq[1]


			searchseq = "".join(fasta.splitlines()[1:])
                        # how different is translation for the translated template?
                        # align in all three frames and pick one which is the most similar to the template 
                        mitochondrial = is_mitochondrial(cursor, gene_id, ensembl_db_name[species])
                        [pepseq, search_start, search_end] = find_the_most_similar_frame (cursor, mitochondrial,
                                                                                          template_pepseq, searchseq, 
                                                                                          search_start, search_end)
                        if len(pepseq)<3: continue
			flanklen    = 60
			left_flank  = searchseq[search_start-flanklen if search_start>flanklen else 0 : search_start]
			right_flank = searchseq[search_end : search_end+flanklen if search_end+flanklen<len(searchseq) else len(searchseq)]
			has_NNN  = 1 if "NNN" in searchseq else 0
			has_stop = 1 if "*" in pepseq else 0

			prev_end_in_gene = max((get_exon(cursor,m.exon_id_2,db_name=ensembl_db_name[species]).end_in_gene
						for m in maps[positions[p-1]] if m.species_2 == species))
			start_in_gene = prev_end_in_gene + search_start + 1 # Do I need this last +1 ?
			end_in_gene   = prev_end_in_gene + search_end

                        dnaseq = searchseq [search_start:search_end]
                        found += 1
                        if not pepseq == translate(searchseq, search_start, search_end, mitochondrial):
                            print " ! "
                            print pepseq
                            print translate(searchseq, search_start, search_end, mitochondrial)
                            continue


			strand = -1 if searchstrand != 1 else 1
                        sw_exon_id  = None
			#sw_exon_id = store_sw_exon (cursor,ensembl_db_name[species], he.exon_id, gene_id, 
			#	       start_in_gene, end_in_gene, strand, dnaseq, left_flank, right_flank, pepseq,
                        #               has_NNN, has_stop, template_name, nearestspecies)

			if verbose:# Print out some debugging info
                            print "============================================"
                            print human_gene_id, human_stable, human_description 
                            print "found sequence for {0}:".format(species)
                            print dnaseq
                            print "left flank:    " + left_flank 
                            print "right flank:   " + right_flank
                            print "translation:   " + pepseq + (" (Stop codon)" if has_stop == 1 else "")
                            print "template:      " + template_pepseq
                            print "human version: " + human_exon_pepseq
                            print "based on {0}, exon {1}, position {2} - {3}:".format(nearestspecies,template_name,
                                                                                       template_start,template_end)
                            print "storing to ", ensembl_db_name[species]
                            print "stored as exon ", sw_exon_id
                            print "exons found: ", found, "out of ", sought, "sought"
                            print 


        exit(1)

#########################################
def main():
    
    no_threads = 1

    local_db = False

    if local_db:
        db = connect_to_mysql()
    else:
        db = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    [all_species, ensembl_db_name] = get_species (cursor)


    species       = 'homo_sapiens'
    switch_to_db (cursor, ensembl_db_name[species])
    gene_list     = get_gene_ids (cursor, biotype='protein_coding', is_known=1)
    cursor.close()
    db.close()

    parallelize (no_threads, search_for_exons, gene_list, [local_db, ensembl_db_name])
    
    return True

#########################################
if __name__ == '__main__':
	main()

'''
       for species in species_sorted_from_human:
            for he in human_exons:
                m =  map_table[species][he]
                if not m:
                    print "%6s " % 'none',
                else:
                    print "%6.2f " % m.similarity,
            print
'''
