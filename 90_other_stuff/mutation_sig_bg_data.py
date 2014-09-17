#!/usr/bin/python -u
# extract all info one needs to estimate the background mutation rate when analyzing tumors

import MySQLdb, commands, re

from   el_utils.mysql       import  *
from   el_utils.ensembl     import  *
from   el_utils.el_specific import  *
from   el_utils.utils   import  erropen
from   el_utils.map     import  Map, get_maps
from   el_utils.tree    import  species_sort
from   el_utils.ncbi    import  taxid2trivial
from   el_utils.almt_cmd_generator import AlignmentCommandGenerator
from   el_utils.config_reader      import ConfigurationReader

from Bio.Seq      import Seq
from bitstring import Bits




category_dict = {}
categories = []

def fill_category ():
    
    categories.extend(["AT_tsvn", "AT_tstn", "CG_tsvn", "CG_tstn", "CpG_tsvn", "CpG_tstn"])
    
    for from_nt in ['A', 'C', 'G', 'T']:
        category_dict[from_nt] = {}
        for to_nt in ['A', 'C', 'G', 'T']:
           category_dict[from_nt][to_nt] = {}
           category_dict[from_nt][to_nt][True]  = ""
           category_dict[from_nt][to_nt][False] = ""
    
    category_dict['A']['A'][False] = ""
    category_dict['A']['T'][False] = "AT_tsvn"
    category_dict['A']['C'][False] = "AT_tsvn"
    category_dict['A']['G'][False] = "AT_tstn"
    
    category_dict['T']['A'][False] = "AT_tsvn"
    category_dict['T']['T'][False] = ""
    category_dict['T']['C'][False] = "AT_tstn"
    category_dict['T']['G'][False] = "AT_tsvn"
    
    category_dict['C']['A'][False] = "CG_tsvn"
    category_dict['C']['T'][False] = "CG_tstn"
    category_dict['C']['C'][False] = ""
    category_dict['C']['G'][False] = "CG_tsvn"
    
    category_dict['C']['A'][True]  = "CpG_tsvn"
    category_dict['C']['T'][True]  = "CpG_tstn"
    category_dict['C']['C'][True]  = ""
    category_dict['C']['G'][True]  = "CpG_tsvn"
    
    category_dict['G']['A'][False] = "CG_tstn"
    category_dict['G']['T'][False] = "CG_tsvn"
    category_dict['G']['C'][False] = "CG_tsvn"
    category_dict['G']['G'][False] = ""
    
    category_dict['G']['A'][True] = "CpG_tstn"
    category_dict['G']['T'][True] = "CpG_tsvn"
    category_dict['G']['C'][True] = "CpG_tsvn"
    category_dict['G']['G'][True] = ""
    
    
#########################################
def mutate(codon, mutation_position, new_nt):
    mutated_codon = ""
    for k in range(3):
        if k==mutation_position:
            mutated_codon += new_nt
        else:
            mutated_codon += codon[k]
    return mutated_codon
   
    
#########################################
def main():

    verbose  = True
    local_db = False

    if local_db:
        db     = connect_to_mysql()
        acg    = AlignmentCommandGenerator()
    else:
        db     = connect_to_mysql         (user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
        acg    = AlignmentCommandGenerator(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()

    # find db ids adn common names for each species db
    [all_species, ensembl_db_name] = get_species (cursor)
    
    logf = erropen("error.log", "w") 
    if not logf: exit(1)
    
    outf = erropen("mut_significance_bg_data.txt", "w") 
    if not outf: exit(1)


    switch_to_db (cursor, ensembl_db_name['homo_sapiens'])
    #gene_ids = get_gene_ids (cursor, biotype='protein_coding', is_known=1, ref_only=True)
    gene_ids = [10093770, 10096272, 10096345, 10096675, 10099748, 10100334, 10102321,
                 10107648, 10109042, 10110817, 10110942, 10112233, 10113096, 10113985,
                 10114644, 10116865, 10116871, 10116960, 10117020, 10117806, 10118088,
                 10118560, 10119226, 10119621, 10119795, 10120755, 10120798, 10122731,
                 10124106, 10124118, 10124151, 10124517, 10124564, 10126667, 10127014,
                 10127114, 10127801, 10127963, 10127980, 10128472, 10128786, 10128793, 10131147, 10132254,
                 10132706, 10132994, 10133724, 10134243, 10135990, 10136396, 10136634, 10136887, 10137107,
                 10137840, 10137951, 10138257, 10140899, 10142006, 10142667, 10142688, 10142698, 10143238,
                 10143462, 10143600, 10143804, 10143858, 10143974, 10144456, 10145879, 10146358, 10147626,
                 10148092, 10149428, 10149492, 10149569, 10150429, 10150941, 10152510, 10152784, 10153063,
                 10153305, 10153773, 10155884]
    
    
    # the categories of mutations for which we will be collecting statistics
    fill_category ()    
    # for each human gene
    #gene_ids = [10093176 ]
    for gene_id in gene_ids:
       
        switch_to_db (cursor,  ensembl_db_name['homo_sapiens'])
        stable_id = gene2stable(cursor, gene_id)


        # find all canonical coding  human exons 
        # get_canonical_coding_exons also sorts exons by the start in the gene
        canonical_human_exons = get_canonical_coding_exons (cursor, gene_id, ensembl_db_name['homo_sapiens'])

        # bail out if there is a problem
        if not canonical_human_exons: continue

        full_reconstituted_cDNA = ""
        prev_codon_piece_plus_right_flank = ""
        for human_exon in canonical_human_exons:
            [exon_seq_id, pepseq, pepseq_transl_start, pepseq_transl_end, left_flank, right_flank, nucseq] = \
                    get_exon_seqs(cursor, human_exon.exon_id, human_exon.is_known)
            # add the split codon
            phase = get_exon_phase (cursor, human_exon.exon_id, human_exon.is_known)
            left_flank_plus_codon_piece = left_flank + nucseq[:pepseq_transl_start]
            split_codon = ""
            if phase > 0 and prev_codon_piece_plus_right_flank and left_flank:
                offset      = (3-phase)%3
                # hedge against the possibility that the translation starts
                # right at the start of the exon, but there is supposed to be a phase
                split_codon = prev_codon_piece_plus_right_flank[:phase] + left_flank_plus_codon_piece[-offset:]
            full_reconstituted_cDNA += split_codon + nucseq[pepseq_transl_start:pepseq_transl_end]
            prev_codon_piece_plus_right_flank = nucseq[pepseq_transl_end:] + right_flank
            
        mitochondrial = is_mitochondrial(cursor, gene_id);
        if (mitochondrial):
            full_reconstituted_seq = Seq(full_reconstituted_cDNA).translate(table="Vertebrate Mitochondrial").tostring()
        else:
            full_reconstituted_seq = Seq(full_reconstituted_cDNA).translate().tostring()
            
        canonical = get_canonical_transl (acg, cursor, gene_id, 'homo_sapiens', strip_X = False)
        if canonical[0] == 'X': #that's some crap apparently wrong transcript is annotated as canonical
            print >> logf,   "warning", gene_id, stable_id,  get_description (cursor, gene_id)
            print >> logf, "the deposited canonical sequence starts with X - is there an alternative (?)"
            canonical = canonical[1:]
            
        if full_reconstituted_seq[-1] == '*' and canonical[-1] != '*':
            canonical += '*'
        if ( len(full_reconstituted_seq) != len(canonical)  or  full_reconstituted_seq != canonical):
            
            print >> logf, "error" , gene_id, stable_id, get_description (cursor, gene_id)
            print >> logf, "error reassembling,  len(full_reconstituted_seq) != len(canonical) ",  len(full_reconstituted_seq) , len(canonical) 
            print >> logf, "canonical:"
            print >> logf, canonical
            print >> logf, "reconstituted:"
            print >> logf, full_reconstituted_seq
            continue

        # nucleotide stats
        count = {'A':0, 'C':0, 'C-CpG':0, 'T':0, 'G':0, 'G-CpG':0} 
        codons = map(''.join, zip(*[iter(full_reconstituted_cDNA)]*3))
        L = len(full_reconstituted_cDNA)
        is_CpG = {}
        for i in range(L):
            is_CpG[i] = False
            if full_reconstituted_cDNA[i] == 'A':
                count['A'] += 1
            elif full_reconstituted_cDNA[i] == 'T':
                count['T'] += 1
            elif full_reconstituted_cDNA[i] == 'C':
                if i + 1 < L and full_reconstituted_cDNA[i + 1] == 'G':
                    count['C-CpG'] += 1
                    is_CpG[i] = True
                else:
                    count['C'] += 1
            elif full_reconstituted_cDNA[i] == 'G':
                if i > 0 and full_reconstituted_cDNA[i - 1] == 'C':
                    count['G-CpG'] += 1
                    is_CpG[i] = True
                else:
                    count['G'] += 1
                    
        # in each category_dict (AT transt, AT transv, CG trans, CG transv, Cpg trans, cpGtransv, how many missense, 
        #  how many nonsense, how many silent  possible    
        silent   = {}
        missense = {}
        nonsense = {}
        for cg in categories:
            silent[cg] = 0
            missense[cg] = 0
            nonsense[cg] = 0
        for i in range(len(codons)):
            codon = codons[i]
            aa = full_reconstituted_seq[i]
            for j in range(3):
                nt_position = i*3 + j
                nt = full_reconstituted_cDNA[nt_position]
                for new_nt in ['A', 'C', 'T', 'G']:
                    if new_nt == nt: continue
                    mutated_codon = mutate(codon, j, new_nt)
                    if (mitochondrial):
                        mutated_aa = Seq(mutated_codon).translate(table="Vertebrate Mitochondrial").tostring()
                    else:
                        mutated_aa = Seq(mutated_codon).translate().tostring()
                    cg = category_dict[codon[j]][new_nt][is_CpG[nt_position]];
                    if not cg or not cg in categories:
                        print >> logf, "category problem in ", gene_id, stable_id, get_description (cursor, gene_id)
                        print >> logf, codon, mutated_codon, j, codon[j], new_nt, is_CpG[nt_position], cg
                        print >> logf, i, j, nt_position, nt
                        print >> logf, aa, mutated_aa
                        continue
                    if (mutated_aa == aa):
                        silent[cg] += 1
                    elif (mutated_aa == "*"):
                        nonsense[cg] += 1
                    else:
                        missense[cg] += 1
                
        print >> outf, stable_id, get_description (cursor, gene_id)
        print >> outf, "CpG nucleotides"
        for i in range(len(codons)):
            if (is_CpG[i]):
                print >> outf," %5d  %s  %s " % (i, full_reconstituted_cDNA[i], codons[i])
        print >> outf,"%10s  %5s  %5s  %5s" % ("category", "silent", "nonsense", "missense")
        for cg in categories:
            print >> outf,"%10s  %5d  %5d  %5d" % (cg, silent[cg], nonsense[cg], missense[cg])
        print >> outf, "done", stable_id
        #print  "canonical:"
        #print canonical
        #print "reconstituted:"
        #print full_reconstituted_seq
        

    logf.close()

#########################################
if __name__ == '__main__':
    main()

