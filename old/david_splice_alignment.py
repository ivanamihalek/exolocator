#!/usr/bin/python

from sys import argv, exit
import re
from el_utils.mysql   import connect_to_mysql, switch_to_db, search_db, store_or_update
from el_utils.ensembl import get_species, gene2stable, exon2stable, get_exon, get_exon_seqs

def main():
  if len(argv) < 2: exit(1) # No gene_id provided
  gene_id = argv[1]

  if   re.match("ENSG\d{11}$",    gene_id): species = "homo_sapiens"
  elif re.match("ENSMUSG\d{11}$", gene_id): species = "mus_musculus"
  else: exit(2) # gene_id is not a valid Ensembl id (die)

  exodb = connect_to_mysql(user='marioot', passwd='tooiram')
  exocur = exodb.cursor()
  switch_to_db(exocur, 'exolocator_db')

  # Get list of transcript_id for this gene_id from splicing table
  transcript_ids = search_db(exocur, "SELECT transcript_id FROM splicing WHERE gene_id = '%s'" % gene_id)
  if not transcript_ids: exit(3) # No transcripts for this gene_id (die)

  alignment = dict()
  for transcript_id in transcript_ids:
    # Get list of exon_id for this transcript_id from splicing_transcript table
    exon_ids = search_db(exocur, "SELECT exon_id, is_covered FROM splicing_transcript WHERE transcript_id = '%s'" % transcript_id)
    if exon_ids:
      alignment[transcript_id] = dict()
      for exon_id in exon_ids:
        # Select the correct table where information on this exon is stored
        if exon_id[1] = "1": dbname = "exon_covered"
        else:                dbname = "exon_" + species
        # Get start, end and dna_seq for this exon
        exon = search_db(exocur, "SELECT start_in_gene, end_in_gene, dna_seq FROM %s WHERE exon_id = '%s' LIMIT 1" % (dbname, exon_id[0]))
        if exon:
          alignment[transcript_id][exon_id[0]] = exon
        else: # No info for this exon (not a fatal error, do nothing)
    else: # No exons for this transcript (not a fatal error, do nothing)
        
    for transcript_id, exons in alignment.iteritems():
      for exon_id, (start, end, dna_seq) in exons.iteritems():
        # Assemble the alignment based on start, end and dna_seq 
        # ***FIXME***

    f.close()

if __name__ == '__main__':
  main()
