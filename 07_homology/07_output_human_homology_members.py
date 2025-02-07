#!/usr/bin/python3 -u
# we are collecting orthologues for human only
# the possibility that "A is ortholog of B, and B is ortholog of C, but A is not ortholog of C"
# is something we don't want to know about right now

from el_utils.ensembl import *


def main():
	# in version 101, this takes 4 CPU hrs, producing a file
	# that has 207M and 10,978,400 lines/entries;
	# the original homology_member file has 945 million entries
	# how does the Ensembl borwser use that - what am I missing?

	cursor = mysql_using_env_creds()

	ensembl_compara_name = get_compara_name(cursor)

	table_name = "homology_member_human"
	outf = open(f"raw_tables/{table_name}.tsv", "w")

	qry  = f"select gene_member_id, stable_id from {ensembl_compara_name}.gene_member "
	qry += "where taxon_id=9606 and biotype_group='coding'"
	for gene_member_id, stable_id in  hard_landing_search(cursor, qry):
		qry = f"select homology_id from {ensembl_compara_name}.homology_member where gene_member_id={gene_member_id}"
		ret =  error_intolerant_search(cursor, qry)
		if not ret: continue
		for line in ret:
			homology_id = line[0]
			qry = f"select homology_id, gene_member_id from {ensembl_compara_name}.homology_member where homology_id={homology_id}"
			ret = error_intolerant_search(cursor, qry)
			if not ret: continue
			print("\n".join(["\t".join([str(field) for field in  line]) for line in ret]), file=outf)

	mysql_server_conn_close(cursor)
	return True


#########################################
if __name__ == '__main__':
	main()

'''
this didn't work
create  table  scratch_gene_mem_subtable  engine=MYISAM  as select gene_member_id, 
stable_id from gene_member where taxon_id=9606  and biotype_group='coding' ;

Query OK, 23471 rows affected (0.66 sec)
Records: 23471  Duplicates: 0  Warnings: 0

'''

