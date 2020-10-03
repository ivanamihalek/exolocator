#!/usr/bin/python3 -u
from el_utils.mysql import  *
from el_utils.ensembl import *
from config import Config
import re
from el_utils.tree import species_sort
from string import ascii_uppercase

# if I go to the depth of 4, as in
# path_pieces = [home] + list(symbol[:4]) + [symbol]
# the size of the empty tree structure is 119M
# depth 3:
# the size of the empty tree structure is 92M
def symbol_path(home, symbol):
	path_pieces = [home] + list(symbol[:3]) + [symbol]
	return "/".join(path_pieces)

def create_directory_tree(home, hgnc_symbols):
	# there are HGNC symbols that have only teo characters
	# print(min([len(s) for s in hgnc_symbols]))
	# we also seem to have some special characters in symbols
	# eg C4B_2, HLA-C, PCDHA
	# but linux appear to be ok with it
	# for symbol in hgnc_symbols:
	# 	if any(not c.isalnum() for c in symbol):
	# 		print(symbol)
	for symbol in hgnc_symbols:
		path = symbol_path(home, symbol)
		if not os.path.exists(path): os.makedirs(path)

#########################################
def main():

	home = "/storage/exolocator/alignments"
	if not os.path.exists(home):
		print(f"{home} not found")
		exit()

	db = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	qry = "select approved_symbol from identifier_maps.hgnc where locus_group='protein-coding gene'"
	hgnc_symbols = [r[0] for r in hard_landing_search(cursor,qry)]
	create_directory_tree(home, hgnc_symbols) # does nothing if all paths exist

	for symbol in hgnc_symbols:
		path = symbol_path(home, symbol)
		# orthologues
		cmd  = f"touch {path}/{symbol}_ortho_notes.txt"
		subprocess.call(["bash","-c", cmd])
		# zip orthologues
		# orthologue notes
		# move to path
		# paralogues
		# zip paralogues
		# move to path
	cursor.close()
	db.close()
	return True


#########################################
if __name__ == '__main__':
	main()
