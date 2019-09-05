#!/usr/bin/python3
# http://phylotree.hyphy.org/

from el_utils.mysql   import connect_to_mysql
from el_utils.ensembl import get_species
from el_utils.tree    import Tree, Node
from config import Config


#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species (cursor)

	tree   = Tree()
	for species in all_species:
		leaf = Node(species)
		tree.leafs.append(leaf)

	tree.build(cursor)

	print()
	print(tree.nhx_string())
	print()

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
