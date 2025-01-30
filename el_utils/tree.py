#!/usr/bin/python3

import os
from .ncbi    import get_ncbi_tax_name, taxid2parentid, taxid2sciname
from .ensembl import get_compara_name, species2taxid, get_species
from .mysql import switch_to_db, connect_to_mysql, error_intolerant_search
import re

#########################################################
class Node:
	###################################
	# nhx string (for output)
	def nhx_string (self):
		ret_string = ""
		# if I'm a leaf, spit out the names and return
		if (self.is_leaf):
			if (self.name):
				ret_string += self.name
		else:
			# otherwise check hat children are doing, recursively
			# left bracket
			ret_string += "("
			# all children separated by a comma
			firstborn = 1
			for child in self.children:
				child_string = child.nhx_string()
				if not firstborn:
					ret_string += ","
				ret_string += child_string
				firstborn = 0
			# right bracket
			ret_string += ")"
			# my own name
			if self.name:
				ret_string += self.name
		return ret_string

	###########################################################
	def from_nhx_string(self, nhxstr, name2node=None, leafs=None):
		nhxstr = nhxstr.replace(' ','') # I don trust anybody to have doe that for me
		if not nhxstr: return
		if not "(" in nhxstr and not ")" in nhxstr:
			# our input consist of a single leaf
			self.name = nhxstr
			return
		# otherwise, we must start with an open bracket:
		if nhxstr[0]!="(":
			print(f"improperly formatted nhx string: {nhxstr[:10]} open bracket expected")
			exit()
		node_payload_end = 0
		for i in range(1,len(nhxstr)):
			if nhxstr[-i] == ")":
				node_payload_end = -i
				break
			elif nhxstr[-i] == "(":
				print(f"improperly formatted nhx string: {nhxstr[i:]} closing bracket expected at the end")
				exit()

		self.name = nhxstr[node_payload_end+1:]
		if not self.name:
			print(nhxstr)
			exit()

		# now check what we have in the payload
		payload = nhxstr[1:node_payload_end]
		
		# at any point of descent I might have a mixed set of leafs and subtrees
		# find the  nests of inner brackets
		bracket_start = []
		bracket_end   = []
		open_brackets = 0
		for i in range(len(payload)):
			if payload[i]=="(":
				if open_brackets==0: bracket_start.append(i)
				open_brackets += 1
			elif payload[i]==")":
				if open_brackets==0:
					print(f"improperly formatted nhx string: {payload}")
					exit()
				open_brackets -= 1
				if open_brackets==0: bracket_end.append(i)


		# the number of bracket nests is the number of subtrees
		number_of_subtrees = len(bracket_start)
		# they should all be matching though
		if len(bracket_end)!=number_of_subtrees:
			print(f"improperly formatted nhx string: {payload}")
			exit()

		# there might be isolated leafs between the bracket nests/subtrees
		# they should be separated by commas
		comma_positions = []
		prev = 0
		for b in range(number_of_subtrees):
			nest_start = bracket_start[b]
			nest_end = bracket_end[b]
			for i in range(prev,nest_start):
				if payload[i] == ",": comma_positions.append(i)
			prev = nest_end

		for i in range(prev, len(payload)):
			if payload[i]==",": comma_positions.append(i)
		comma_positions.append(len(payload))


		prev_pos = 0
		for comma_position in comma_positions:
			# everything between the two commas is a leaf or a subtree
			chunk = payload[prev_pos:comma_position]
			if "(" in chunk: # this is subtree
				new_node = Node()
				new_node.from_nhx_string(chunk, name2node, leafs)
				if name2node is not None: name2node[new_node.name] = new_node
				self.children.append(new_node)

			else: # this is leaf name
				new_node = Node(chunk)
				new_node.is_leaf=True

				if leafs is not None:
					leafs.append(new_node)
				if name2node is not None: name2node[new_node.name] = new_node
				self.children.append(new_node)
			prev_pos = comma_position+1


		if not self.children:
			self.is_leaf = True
			if leafs: leafs.append(self)

		return

	def subtree_leafs(self):

		leafs = []
		if self.is_leaf:
			leafs.append(self.name)
		else:
			for child in self.children:
				leafs += child.subtree_leafs()
		return leafs

	###################################
	# this should really be operation on the
	# tree, but python wouldn't let me define it there
	def  __cleanup__ (self):
		if  (self.is_leaf):
			return None
		for i in range(len(self.children)):
			child = self.children[i]
			ret   = child.__cleanup__()
			if ret:
				self.children[i] = ret

		if len(self.children) == 1:
			firstborn = self.children[0]
			self.children = []
			return firstborn

	###################################
	# when something is defined as a Node ....
	def __init__ (self, name="anon"):

		self.name      = name

		self.tax_id    = None
		self.parent    = None
		self.parent_id = None
		self.is_leaf   = False
		self.is_root   = False
		self.children  = []
		self.payload   = {}


#########################################################
class Tree:

	def print_leafs(self):
		for leaf in self.leafs:
			print("leaf:", leaf.name, leaf.tax_id, leaf.parent_id)

	###################################
	def nhx_string(self):
		return self.root.nhx_string()


	def from_nhx_string(self, nhxstr):
		self.root = Node()
		self.root.is_root = True
		self.root.from_nhx_string(nhxstr.strip(), self.name2node, self.leafs)
		self.name2node[self.root.name] = self.root
		self.__set_parent_ids__(self.root)

	###################################
	def add(self, cursor, name):
		leaf = Node(name)
		# get tax ids from the gene_db table in compara database
		switch_to_db(cursor, get_compara_name(cursor))
		leaf.tax_id = species2taxid (cursor, leaf.name)
		leaf.is_leaf = True
		self.leafs.append(leaf)
		self.node[leaf.tax_id] = leaf

	###################################
	# construct tree using the ncbi tree and the given leafs
	def build(self, cursor):

		for leaf in self.leafs:
			leaf.tax_id = species2taxid(cursor, leaf.name)
			leaf.is_leaf = True
			self.node[leaf.tax_id] = leaf

		# build the tree using ncbi_tax.nodes
		# fill in the names using ncbi_tax.names
		for leaf in self.leafs:
			parent_id = taxid2parentid (cursor, leaf.tax_id)
			leaf.parent_id = parent_id
			current_id     = leaf.tax_id
			# move to the root
			while current_id:
				current_node = self.node[current_id]
				parent_id    = taxid2parentid (cursor, parent_id)
				if not parent_id or  current_id == parent_id:
					current_node.is_root = True
					self.root = self.node[current_id]
					current_id = None

				else:
					# does parent exist by any chance
					if parent_id in self.node:
						parent_node = self.node[parent_id]
						parent_node.children.append(current_node)
						# we are done here
						current_id = None
					else: # make new node
						parent_name    = taxid2sciname(cursor, parent_id)
						parent_node    = Node(parent_name)
						parent_node.tax_id = parent_id
						# grampa:
						parent_node.parent_id = taxid2parentid (cursor, parent_id)
						parent_node.children.append(current_node)
						self.node[parent_id]  = parent_node
						# attach the current node to the parent
						current_id = parent_id

		# shortcircuit nodes with a single child
		new_root = self.root.__cleanup__()
		if (new_root):
			new_root.is_root = True
			self.root = new_root

		del_ids  = []
		for node_id, node in self.node.items():
			if node.is_leaf:
				continue
			if (not node.children):
				del_ids.append(node_id)

		for node_id in del_ids:
			del self.node[node_id]

		self.__set_parent_ids__ (self.root)

	def get_node(self, name):
		return self.name2node.get(name, None)


	def invert_node_map(self):
		self.name2node = dict([(node.name, node) for node in self.node.values()])

	###################################
	def __set_parent_ids__ (self, node):
		if node.is_leaf:
			return
		for child in node.children:
			child.parent_id = node.tax_id
			child.parent    = node
			self.__set_parent_ids__(child)
		return

	###################################
	# when something is defined as a Tree ....
	def __init__ (self, nhxstr=None):
		self.root  = None
		self.node  = {}
		self.leafs = []
		self.name2node = {}
		if nhxstr:
			self.from_nhx_string(nhxstr)




#########################################
def find_cousins (qry_node):

	cousins = []
	if not qry_node.parent:
		return cousins
	else:
		for sibling in qry_node.parent.children:
			if sibling == qry_node: continue
			cousins += sibling.subtree_leafs()

	cousins += find_cousins (qry_node.parent)

	return cousins


def species_tree(cursor, all_species):
	tree = Tree()
	for species in all_species:
		leaf = Node(species)
		tree.leafs.append(leaf)
	tree.build(cursor)
	tree.invert_node_map()
	return tree


def species_sort(cursor, species_list, qry_species):
	tree = species_tree(cursor, species_list)
	qry_leaf = None
	for leaf in tree.leafs:
		if leaf.name == qry_species:
			qry_leaf = leaf
			break
	assert qry_leaf, f"in species_sort(), leaf not found for {qry_species}"
	#find cousins for the qry leaf (recursively)
	cousins = find_cousins(qry_leaf)

	return [qry_species]+cousins


#####################################        
if __name__ == '__main__':

	local_db = False
	db = connect_to_mysql()
	cursor = db.cursor()

	[all_species, ensembl_db_name] = get_species (cursor)

	tree   = Tree()
	for species in all_species:
		tree.add(cursor, species)
	tree.build(cursor)

	print(tree.nhx_string())
