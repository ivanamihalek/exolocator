
import os
from ncbi    import get_ncbi_tax_name, taxid2parentid, taxid2name
from ensembl import get_compara_name, species2taxid
from mysql   import switch_to_db

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
                if (not firstborn):
                    ret_string += ","
                ret_string += child_string
                firstborn = 0
            # right bracket
            ret_string += ")"    
            # my own name
            if ( self.name):
                ret_string += self.name

        return ret_string
        
    ###################################
    # this should really be operation on the
    # tree, but python woudn't let me define it there
    def  __cleanup__ (self):
        if  (self.is_leaf):
            return None
        for i in range(len(self.children)):
            child = self.children[i] 
            ret   = child.__cleanup__()
            if ret:
                self.children[i] = ret
        if (len(self.children) == 1):
            firstborn = self.children[0]
            self.children = []
            return firstborn


    ###################################
    # when something is defined as a Node ....
    def __init__ (self, name):

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
            print "lf:", leaf.name, leaf.tax_id, leaf.parent_id
    
    ###################################
    def nhx_string(self):
        return self.root.nhx_string()



    ###################################
    # construct tree using the ncbi tree and the given leafs
    def build (self, cursor):

        # get tax ids from the gene_db table in compara database
        switch_to_db(cursor, get_compara_name(cursor))
        for leaf in self.leafs:
            taxid = species2taxid (cursor, leaf.name)
            leaf.tax_id =  taxid
        # add leafs to the rest of the nodes
        for leaf in self.leafs:
            leaf.is_leaf = True
            self.node[leaf.tax_id] = leaf

        # build the tree using ncbi_tax.nodes
        # fill in the names using ncbi_tax.names
        switch_to_db(cursor, get_ncbi_tax_name(cursor))
        for leaf in self.leafs:
            parent_id = taxid2parentid (cursor, leaf.tax_id)
            leaf.parent_id = parent_id
            current_id     = leaf.tax_id
            # move to the root
            while current_id:
                current_node = self.node[current_id]
                parent_id    = taxid2parentid (cursor, parent_id)
                if (not parent_id or  current_id == parent_id):
                    current_node.is_root = True
                    self.root = self.node[current_id]
                    current_id = None

                else:

                    # does parent exist by any chance
                    if self.node.has_key(parent_id):
                        parent_node = self.node[parent_id]
                        parent_node.children.append(current_node)
                        # we are done here
                        current_id = None
                    else: # make new node
                        parent_name    = taxid2name(cursor, parent_id)
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
        for node_id, node in self.node.iteritems():
            if node.is_leaf:
                continue
            if (not node.children):
                del_ids.append(node_id)

        for node_id in del_ids:
            del self.node[node_id]
                              
        self.__set_parent_ids__ (self.root)

    ###################################
    def __set_parent_ids__ (self, node):
        if node.is_leaf:
            return
        for child in node.children:
            child.parent_id = node.tax_id
            child.parent    = node
            self.__set_parent_ids__ (child)
        return


    ###################################
    # when something is defined as a Tree ....
    def __init__ (self):
        
        self.root  = None
        self.node  = {}
        self.leafs = []

    


#########################################
def subtree_leafs (node):

    leafs = []
    if node.is_leaf:
        leafs.append (node.name)
    else:
        for child in node.children:
            leafs += subtree_leafs(child)
    return leafs
    
#########################################
def  find_cousins (qry_node):
    
    cousins = []
    if not qry_node.parent:
        return cousins
    else:
        for sibling in qry_node.parent.children:
            if sibling == qry_node: continue
            cousins += subtree_leafs(sibling)

    cousins += find_cousins (qry_node.parent)
        
    return cousins

#########################################
def  species_sort(cursor, all_species, qry_species):

    cousins = []
    tree   = Tree()
    for species in all_species:
        leaf = Node(species)
        tree.leafs.append(leaf)
        if leaf.name == qry_species:
            qry_leaf = leaf

    tree.build(cursor)
    #find cousins for the qry leaf (recursively)
    cousins = find_cousins (qry_leaf) 
    
    return [qry_species]+cousins

