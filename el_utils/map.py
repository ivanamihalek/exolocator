
import os
from ensembl import genome_db_id2species
from mysql   import switch_to_db, search_db

#########################################################
class Map:    # this particular map is between exons

    #################################################
    # print
    def __str__ (self):

        printstr = ""

        for attr, value in self.__dict__.iteritems():
            if ( not value is None):
                printstr += attr+ "    "+ str(value)
            else:
                printstr += attr+ "    None"
                
            printstr += "\n"
        return printstr

    ###################################
    # when something is defined as an exon ....
    def __init__ (self):
        
        self.species_1          = 'homo_sapiens'
        self.species_2          = None
        self.exon_id_1          = None
        self.exon_id_2          = None
        self.exon_known_1       = None
        self.exon_known_2       = None
  
        self.cigar_line         = None
        self.similarity         = None
        self.bitmap             = None

    ###################################
    def load_from_db (self, db_row, cursor, paralogue=False):
        
        if (paralogue):
            [exon_map_id, exon_id, cognate_exon_id,  
             exon_known, cognate_exon_known,  
             cigar_line, similarity, source, msa_bitmap] = db_row
        else:
            [exon_map_id, exon_id, cognate_exon_id,  
             exon_known, cognate_exon_known, cognate_genome_db_id, 
             cigar_line, similarity, source, msa_bitmap] = db_row
            self.species_2      = genome_db_id2species (cursor,  cognate_genome_db_id)

        self.exon_id_1          = exon_id
        self.exon_id_2          = cognate_exon_id
        self.exon_known_1       = exon_known
        self.exon_known_2       = cognate_exon_known  
        self.cigar_line         = cigar_line
        self.similarity         = similarity
        self.bitmap             = msa_bitmap
       
#########################################
def get_maps(cursor, ensembl_db_name, exon_id, is_known, species = 'homo_sapiens', table='exon_map'):
    
    maps = []

    switch_to_db (cursor,  ensembl_db_name[species])
    qry  = "select * from %s where exon_id = %d " % (table, exon_id)
    qry += " and exon_known = %d " % is_known
    rows = search_db (cursor, qry)
    if not rows or "ERROR" in rows[0]:
        return []

    for row in rows:
        map = Map()
        paralogue = (table=='para_exon_map')
        map.load_from_db(row, cursor, paralogue)
        maps.append(map)

    return maps
