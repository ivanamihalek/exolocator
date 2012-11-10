
import os

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
  
        self.cigar_line         = None
        self.similarity         = None
