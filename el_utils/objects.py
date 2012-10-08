
#########################################################
class Exon:


    ###################################
    # load in from gene2exon table 
    def load (self, in_list):

        if ( len(in_list) < 13):
            print "error loading exon: the in list must be",
            print " at least 13 elements long"
            return False
        
        self.gene_id           = in_list[0]
        self.exon_id           = in_list[1]
        self.start_in_gene     = in_list[2]
        self.end_in_gene       = in_list[3]
        self.exon_seq_id       = in_list[4]
        self.strand            = in_list[5]
        self.phase             = in_list[6]
        self.is_known          = in_list[7]
        self.is_coding         = in_list[8]
        self.is_canonical      = in_list[9]
        self.covering_exon     = in_list[10]
        self.covering_is_known = in_list[11]
        self.analysis_id       = in_list[12]

        return True

    ###################################
    # print
    def __str__ (self):
        printstr = ""
        printstr += "gene_id  ",       self.gene_id
        printstr += "exon_id  ",       self.exon_id
        printstr += "start_in_gene  ", self.start_in_gene
        printstr += "end_in_gene  ",   self.end_in_gene
        printstr += "exon_seq_id  ",   self.exon_seq_id
        printstr += "strand  ",        self.strand
        printstr += "phase  ",         self.phase
        printstr += "is_known  ",      self.is_known
        printstr += "is_coding  ",     self.is_coding
        printstr += "is_canonical  ",  self.is_canonical
        printstr += "covering_exon  ", self.covering_exon
        printstr += "covering_is_known  ", self.covering_is_known
        printstr += "analysis_id  ",   self.analysis_id

        return True


    ###################################
    # when somrthing is defined as an exon ....
    def __init__ (self):
        
        self.gene_id           = None
        self.exon_id           = None
        self.start_in_gene     = None
        self.end_in_gene       = None
        self.exon_seq_id       = None
        self.strand            = None
        self.phase             = None
        self.is_known          = None
        self.is_coding         = None
        self.is_canonical      = None
        self.covering_exon     = None
        self.covering_is_known = None
        self.analysis_id       = None

        return True

