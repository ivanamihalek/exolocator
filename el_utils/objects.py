
#########################################################
class Exon:

    #################################################
    # the ensembl row is assumed to have been 
    # obtained by select * from exon
    def load_from_ensembl_prediction(self, gene_start, ensembl_row):
        
        self.start_in_gene     = ensembl_row[4]-gene_start
        self.end_in_gene       = ensembl_row[5]-gene_start
        self.strand            = ensembl_row[6]
        self.phase             = ensembl_row[7]

        return True


    #################################################
    # load in from gene2exon table (select * from gene2exon)
    def load_from_gene2exon (self, in_list):

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

    #################################################
    # print
    def __str__ (self):

        printstr = ""
        if (self.gene_id):
            printstr += "gene_id  "+       str(self.gene_id)
        else:
            printstr += "gene_id  "+       "None"
        printstr += "\n"

        if (self.exon_id):
            printstr += "exon_id  "+       str(self.exon_id)
        else:
            printstr += "exon_id  "+       "None"
        printstr += "\n"

        if (self.start_in_gene):
            printstr += "start_in_gene  "+  str(self.start_in_gene)
        else:
            printstr += "start_in_gene  "+ "None"
        printstr += "\n"


        if (self.end_in_gene):
            printstr += "end_in_gene  "+   str(self.end_in_gene)
        else:
            printstr += "end_in_gene  "+    "None"
        printstr += "\n"

        if (self.exon_seq_id):
            printstr += "exon_seq_id  "+   str(self.exon_seq_id)
        else:
            printstr += "exon_seq_id  "+    "None"
        printstr += "\n"

        if (self.strand):
            printstr += "strand  "+        str(self.strand)
        else:
            printstr += "strand  "+         "None"
        printstr += "\n"

        if (self.phase):
            printstr += "phase  "+         str(self.phase)
        else:
            printstr += "phase  "+          "None"
        printstr += "\n"

        if (self.is_known):
            printstr += "is_known  "+      str(self.is_known)
        else:
            printstr += "is_known  "+       "None"
        printstr += "\n"

        if (self.is_coding):
            printstr += "is_coding  "+     str(self.is_coding)
        else:
            printstr += "is_coding  "+      "None"
        printstr += "\n"

        if (self.is_canonical):
            printstr += "is_canonical  "+  str(self.is_canonical)
        else:
            printstr += "is_canonical  "+   "None"
        printstr += "\n"

        if (self.covering_exon):
            printstr += "covering_exon  "+ str(self.covering_exon)
        else:
            printstr += "covering_exon  "+  "None"
        printstr += "\n"

        if (self.covering_is_known):
            printstr += "covering_is_known  "+ str(self.covering_is_known)
        else:
            printstr += "covering_is_known  "+  "None"
        printstr += "\n"

        if (self.analysis_id):
            printstr += "analysis_id  "+   str(self.analysis_id)
        else:
            printstr += "analysis_id  "+    "None"
        printstr += "\n"


        return printstr

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

