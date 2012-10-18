
#########################################################
class Exon:

    #################################################
    # the ensembl row is assumed to have been 
    # obtained by select * from exon
    def load_from_ensembl_exon(self, gene_start, ensembl_row):
 
        self.exon_id           = ensembl_row[0]
        self.start_in_gene     = ensembl_row[2]-gene_start
        self.end_in_gene       = ensembl_row[3]-gene_start
        self.strand            = ensembl_row[4]
        self.phase             = ensembl_row[5]
        self.end_phase         = ensembl_row[6]
        self.is_constitutive   = ensembl_row[8]
        #  known (that's the source indicator - we 
        # got it from table called 'exon')
        self.is_known          = 1
        
        
        return True

    #################################################
    # the ensembl row is assumed to have been 
    # obtained by select * from exon_pprediction
    def load_from_ensembl_prediction(self, gene_start, ensembl_row):

        self.exon_id           = ensembl_row[0]
        self.start_in_gene     = ensembl_row[4]-gene_start
        self.end_in_gene       = ensembl_row[5]-gene_start
        self.strand            = ensembl_row[6]
        self.phase             = ensembl_row[7]
        # not known (that's the source indicator - we 
        # got it from table called 'prediction_exon')
        self.is_known          = 0
        self.is_canonical      = 0
        self.is_constitutive   = 0
        
        return True


    #################################################
    # load in from gene2exon table (select * from gene2exon)
    def load_from_gene2exon (self, in_list):

        if ( len(in_list) < 13):
            print "error loading exon: the in list must be",
            print " at least 13 elements long"
            return False
        
        self.gene_id             = in_list[0]
        self.exon_id             = in_list[1]
        self.start_in_gene       = in_list[2]
        self.end_in_gene         = in_list[3]
        self.exon_seq_id         = in_list[4]
        self.strand              = in_list[5]
        self.phase               = in_list[6]
        self.is_known            = in_list[7]
        self.is_coding           = in_list[8]
        self.is_canonical        = in_list[9]
        self.is_constitutive     = in_list[10]
        self.covering_exon       = in_list[11]
        self.covering_exon_known = in_list[12]
        self.analysis_id         = in_list[13]

        return True

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
    # when somrthing is defined as an exon ....
    def __init__ (self):
        
        self.exon_id             = 0
        self.gene_id             = 0
        self.start_in_gene       = 0
        self.end_in_gene         = 0
        self.exon_seq_id         = 0
        self.strand              = 0
        self.phase               = 0
        self.end_phase           = 0
        self.translation_starts  = 0
        self.translation_ends    = 0
        self.is_known            = 0
        self.is_coding           = 0
        self.is_canonical        = 0
        self.covering_exon       = 0
        self.covering_exon_known = 0
        self.analysis_id         = 0

