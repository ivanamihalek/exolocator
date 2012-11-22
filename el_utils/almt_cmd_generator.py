#!/usr/bin/python
'''
Created on Apr 12, 2012

@author: Mario OT
'''

import os, re, commands
from   config_reader   import ConfigurationReader


###########
def isinteger(x):
    try:
        int(x)
    except:
        return False
    try:
        int(x) == x
    except:
        return False
    return True


###########
class AlignmentCommandGenerator(object):
    '''
    Generates commands for utilities that are used (blast, sw, genewise, fastacmd, formatdb)
    '''

    def __init__(self, user=None, passwd=None, host=None, port =None, check=True):
        '''
        Loads the utils configuration (utils.cfg)
        
        if (not os.path.isfile("../utils.cfg")):
            raise IOError("There is no utils.cfg file present in the project directory.")
        
        config = ConfigParser.RawConfigParser()
        config.read("../utils.cfg")
        '''
        
        self.configReader = ConfigurationReader(user, passwd, host, port, check)
        # blast tools
        
        
        self.fastacmd  = self.configReader.get_path('fastacmd')
        self.blastall  = self.configReader.get_path('mafft')

        self.blastp    = self.configReader.get_path('blastall')
        self.blastp    += " -p blastp -e "+ self.configReader.get_value('blastp_e_value')
        

        # ensembl database
        self.ensembldb      = self.configReader.get_path('ensembl_fasta')
        
        # Smith-Waterman
        self.sw_sharp       = self.configReader.get_path('sw#')
        # hacked blast, to be used by swsharp to align exons while respecting the boundaries
        self.blosum_matrix  = "{0}/{1}".format(self.configReader.get_path('resources'),
                                               self.configReader.get_value('blosum_hacked'))
        # mafft
        self.mafft          = self.configReader.get_path('mafft')
       
        
        
    def generate_fastacmd_gene_command (self, species, seq_name, fasta_db_file, 
                                   strand = None, 
                                   sequence_start = None, sequence_stop = None, 
                                   output_file_path = None):
        
        if (isinteger(seq_name)):
            seq_id_cmd = "-s 'lcl|%s' " % seq_name
        else:
            seq_id_cmd = "-s %s" % seq_name
        

        data_type_cmd = "-p F"
        database = "-d {0}/{1}/dna/{2}".format(self.ensembldb, species, fasta_db_file)
            
        if (strand == None or int(strand) == 1):
            strand_cmd = "-S 1"
        else:
            strand_cmd = "-S 2"
            
        if (sequence_start and sequence_stop):
            location_cmd = "-L %s,%s" % (sequence_start, sequence_stop)
        else:
            location_cmd = ""
            
        if output_file_path:
            output_cmd = "-o %s" % output_file_path
        else:
            output_cmd = ""


        fastacmd_cmd_line  = self.fastacmd
        fastacmd_cmd_line +=  " {0} {1} {2} {3} {4} {5}".format(database, seq_id_cmd, 
                                                               data_type_cmd, strand_cmd,
                                                               location_cmd, output_cmd)
        return fastacmd_cmd_line
        
    
    def generate_fastacmd_plain (self, database, seq_id, output_file_path = ''):
        
        database_cmd = "-d %s" % database
        if (isinteger(seq_id)):
            seq_id_cmd = "-s 'lcl|%s' " % seq_id
        else:
            seq_id_cmd = "-s %s" % seq_id
        if (output_file_path):
            output_cmd  = "-o %s" % output_file_path
        else:
            output_cmd  = ""

        fastacmd_cmd_line  = self.fastacmd
        fastacmd_cmd_line +=  " {0} {1} {2} ".format(database_cmd, seq_id_cmd, output_cmd)
        return fastacmd_cmd_line
 
    def generate_fastacmd_protein_command (self, protein_id, species_name, protein_type, output_file_path):
        
        data_type_cmd = "-p T"
        prot_id_cmd = "-s %s" % protein_id
        database = "-d %s" % self._generate_proteindb_file_name(species_name, protein_type)
        if (output_file_path):
            output_cmd = "-o %s" % output_file_path
        else:
            output_cmd = ""
        return "fastacmd {0} {1} {2} {3}".format(prot_id_cmd, data_type_cmd, database, output_cmd)
    
    
    def generate_blastp_plain (self, database, input_file, output_file):
        
        # sombbody hardcoded the  output format in the self.blastp
        blastp = self.blastp + " -m 8 "
        if output_file:
            cmd = "{0} -d {1} -i {2} -o {3}".format(blastp, database, input_file, output_file)
        else:
            cmd = "{0} -d {1} -i {2} ".format(blastp, database, input_file)
            
        return cmd

    def generate_blastp_command (self, database, input_file, output_file):
        if output_file:
            cmd = "{0} -d {1} -i {2} -o {3}".format(self.blastp, database, input_file, output_file)
        else:
            cmd = "{0} -d {1} -i {2} ".format(self.blastp, database, input_file)
            
        return cmd

    
    def generate_SW_nt (self, query_sequence_file, target_fasta_db_file, 
                             output_file, supress_stdout = True):
        # Matija's current implementation is switching the order
        cmd = "{0} -j {1} -i {2} --out {3}".format(self.sw_sharp, query_sequence_file, 
                                                   target_fasta_db_file, output_file)
        if supress_stdout:
            cmd += " > /dev/null"
        return cmd
    
    def generate_SW_peptide (self, query_sequence_file, target_fasta_db_file, 
                             output_file = None):
        cmd  = "{0} --verbose 0 --matrix-file {1}  ".format(self.sw_sharp, self.blosum_matrix)
        cmd += " -j {0} -i {1} ".format(query_sequence_file, target_fasta_db_file)
        cmd += " --out-type 1 --gap-open 3.0   "
        if output_file:
            cmd += " --out {0} ".format(output_file)
        return cmd
    
    
    def generate_formatdb_command (self, input_db_file, sequence_type):
        
        if sequence_type == "protein" or sequence_type == "P":
            cmd = "formatdb -i {0} -p T".format(input_db_file)
        else:
            cmd = "formatdb -i {0} -p F".format(input_db_file)
            
        return cmd
    
    def generate_mafft_command (self, input_file, output_file=None):
        '''
        If there are unusual characters (e.g., U as selenocysteine in protein sequence), 
        use the --anysymbol option. It accepts any printable characters (U, O, #, $, %, etc.; 
        0x21-0x7e in the ASCII code), execpt for > (0x3e).  
        They are scored equivalently to X.  Gap is - (0x2d), as in the default mode. 
        '''
        if output_file:
            return "{0} --quiet --anysymbol {1} > {2}".format(self.mafft, input_file, output_file)
        else:
            return "{0} --quiet --anysymbol {1} ".format(self.mafft, input_file)
    
    def _generate_genedb_file_name (self, species, sequence_type, sequence_id, masked):
        '''
        @param species: species name
        @param sequence_type: scaffold / chromosome...
        @param sequence_id: ensembl sequence ID
        @param masked: 0 if dna should not be masked, 1 if it should
        '''
        file_name = "{0}/{1}/dna".format(self.ensembldb, species.lower())
        # get the template name (dependent on the assembly)
        tmp_file=""
        for f in os.listdir(file_name):
            # check if we've stumbled on the protein file
            m = re.findall(".dna", f)
            if not m:
                continue
            tmp_file = f
        m = re.findall ('(.*).dna', tmp_file)   
        if (masked != 0):
            file_name = "%s/%s.dna_rm." % (file_name, m[0])
        else :
            file_name = "%s/%s.dna." % (file_name, m[0])
        if (sequence_type == 'chromosome'):
            file_name = "%schromosome.%s.fa" % (file_name, sequence_id)
        else :
            file_name = "%stoplevel.fa" % (file_name)

        if (not os.path.exists(file_name)):
            # find *dna.toplevel.fa file
            path =  "{0}/{1}/dna".format(self.ensembldb, species.lower())
            cmd = "ls "+path+"/*dna.toplevel.fa"
            retval =  commands.getoutput(cmd)
            if (retval):
                file_name = retval.rstrip()
                # otherwise just let the thing fail, I have no better idea
        return file_name
       
       
    def _generate_proteindb_file_name (self, species, protein_type):
        '''
        @param species: species name (ensembl)
        @param protein_type: all / abinitio
        @return: protein database name
        '''
        file_name = "%s/%s/pep" % (self.ensembldb, species.lower())
        tmp_file=""
        for f in os.listdir(file_name):
            # check if we've stumbled on the protein file
            m = re.findall(".pep", f)
            if not m:
                continue
            tmp_file = f
        m = re.findall ('(.*).pep', tmp_file)
        if (protein_type == "all"):
            file_name = "%s/%s.pep.all.fa" % (file_name, m[0])
        else:
            file_name = "%s/%s.pep.abinitio.fa" % (file_name, m[0])
            
        return file_name
    
def main():
    acg = AlignmentCommandGenerator()
    cmd = acg.generate_SW_command("query.fa", "target.fa", "output", True)
    print cmd
    print
    for [key, val] in acg.__dict__.iteritems():
        print key, val
    
    
if __name__ == '__main__':
    main()
