# in ensebl  101 says is_current is always set to 1 in ensembl dbs, but needed for otterlace dbs

from el_utils.mysql import get_column_names, error_intolerant_search


#########################################################
class Exon:

	#################################################
	# object with provided attribute names:
	def fields2attributes(self, field_names, fields):
		for i in range(len(field_names)):
			setattr(self, field_names[i], fields[i])

	# the ensembl row is assumed to have been
	# obtained by select * from exon
	def load_from_db(self, cursor, db_name, table, exon_id, column_names=None):
		if not column_names: column_names = get_column_names(cursor, db_name=db_name, table_name=table)
		qry = f"select * from {db_name}.{table} where exon_id={exon_id}"
		ret = error_intolerant_search(cursor, qry)
		if not ret:
			print(f"exon_id {exon_id} not found in   {db_name}.{table} ")
			return None
		elif len(ret)>1:
			print(f"multiple entries for exon_id {exon_id} found in {db_name}.{table} ")
			return None
		if len(column_names) != len(ret[0]):
			print(f"provided column names do not match {db_name}.{table} ")
			return None
		self.fields2attributes(column_names, ret[0])

	def set_within_gene_coords(self, gene_start):
		setattr(self, "start_in_gene", self.seq_region_start-gene_start)
		setattr(self, "end_in_gene", self.seq_region_end-gene_start)

	def set_gene_id(self, gene_id):
		setattr(self, "gene_id", gene_id)

	def start(self):
		return self.start_in_gene if self.canon_transl_start < 0 else self.canon_transl_start

	def end(self):
		return self.end_in_gene if self.canon_transl_end < 0 else self.canon_transl_end


	#################################################
	# print
	def __str__(self):
		printstr = ""
		for attr, value in self.__dict__.items():
			if value is None:
				printstr += " %-20s    None"  % attr
			else:
				printstr += " %-20s    %s" % (attr,  str(value))
			printstr += "\n"
		return printstr

	###################################
	# when something is defined as an exon ....
	def __init__(self):
		self.strand           = None
		self.gene_id          = None
		self.exon_id          = None
		self.pepseq           = None
		self.is_canonical     = False
		self.is_coding        = False
		self.seq_region_start = -1
		self.seq_region_end   = -1
		self.start_in_gene    = -1
		self.end_in_gene      = -1
		self.canon_transl_start  = -1
		self.canon_transl_end  = -1
