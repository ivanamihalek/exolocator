from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# Create a sample MultipleSeqAlignment object
seq1 = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha")
seq2 = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta")
seq3 = SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma")

alignment = MultipleSeqAlignment([seq1, seq2, seq3])

# Define the list of columns you want to extract
columns_to_extract = [2, 5, 8]  # Example: extract columns 2, 5, and 8 (0-based indexing)

# Extract the specified columns
sub_alignment = alignment[:, [2, 3]]

# Print the sub-alignment
print(sub_alignment)
