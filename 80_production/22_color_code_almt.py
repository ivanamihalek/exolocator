#!/usr/bin/env python

import sys
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

from el_utils.mysql import mysql_using_env_creds, mysql_server_conn_close
from el_utils.tree import species_tree

# IBM color blind palette
# https://www.color-hex.com/color-palette/1044488
'''
yello   #ffb000 	(255,176,0)
orange	#fe6100 	(254,97,0)
magenta	#dc267f 	(220,38,127)
purplish #785ef0 	(120,94,240)
powder  #648fff 	(100,143,255)

'''
def color_code(alignment:MultipleSeqAlignment):
    pass

def main():
    """
    Main function to read alignment, restricts the alignment to the target sequence,
    remove gappy sequences, and  write the output.
    """
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_file> <output_file> ")
        sys.exit(1)

    [input_file, output_file] = sys.argv[1:3]

    try:
        # Read the alignment from the input file
        alignment = AlignIO.read(input_file, "fasta")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)

    species_list = [record.id for record in alignment]
    print(species_list)
    exit()
    cursor = mysql_using_env_creds()
    tree = species_tree(cursor, species_list)
    mysql_server_conn_close(cursor)

    # assign color code to each position in the alignment according to conservation
    color_code(alignment)

    print(f"Successfully processed {input_file} and wrote output to {output_file}")



if __name__ == "__main__":
    main()
