#!/usr/bin/env python

import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import List, Dict

from el_utils.bioinfo import get_seq_record_by_id


def restrict_to_target(alignment: MultipleSeqAlignment, target_sequence_id: str) -> MultipleSeqAlignment:
    """
    Removes columns from an alignment where the target sequence has a gap.

    Args:
        alignment (Bio.Align.MultipleSeqAlignment): An alignment object.
        target_sequence_id (str): The ID of the sequence to use as a reference for gap removal.

    Returns:
        Bio.Align.MultipleSeqAlignment: A new alignment with gapped columns removed.
    """

    # Find the target sequence
    target_sequence =  get_seq_record_by_id(alignment, target_sequence_id)

    # Identify columns to keep (where the target sequence has no gap)
    columns_to_keep = [i for i, char in enumerate(target_sequence) if char != '-' and char != '.']

    # Create a new alignment with only the selected columns
    # Build new records with only the selected columns
    new_records: List[SeqRecord] = []
    for record in alignment:
        new_seq_str: str = ''.join(str(record.seq)[i] for i in columns_to_keep)
        new_seq: Seq = Seq(new_seq_str)
        new_record: SeqRecord = SeqRecord(new_seq, id=record.id, name=record.name, description=record.description)
        new_records.append(new_record)

    return MultipleSeqAlignment(new_records)


def remove_gappy_sequences(alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
    """
   Removes sequences with more than 33% gaps ('X', '.', or '-') from a MultipleSeqAlignment object.

   Args:
       alignment: A Bio.Align.MultipleSeqAlignment object.

   Returns:
       A new Bio.Align.MultipleSeqAlignment object with gappy sequences removed.
   """
    filtered_sequences = []
    for seq_record in alignment:
        sequence = str(seq_record.seq)
        # Count gaps represented by 'X', '.', or '-'
        gap_count = sum([sequence.count(gap) for gap in ['X', '.',  '-']])
        if gap_count / len(sequence) > 0.33: continue
        filtered_sequences.append(seq_record)

    return MultipleSeqAlignment(filtered_sequences)


def prune_alignment(alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
    """
    Prunes a MultipleSeqAlignment by removing sequences that deviate significantly
    from the most common characters in each column.

    Args:
        alignment: The input MultipleSeqAlignment object.

    Returns:
        A new MultipleSeqAlignment object with pruned sequences.
    """

    def find_top_characters(alignment_column: str) -> List[str]:
        """
        Finds the most frequent character(s) in an alignment column.

        Args:
            alignment_column: A string representing a single column from the alignment.

        Returns:
            A list of characters that appear most frequently in the column.
            If there's a tie, all characters with the maximum count are returned.
        """
        char_counts: Dict[str, int] = {}
        for char in alignment_column:
            char_counts[char] = char_counts.get(char, 0) + 1
        max_count = max(char_counts.values())

        top_characters: List[str] = [c for c, count in char_counts.items() if count == max_count]
        return top_characters

    # 1. Find top characters for each column
    top_chars_per_column: List[List[str]] = []
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        top_chars_per_column.append(find_top_characters(column))

    # 2. Identify sequences to remove
    sequences_to_keep = []
    num_columns = alignment.get_alignment_length()

    for i in range(len(alignment)):
        sequence: SeqRecord = alignment[i]
        deviant_column_count = 0
        for j in range(num_columns):
            if sequence[j] not in top_chars_per_column[j]:
                deviant_column_count += 1

        if deviant_column_count / num_columns < 0.6:
            sequences_to_keep.append(i)
        else:
            print(f"dropping the oddball seq: {sequence.id}")

    # 3. Create a new alignment with the filtered sequences
    pruned_alignment = MultipleSeqAlignment([alignment[i] for i in sequences_to_keep])

    return pruned_alignment


def main():
    """
    Main function to read alignment, restricts the alignment to the target sequence,
    remove gappy sequences, and  write the output.
    """
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_file> <output_file> <target_sequence_id>")
        sys.exit(1)

    [input_file, output_file, target_sequence_id] = sys.argv[1:4]

    try:
        # Read the alignment from the input file
        alignment = AlignIO.read(input_file, "fasta")

        # Remove columns where the target sequence has gaps
        new_alignment =  remove_gappy_sequences(restrict_to_target(alignment, target_sequence_id))
        # remove oddballs
        new_alignment = prune_alignment(new_alignment)

        # Write the new alignment to the output file
        AlignIO.write(new_alignment, output_file, "fasta")

        print(f"Successfully processed {input_file} and wrote output to {output_file}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
