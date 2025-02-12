from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord


def get_seq_record_by_id(alignment: MultipleSeqAlignment,target_sequence_id: str ) -> SeqRecord:
    target_sequences = [seq_record for seq_record in alignment if seq_record.id == target_sequence_id]

    if len(target_sequences) == 0:
        raise ValueError(f"Sequence with ID '{target_sequence_id}' not found in alignment.")
    if len(target_sequences) >  1:
        raise ValueError(f"Multiple sequences with ID '{target_sequence_id}' found in alignment.")
    return target_sequences[0]
