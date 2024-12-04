from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def split_fasta(file_in, file_out, dict_positions):
    sequences = []

    record_dict = SeqIO.to_dict(SeqIO.parse(file_in, "fasta"))

    for id_seq, pos in dict_positions.items():
        print(id_seq)
        pos = [0] + pos
        seq_record = record_dict[id_seq]
        indices = sorted(pos)
        print(indices)
        parts = [seq_record.seq[i:j] for i,j in zip(indices, indices[1:]+[None])]
        for i, seq_part in enumerate(parts):
            seq_part_id = seq_record.id + "_misjoin{}".format(i)
            seq_part_record = SeqRecord(seq_part, id=seq_part_id, description="")
            sequences.append(seq_part_record)

    for id_seq, rec in record_dict.items():
        if id_seq not in dict_positions.keys():
            sequences.append(rec)

    with open(file_out, "w") as outgenome:
        SeqIO.write(sequences, outgenome, "fasta")


if __name__ == "__main__":
    file_in = sys.argv[1]
    file_out = sys.argv[2]
    dict_positions = {'NC_007429.1': [5892, 104006, 461731, 769196, 820083, 881521,938156,975938,1018520]}
    split_fasta(file_in, file_out, dict_positions)