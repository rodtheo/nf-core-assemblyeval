from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if __name__ == "__main__":
    file_in = sys.argv[1]
    file_splited = sys.argv[2]
    rec_original = SeqIO.to_dict(SeqIO.parse(file_in, "fasta"))
    rec_splited = SeqIO.to_dict(SeqIO.parse(file_splited, "fasta"))

    for rr_id, rr in rec_original.items():
        original_seq = str(rr.seq)
        splited_seq = ''
        for splited_id, r in sorted(rec_splited.items()):
            if splited_id.startswith(rr_id):
                splited_seq += str(r.seq)
        
        if original_seq == splited_seq:
            print('{} == {}'.format(rr_id, splited_id))
            print('Sequences are the same')