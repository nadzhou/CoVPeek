from identical_sequence_parser import IdenticalSequencesParser
import pytest

from Bio import SeqIO
import numpy as np


@pytest.fixture
def make_seq_inst(): 
    seq_record = list(SeqIO.parse("test.fasta", "fasta"))
    ref_seq = seq_record[0]
    tar_seq = seq_record[1]

    id_inst = IdenticalSequencesParser(ref_seq, tar_seq)

    return id_inst

def np_seq(): 
    seqs = list(SeqIO.parse("test.fasta", "fasta"))

    return np.asarray((seqs[0].seq, seqs[1].seq), dtype="S1")


def test_seqs2np(make_seq_inst):
    new_seq = np_seq()
    test_np_seq = make_seq_inst.seqs2np()

    print(new_seq[:, 0])
    print(test_np_seq[:, 0])

    assert (test_np_seq[:, 0] == new_seq[:, 0]).all()



def test_trim_target_seq(make_seq_inst): 
    sample_seq = ['S', 'S']

    result = make_seq_inst._trim_target_seq(sample_seq)

    assert result == 'S'





