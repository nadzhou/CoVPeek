#! usr/bin/env python3
from typing import List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

import re
import numpy as np
from pathlib import Path

from Bio.SeqRecord import SeqRecord

from emboss import emboss_needle
from mafft import mafft

from dataclasses import astuple, dataclass


def pad_seq(sequence: Seq) -> Seq:
    """ Pad sequence to multiple of 3 with N

        Args:
            Sequence [str]: Amino acid sequence

        Returns:
            sequence [str]: Padded sequence
    """
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))


def main():
    """ 1. It takes the GISAID DNA FASTA file and translates it to protein, given the open reading frames.
        2. The results are written to a file in FASTA format
    """
    translate_genome()

    align_sequences_globally()

    # Now for the identitiy calculation bit.
    results_record = identity_calculation()
    results_filename = "gisaid_results/trimmed_seqs.fasta"
    if results_record:
        write_records_to_file(results_record, filename=results_filename)
        # Run a MAFFT alignment to pad the sequences into equal length
        mafft("gisaid_results/trimmed_seqs.fasta")
    # After this go to variation_parser.py


def translate_genome():
    genome_path = "/home/nadzhou/SEQs/CoV/raw_seqs/gisaid19.fasta"
    # unused variable, probably needs to be deleted
    uniprot_record = SeqIO.read("spike_uniprot.fasta", "fasta")
    genome_record = list(SeqIO.parse(genome_path, "fasta"))
    # path is relative to scripts folder here, needs to be different in notebooks
    gisaid_results_path = Path("../notebooks/gisaid_results/")
    results = extract_dna(genome_record, gisaid_results_path)
    with open(gisaid_results_path / "translated.fasta", "w") as file:
        for record in results:
            SeqIO.write(record, file, "fasta")


def extract_dna(genome_records: List[SeqRecord], out_file_path: Path) -> List[List[SeqRecord]]:
    """ Extract GISAID genome DNA and then translate each record
        into ORFs, and output the list.

        Args:
            genome_records: GISAID genomes in a SeqRecord

            out_file_path: Output file path


        Returns:
            records: Translated ORFs in a list
    """
    out_file_path.mkdir(parents=True, exist_ok=True)
    records = []

    print("Initiating translation...")

    for genome_record, num in zip(genome_records,
                                  range(len(genome_records) + 1)):

        genome_record.seq = Seq(re.sub(r'(\W)', '', str(genome_record.seq)))

        saved_records = []

        orf_proteins = find_orfs_with_trans(genome_record.seq)

        for protein, i in zip(orf_proteins, range(len(orf_proteins) + 1)):
            translated_pt_record = SeqRecord(Seq(protein), id=f"{num}gen_{i}_orf_{genome_record.id}",
                                             description=f"{num}gen_{i}{genome_record.description}",
                                             name=f"{num}gen_{i}{genome_record.name}")

            saved_records.append(translated_pt_record)

        if len(saved_records) > 2:
            records.append(saved_records)

    print("Genome record translated.")
    return records


def find_orfs_with_trans(seq: Seq, trans_table: int = 1, min_protein_length: int = 100) -> List[str]:
    """ Code from the Biopython library for finding open reading frames

        Args:
            seq: Protein sequence
            trans_table: Look-up table
            min_protein_length: Minimum protein length

        Returns:
            answer: List of protein sequences with different reading frames
    """
    answer = []
    seq = pad_seq(seq)
    seq_len = len(seq)
    print(seq_len)

    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append(trans[aa_start:aa_end])
                aa_start = aa_end + 1
    return answer


def align_sequences_globally():
    seq_a_file = "spike_uniprot.fasta"
    seq_b_file = "gisaid_results/translated.fasta"
    out_file = "gisaid_results/needle.fasta"
    emboss_needle(seq_a_file, seq_b_file, out_file)


def identity_calculation():
    path_to_needle = "/home/nadzhou/DEVELOPMENT/tmp/gisaid_results/needle.fasta"

    print("Intializing trimming of the aligned sequences...")
    result_record = trim_sequences(path_to_needle)
    return result_record


def trim_sequences(filepath: str) -> List:
    record = list(AlignIO.parse(filepath, "msf"))
    result_record = []
    for rec, num in zip(record[1:], range(len(record[1:]))):
        canonical_ = rec[0]
        ref_seq_len = len(canonical_.seq)

        result = percent_id_calc(canonical_, rec[1], f"{rec[0].id[:7]}_{num}", ref_seq_len)

        if result:
            result_record.append(result)
    return result_record



class IdenticalSeqSeparator: 

    def __init__(self, ref_seq, tar_seq): 
        self.ref_seq = ref_seq
        self.tar_seq = tar_seq
        self.ref_seq_len = len(self.ref_seq)

                
    def percent_id_calc(self):
        """ Calculate percent identity of given two sequences
        """

        np_seq = np.asarray((self.ref_seq.seq, 
                            self.tar_seq.seq))

        identical_aa_freq, trimd_tar_seq = np.apply_along_axis(_identity_calc, 0, np_seq)

        self.identical_aa_freq = identical_aa_freq.astype(np.int64)
        self.trimd_tar_seq = trimd_tar_seq.astype(np.str)


    def count_checker(self):     
        if not np.all(self.identical_aa_freq == 0):
            self.trimd_tar_seq = "".join(item for item in trimd_tar_seq)
            self.identity_score = np.true_divide(sum(identical_aa_freq), 
                                                        self.ref_seq_len)

            return self.write_trimmed_seqs()

    @staticmethod
    def _identity_calc(array):
        """ Calculate idenity if both sequences have a match identical_aa_freq it up.

            Args:
                array [nd array]: Pair of protein sequence characters
        """
        identical_aa_freq = 0
        seq2 = ""

        if array[0] == array[1] and array[0] != "-":
            seq2 = array[1]
            identical_aa_freq += 1

        return identical_aa_freq, seq2


    def write_trimmed_seqs(self):
        """ Write the sequences that have identity scores greater than 80 percent.
        """
        try: 
            if self.identity_score > 0.7 and len(self.trimd_tar_seq) > 50:
                print(f"hit hit {identity_score}")

                seq_record2 = SeqRecord(Seq(self.trimd_tar_seq),
                                        id=self.tar_seq.id,
                                        name=self.tar_seq.name,
                                        description=self.tar_seq.description)

                return seq_record2


def write_records_to_file(result_record, filename: str):
    with open(filename, "w") as file:
        for rec in result_record:
            SeqIO.write(rec, file, "fasta")
    print("Trimmed sequences written to file.")


if __name__ == '__main__':
    main()
