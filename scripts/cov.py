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
from identical_sequence_parser import IdenticalSequencesParser

from dataclasses import astuple, dataclass


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


def pad_seq(sequence: Seq) -> Seq:
    """ Pad sequence to multiple of 3 with N

        Args:
            Sequence [str]: Amino acid sequence

        Returns:
            sequence [str]: Padded sequence
    """
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))



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



def main():
    """ 1. It takes the GISAID DNA FASTA file and translates it to protein, given the open reading frames.
        2. The results are written to a file in FASTA format
    # """
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
    """Translate the GISAID genome into ORFs. 
    """
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



def align_sequences_globally():
    seq_a_file = "spike_uniprot.fasta"
    seq_b_file = "gisaid_results/translated.fasta"
    out_file = "gisaid_results/needle.fasta"
    emboss_needle(seq_a_file, seq_b_file, out_file)


def identity_calculation() -> List:
    path_to_needle = "/home/nadzhou/DEVELOPMENT/tmp/gisaid_results/needle.fasta"

    print("Intializing trimming of the aligned sequences...")
    result_record = trim_sequences(path_to_needle)
    return result_record


def trim_sequences(filepath: str) -> List:
    """ Parse the Needle file and look for highly identical matches

        Args: 
            filepath [str]: File path

        Returns: 
            result_record [list]: List of records of best matching sequences
    """
    needle_record = list(AlignIO.parse(filepath, "msf"))

    result_record = []

    for rec in needle_record[1:]: 
        reference_seq = rec[0]
        seq_parser = IdenticalSequencesParser(reference_seq, rec[1])

        result = seq_parser.highly_identical_seqs()

        if result:
            result_record.append(result)
    return result_record


def write_records_to_file(result_record: List, filename: str):
    """ Write records to file 

        Args: 
            result_record [list]: List of SeqRecords
            filename [str]: File path
    """ 
    with open(filename, "w") as file:
        for rec in result_record:
            SeqIO.write(rec, file, "fasta")
    print("Trimmed sequences written to file.")


if __name__ == '__main__':
    main()
