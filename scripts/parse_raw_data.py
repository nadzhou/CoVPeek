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
from parser import parse_arguments


def remove_non_chars_from_seq(seq: str) -> str: 
    return Seq(re.sub(r'(\W)', '', str(seq)))

def find_patient_country(genome_id): 
    return str((re.findall(f'19/(\w*)/', genome_id))).strip("[]'")

def make_seqrecord(protein, patient_country, i, num): 
    return SeqRecord(Seq(protein), 
                    id=f"{num}-{i}-{patient_country}",
                    name=f"{num}-{i}-{patient_country}",
                    description=f"{num}-{i}-{patient_country}")


def extract_dna(genome_records: List[SeqRecord], out_file_path: Path) -> List[List[SeqRecord]]:
    """ Extract GISAID genome DNA and then translate each record
        into ORFs, and output the list.

        Args:
            genome_records: GISAID genomes in a SeqRecord

            out_file_path: Output file path

        Returns:
            records: Translated ORFs in a list
    """
    records = []

    print("Initiating translation...")

    for num, genome_record in enumerate(genome_records, start=1): 
        genome_record.seq = remove_non_chars_from_seq(genome_record.seq)
        saved_records = []

        orf_proteins = find_orfs_with_trans(genome_record.seq)
        patient_country = find_patient_country(genome_record.id)

        print(patient_country)

        for i, protein in enumerate(orf_proteins, start=1): 
            translated_pt_record = make_seqrecord(protein, 
                                                patient_country, 
                                                i, 
                                                num)

            saved_records.append(translated_pt_record)

        if len(saved_records) > 2:
            records.append(saved_records)

    print("Genome record translated.")
    return records


def pad_seq(seq: Seq) -> Seq:
    """ Pad seq to multiple of 3 with N

        Args:
            seq [str]: Amino acid seq

        Returns:
            seq [str]: Padded seq
    """
    remainder = len(seq) % 3
    return seq if remainder == 0 else seq + 'N' * (3 - remainder)


def find_orfs_with_trans(seq: Seq, trans_table: int = 1, min_protein_length: int = 100) -> List[str]:
    """ Code from the Biopython library for finding open reading frames

        Args:
            seq: Protein seq
            trans_table: Look-up table
            min_protein_length: Minimum protein length

        Returns:
            answer: List of protein seqs with different reading frames
    """
    answer = []
    seq = pad_seq(seq)
    seq_len = len(seq)

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
    """

    parser = parse_arguments()
    output_path = Path(parser.output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    
    emboss_out_file = output_path/"needle.fasta"
    translate_genome(parser.cov_genome_path, output_path)

    
    emboss_needle("../operations/gisaid_results/translated.fasta", 
                    f"{str(parser.uniprot_refseq_path)}", 
                    emboss_out_file)

    # Now for the identitiy calculation bit.
    results_record = identity_calculation(emboss_out_file)
    results_filename = output_path/"trimmed_seqs.fasta"

    if results_record:
        write_records_to_file(results_record, filename=results_filename)
        # Run a MAFFT alignment to pad the seqs into equal length
        mafft(output_path/"trimmed_seqs.fasta", output_path/"aligned.fasta")
    # # After this go to variation_parser.
    # else:   
    #     print("No hits found")


def translate_genome(genome_path, output_path):
    """Translate the GISAID genome into ORFs. 
    """
    # unused variable, probably needs to be deleted
    genome_record = list(SeqIO.parse(genome_path, "fasta"))
    # path is relative to scripts folder here, needs to be different in notebooks

    results = extract_dna(genome_record, gisaid_results_path)
    with open(output_path/"translated.fasta", "w") as file:
        for record in results:
            SeqIO.write(record, file, "fasta")


def align_seqs_globally():
    seq_a_file = "spike_uniprot.fasta"
    out_file = "gisaid_results/needle.fasta"

    seq_b_file = "gisaid_results/translated.fasta"
    emboss_needle(seq_a_file, seq_b_file, out_file)


def identity_calculation(path_to_needle: str) -> List:
    """ Calculate the identites for the seqs

        Returns: 
            result_record [list]
    """
     
    print("Intializing trimming of the aligned seqs...")
    result_record = trim_seqs(path_to_needle)
    return result_record


def trim_seqs(filepath: str) -> List:
    """ Parse the Needle file and look for highly identical matches

        Args: 
            filepath [str]: File path

        Returns: 
            result_record [list]: List of records of best matching seqs
    """
    needle_record = list(AlignIO.parse(filepath, "msf"))

    result_record = []

    for rec in needle_record[1:]: 
        reference_seq = rec[1]

        seq_parser = IdenticalSequencesParser(reference_seq, rec[0])

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

    print("Trimmed seqs written to file.")


if __name__ == '__main__':
    main()
