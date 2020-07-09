#! usr/bin/env python3
import logging
from pathlib import Path
from typing import List, Union

from Bio import SeqIO, AlignIO

from cli_parser import parse_arguments
from emboss import emboss_needle
from identical_sequence_parser import IdenticalSequencesParser
from mafft import mafft
from translate_genome import translate_genome


def main():
    """ 1. It takes the GISAID DNA FASTA file and translates it to protein, given the open reading frames.
        2. The results are written to a file in FASTA format
    """

    logging.basicConfig(level=logging.INFO)

    args = parse_arguments()
    output_path = Path(args.output_directory)
    output_path.mkdir(parents=True, exist_ok=True)

    translated_file = output_path / "translated.fasta"
    # translate_genome(args.cov_genome_path, translated_file)
    emboss_out_file = output_path / "needle.fasta"
    # pairwise alignment
    # emboss_needle(translated_file,
    #               args.uniprot_refseq_path,
    #               emboss_out_file)

    # Now for the identity calculation bit.
    results_record = identity_calculation(emboss_out_file)
    results_filename = output_path / "trimmed_seqs.fasta"

    if results_record:
        write_records_to_file(results_record, filename=results_filename)
        # Run a MAFFT alignment to pad the seqs into equal length
    #     mafft(results_filename, output_path / "aligned.fasta")
    # # # After this go to variation_parser.
    # else:
    #     print("No hits found")


def align_seqs_globally():
    seq_a_file = Path("spike_uniprot.fasta")
    seq_b_file = Path("gisaid_results/translated.fasta")
    out_file = Path("gisaid_results/needle.fasta")
    emboss_needle(seq_a_file, seq_b_file, out_file)


def identity_calculation(path_to_needle: str) -> List:
    """ Calculate the identites for the seqs

        Returns: 
            result_record [list]
    """

    logging.info("Intializing trimming of the aligned seqs...")
    result_record = trim_seqs(path_to_needle)
    return result_record


def trim_seqs(filepath: Union[Path, str]) -> List:
    """ Parse the Needle file and look for highly identical matches

        Args: 
            filepath: File path to msf file

        Returns: 
            result_record: List of records of best matching seqs
    """
    needle_record = list(AlignIO.parse(str(filepath), "msf"))

    result_record = []

    for rec in needle_record[1:]:
        reference_seq = rec[1]

        seq_parser = IdenticalSequencesParser(reference_seq, rec[0])

        result = seq_parser.highly_identical_seqs()

        if result:
            result_record.append(result)

    return result_record


def write_records_to_file(result_records: List, filename: str):
    """ Write records to file 

        Args: 
            result_records: List of SeqRecords
            filename: File path
    """
    with open(filename, "w") as file:
        for rec in result_records:
            SeqIO.write(rec, file, "fasta")
    logging.info("Trimmed seqs written to file.")


if __name__ == '__main__':
    main()
