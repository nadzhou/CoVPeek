import re
from typing import List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from pathlib import Path
import logging

from multiprocessing import Pool


def main():
    logging.basicConfig(level=logging.INFO)
    genome = Path('../localdata/gisaid_partial.fasta')
    result = Path('../localdata/translated_partial.fasta')
    translate_genome(genome_path=genome, result_file=result)


def translate_genome(genome_path: Path, result_file: Path):
    """Translate the GISAID genome into ORFs.
    """
    genome_records = read_genome(genome_path)
    results = extract_dna(genome_records)
    with result_file.open("w") as file:
        for result in results:
            SeqIO.write(result, file, "fasta")


def read_genome(genome_path: Path) -> List[SeqRecord]:
    return list(SeqIO.parse(str(genome_path), "fasta", alphabet=IUPAC.ambiguous_dna))


def extract_dna(genome_records: List[SeqRecord]) -> List[List[SeqRecord]]:
    """ Extract GISAID genome DNA and then translate each record
        into ORFs, and output the list.

        Args:
            genome_records: GISAID genomes in a SeqRecord

        Returns:
            records: Translated ORFs in a list
    """
    logging.info("Initiating translation...")

    with Pool() as pool:
        results = pool.map(translate_record, enumerate(genome_records, start=1), chunksize=25)
        records = [record for record in results if len(record) > 2]

    logging.info("Genome record translated.")
    return records


def translate_record(numbered_record: Tuple[int, SeqRecord]):
    num, genome_record = numbered_record
    # i use tuple because of map
    genome_record.seq = remove_non_chars_from_seq(genome_record.seq)
    saved_records = []
    orf_proteins = find_orfs_with_trans(genome_record.seq)
    patient_country = find_patient_country(genome_record.id)
    logging.debug(patient_country)
    for i, protein in enumerate(orf_proteins, start=1):
        translated_pt_record = make_seqrecord(protein,
                                              patient_country,
                                              i,
                                              num)

        saved_records.append(translated_pt_record)
    return saved_records


def remove_non_chars_from_seq(seq: Seq) -> Seq:
    return Seq(re.sub(r'(\W)', '', str(seq)))
    # TODO - get rid of it


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
                    answer.append(trans[aa_start:aa_end])
                aa_start = aa_end + 1
                
    return answer


def pad_seq(seq: Seq) -> Seq:
    """ Pad seq to multiple of 3 with N

        Args:
            seq [str]: Amino acid seq

        Returns:
            seq [str]: Padded seq
    """
    remainder = len(seq) % 3
    padded = seq if remainder == 0 else seq + 'N' * (3 - remainder)
    assert len(padded) % 3 == 0
    return padded


def find_patient_country(genome_id: str):
    match = re.search(r'19/(\w*)/', genome_id)
    return match.group(1) if match else ''


def make_seqrecord(protein, patient_country, i, num) -> SeqRecord:
    return SeqRecord(Seq(protein),
                     id=f"{num}-{i}-{patient_country}",
                     name=f"{num}-{i}-{patient_country}",
                     description=f"{num}-{i}-{patient_country}")


def custom_translation(sequence) -> str:
    pass


if __name__ == '__main__':
    main()
