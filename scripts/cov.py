#! usr/bin/env python3

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SubsMat import MatrixInfo
from Bio import SeqIO

import numpy as np
from pathlib import Path
import subprocess


from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO

from alignment.sequence import Sequence
from alignment.vocabulary import Vocabulary
from alignment.sequencealigner import SimpleScoring, LocalSequenceAligner


from Bio import AlignIO
import tables
from Bio import Align

from emboss import emboss_needle
from io import StringIO

import re
from mafft import mafft

def extract_dna(genome_records, out_file_path):

    out_file_path.mkdir(parents=True, 
                        exist_ok=True)

    records = []

    print("Initiating translation...")
    
    for genome_record, num in zip(genome_records, 
                                range(len(genome_records) + 1)): 

        genome_record.seq = Seq(re.sub(r'(\W)', '', str(genome_record.seq)))
                
        saved_records = []

        orf_proteins = find_orfs_with_trans(genome_record.seq)

        for protein, i in zip(orf_proteins, range(len(orf_proteins) + 1)):
            translated_pt_record = SeqRecord(Seq(protein), id=f"{num}gen_{i}_orf_{genome_record.id}", 
                                            description=f"{num}gen_{i}_orf_{genome_record.description}",
                                            name=f"{num}gen_{i}_orf_{genome_record.name}")

            saved_records.append(translated_pt_record)

        if len(saved_records) > 2: 
            records.append(saved_records)
    print("Genome record translated.")
    return records





def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """

    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))




def find_orfs_with_trans(seq, trans_table=1, min_protein_length=100):
    """Code from the Biopython library for finding open reading frames 

    Args: 
        seq [str]: Protein sequence
        trans_table [int]: Look-up table 
        min_protein_length [int]: Minimum protein length

    Returns: 
        answer [list]: List of protein sequences with different reading frames

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




def _identity_calc(array):
    """Calculate idenity if both sequences have a match my_count it up. 
    
    Args: 
        array [nd array]: Pair of protein sequence characters

    Returns: 
        my_count [int]: Hit value, 1 for identical, 0 for else
        trimmed_seq1 [char]: Character from the first sequence
        seq2 [char]: Character formt he second sequence

    """
    my_count = 0
    seq2 = ""

    if array[0] == array[1] and array[0] != "-":
        seq2 = array[1]
        my_count += 1

    return my_count, seq2



def percent_id_calc(whole_seq1, whole_seq2, seq_tag, canonical_len):
    """Calculate percent identity of given two sequences
    
    Args: 
        whole_seq1 [seqrecord object]: First protein sequence
        whole_seq2 [seqrecord]: Second protein sequence
        seq_tag [int]: Tag for the protein sequence used 
                        for writing the file name
    
    """
    np_seq = np.asarray((whole_seq1.seq, whole_seq2.seq))

    my_count, trimmed_seq2 = np.apply_along_axis(_identity_calc, 0, np_seq)
    my_count = np.asarray(my_count, int)


    if not np.all(my_count==0): 
        trimmed_seq2 = "".join(item for item in trimmed_seq2.astype(str))
        identity_score = np.true_divide(sum(my_count), canonical_len)


        return write_trimmed_seqs(my_count, trimmed_seq2, whole_seq1, 
                                whole_seq2, identity_score, seq_tag)

    

def write_trimmed_seqs(my_count, trimmed_seq2, whole_seq1,
                        whole_seq2, identity_score, seq_tag): 

    """Write the sequences that have identity scores greater than 80 percent. 
    
    Args: 
        my_count [int]: Sum of identity scores 
        trimmed_seq1 [str]: First trimmed sequence
        trimmed_seq2 [str]: Second trimmed sequence
        whole_trimmed_seq1 [str]: First wholeprotein sequence
        whole_seq2 [str]: Second protein sequence

    """
    if identity_score > 0.7 and len(trimmed_seq2) > 50:
        print(f"hit hit {identity_score}")

        seq_record2 = SeqRecord(Seq(trimmed_seq2), 
                                    id=whole_seq2.id, 
                                    name=whole_seq2.name,
                                    description=whole_seq2.description)

        return seq_record2





def main():
# Retrieve GISAID genomes and translte according to ORFs
    genome_path = "/home/nadzhou/SEQs/CoV/raw_seqs/gisaid19.fasta"
    uniprot_record = SeqIO.read("spike_uniprot.fasta", "fasta")
    genome_record = list(SeqIO.parse(genome_path, "fasta"))

    gisaid_results_path = Path("gisaid_results/")

    results = extract_dna(genome_record, gisaid_results_path)

    with open(gisaid_results_path/"translated.fasta", "w") as file: 
        for record in results: 
            SeqIO.write(record, file, "fasta")

# Now globally align my sequences
    seq_a = "spike_uniprot.fasta"
    seq_b = "gisaid_results/translated.fasta"
    out_file = "gisaid_results/needle.fasta"

    emboss_needle(seq_a, seq_b, out_file)

# Now for the identitiy calculation bit. 
    record = list(AlignIO.parse("gisaid_results/needle.fasta", "msf"))
    result_record = []

    print("Intializing trimming of the aligned sequences...")
    for rec, num in zip(record[1:], range(len(record[1:]))): 
        canonical_ = rec[0]
        canonical_len = len(canonical_.seq)

        result = percent_id_calc(canonical_, rec[1], f"{rec[0].id[:7]}_{num}", canonical_len)

        if result: 
            result_record.append(result)

    if result_record: 
        with open("gisaid_results/trimmed_seqs.fasta", "w") as file: 
            for rec in result_record:
                SeqIO.write(rec, file, "fasta")
        print("Trimmed sequences written to file.")


# Run a MAFFT alignment to pad the sequences into equal length 
    mafft("gisaid_results/trimmed_seqs.fasta")

# After this go to variation_parser.py



if __name__ == '__main__':
    main()


