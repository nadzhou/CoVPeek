import numpy as np
from pathlib import Path

from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess
from Bio.SubsMat import MatrixInfo
from Bio import SeqIO
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


def extract_dna(genome_records,
                uniprot_record, out_file_path):

    out_file_path.mkdir(parents=True, exist_ok=True)

    records = [uniprot_record]

    for genome_record, num in zip(genome_records, range(len(genome_records) + 1)): 
        genome_record.seq = Seq(re.sub(r'(\W)', '', str(genome_record.seq)))
        
        print(f"Executing job {num}")
        
        saved_records = []

        orf_proteins = find_orfs_with_trans(genome_record.seq)

        for protein, i in zip(orf_proteins, range(len(orf_proteins) + 1)):
            translated_pt_record = SeqRecord(Seq(protein), id=f"{num}gen_{i}_orf_{genome_record.id}", 
                                            name=genome_record.name,
                                            description=genome_record.description)

            saved_records.append(translated_pt_record)

        if len(saved_records) > 2: 
            records.append(saved_records)
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



def genome_retrieve(genome_file_path):
    """Retrieve a given genome from a FASTA file
    
    Args: 
        genome_file_path [str]: Address to the genome sequence

    Returns record [seqrecord object]: Returns the genome as a SeqRecord object
    """
    return SeqIO.read(genome_file_path, "fasta")


def genome_orf_find(genome_record,
                    table=11, min_pro_len=100):

    """Find the open reading frames from a given genome record

    Args: 
        genome_record [seqrecord object]: CoV2 genome record
        table [int]: Look up table for translation
        min_pro_len [int]: Minimum protein length

    
    Returns: 
        find_orfs [function]: Function that finds the open reading frames 
        in a given genome record. 
    """

    return find_orfs_with_trans(genome_record.seq,
                                table, min_pro_len)




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

        seq_record2 = SeqRecord(Seq(trimmed_seq2), id=whole_seq2.id, 
                                description=whole_seq2.description, name=whole_seq2.name)

        return seq_record2




def seq_writer(): 
    out_file = Path("/home/nadzhou/SEQs/CoV/polished_seqs/")
    genome_record, uniprot_record = retrieve_canonical_seqs()

    genomes = extract_dna(genome_record, uniprot_record, out_file)

    with open(out_file/"gisaid5-6_polished.fasta", "w") as file: 
        for genome in genomes: 
            SeqIO.write(genome, file, "fasta")


def retrieve_canonical_seqs(): 
    """Retrieve the GISAID and Uniprot datasets
    
    
    Returns: 
        genome_record [seqrecord object]: GISAID data
        uniprot_record [seqrecord object]: Uniprot record
    """

    genome_path = "/home/nadzhou/SEQs/CoV/raw_seqs/gisaid5-6.fasta"
    genome_record = list(SeqIO.parse(genome_path, "fasta"))

    uniprot_path = "/home/nadzhou/SEQs/CoV/raw_seqs/spike_uniprot.fasta"
    uniprot_record = SeqIO.read(uniprot_path, "fasta")

    return genome_record, uniprot_record








def main():
    out_file = Path("/home/nadzhou/SEQs/")


    # uniprot_path = "/home/nadzhou/ml/CoV/raw_seqs/spike_uniprot.fasta"
    # uniprot_record = SeqIO.read(uniprot_path, "fasta")
    # canonical_len = len(uniprot_record.seq)

    # # msa_dir = Path("/home/nadzhou/Desktop/polished")
    # record = list(AlignIO.read("/home/nadzhou/Desktop/divided/out.fasta", "fasta"))
    # # result_record = [uniprot_record]
    # # for file in msa_dir .iterdir(): 
    # #     record = list(SeqIO.parse(file, "fasta"))s
    # #     canonical_ = record[0]

    record = list(AlignIO.parse(out_file/"needlaa", "msf"))
    result_record = []

    for rec, num in zip(record[1:], range(len(record[1:]))): 
        canonical_ = rec[0]
        canonical_len = len(canonical_.seq)


        print(canonical_.name, rec[1].name)
        result = percent_id_calc(canonical_, rec[1], f"{rec[0].id[:7]}_{num}", canonical_len)

        if result: 
            result_record.append(result)


    if result_record: 
        with open(out_file/"5-6_MAFFT0.fasta", "w") as file: 
            for rec in result_record:
                SeqIO.write(rec, file, "fasta")




if __name__ == '__main__':
    main()


    # out_file = Path("/home/nadzhou/Desktop")
    # genome_record, uniprot_record = retrieve_canonical_seqs()

    # records = extract_dna(genome_record, uniprot_record, out_file)

    # with open(f"{out_file}/results.fasta", "w") as file: 
    #     for record in range( 0, len(records), 500): 
    #         print(record)
    #         SeqIO.write(record, file, "fasta")

    #extract_dna(genome_record, uniprot_record, out_file)

    # translated_seq_path = Path("/home/nadzhou/DEVELOPMENT/CoV/scripts/tmp/gisaid.fasta")
    # uniprot_path = "/home/nadzhou/ml/CoV/raw_seqs/spike_uniprot.fasta"

