
from Bio import SeqIO
from Bio import AlignIO

def trim_seqs(filepath: str):
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

        print(f"reference: {reference_seq.id}")
        print(f"target {rec[0].id}")



path = "../operations/gisaid_results/needl2.fasta"

trim_seqs(path)