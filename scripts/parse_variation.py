from typing import List

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio import SeqIO
from collections import Counter

sns.set()


class DivergenceParser:
    """Parse the globally aligned sequences to look for divergence
    """
    def __init__(self, sequences: List[List[str]]):
        """
        Args:
            sequences: List of sequences represented by list of aminoacids
        """
        self.seqs: List[List[str]] = sequences
        self.npseqs: np.ndarray = self.seq2np(sequences)

    @staticmethod
    def seq2np(seq: List[List[str]]) -> np.ndarray:
        """"Turn the sequence into numpy S1 array for calculations later.

            Args:
                seq [2d list]: List of lists that contain sequences

            Returns:
                np array [2d np array]: Np array that turns the chars into bytes

        """

        return np.asarray(seq, dtype='S1')

    @staticmethod
    def normalize_data(ent_list: np.ndarray) -> np.ndarray:
        """ Takes the entropy array and normalizes the data.

            Args:
                ent_list [Nd array]: Entropy float array

            Returns:
                Normalized list [nd array]: Values between -1 and 1

        """
        normalized = -(ent_list - np.mean(ent_list, axis=0)) / np.std(ent_list, axis=0)
        return normalized

    @staticmethod
    def _shannon(array: np.ndarray) -> np.ndarray:
        """Calculate Shannon Entropy vertically via loop.

            Args:
                array [nd array]: 1d array of sequences from all the species

            Returns:
                entropy [nd float array]: Per column value for each position vertically

        """

        aa_count = Counter(array)

        pA = 1
        for k, v in aa_count.items():
            pA *= (v / 21)

        return -np.sum(pA * np.log2(pA))

    def conservation_score(self, normalized: bool):
        """Calculate the Shannon Entropy vertically
            for each position in the amino acid msa sequence.

            Args:
                normalized - if conservation scored should be normalized

            Returns:
                np apply array [nd float array]: Calculate conservation
                scores vertically into a float nd array
            
        """
        conservation_score = np.apply_along_axis(self._shannon, 0, self.npseqs)
        if normalized:
            conservation_score = self.normalize_data(conservation_score)
        return conservation_score

    def moving_average(self, data, n=3):
        """Calculated the rolling average of teh data

            Args:
                data [numpy nd array]: Float array of the conservation score calculated
                n [int]: Rolling average wegith

            Returns:
                avg_data [numpy nd array]: Float array of rolling average

        """

        avg_data = np.cumsum(data, dtype=float)
        avg_data[n:] = avg_data[n:] - avg_data[:-n]

        return avg_data[n - 1:] / n

    def label_plot(self, norm_list, norm_list_len, val, ax):
        """Label the amino acids that are diverging from the aligned sequences
        """

        a = np.concatenate({'x': norm_list_len, 'y': norm_list, 'val': val}, axis=1)
        for i, point in a.iteritems():
            ax.text(point['x'] + .02, point['y'], str(point['val']))

    def aminoacids_in_variable_positions(self, minimum_var=5):
        variation_score = self.conservation_score(normalized=True)
        positions = {}
        for index, (variation_score, aminoacids) in enumerate(zip(variation_score, self.npseqs)):
            if variation_score > minimum_var:
                sorted_aa = sorted((aminoacids.astype(str)))
                positions[index + 1] = sorted_aa
        return positions


def main():
    """
    Loads the aligned FASTA file and plots the variation score for each AA position
    """
    aligned_path = "../notebooks/gisaid_results/mafft_aligned.fasta"
    aa_counts = plot_variation(aligned_path)

    fig = make_countplot(aa_counts, aa_counts)
    fig.show()

def retrieve_sequence(aligned_path: str) -> List[List[str]]:
    """Extract the protein sequence from the FASTA file
        Args
            aligned_path: Path to a fasta file with aligned sequences
        Returns
            seqs: A list of sequences represented by lists of aminoacids
    """
    records = list(SeqIO.parse(aligned_path, "fasta"))
    seqs = [[aa for aa in record] for record in records]
    return seqs


def plot_variation(aligned_path: str):
    """Make an instance of the DivergenceParser on a given FASTA file 
        and then plot the results. 

        Args: 
            aligned_fasta_file [str]: Address of the FASTA file

        Returns: 
            aa_count [list]: Frequency of amino acids at the most divergent position

    """
    sequences = retrieve_sequence(aligned_path)
    # Call an instance of the class that converts it to ndarray then run through the functions, calculate entropy etc.
    parser = DivergenceParser(sequences)
    norm_list = parser.conservation_score(normalized=True)
    print("Positions with variation score > 5")
    variable_positions = parser.aminoacids_in_variable_positions()
    for position, aminoacids in variable_positions.items():
        print(f"Position {position}: {dict(Counter(aminoacids))}")
    aa_counts = plot_data(norm_list, parser.npseqs)
    return aa_counts


def plot_data(norm_list: np.ndarray, np_seq: np.ndarray) -> List[str]:
    """
    Args:
        norm_list:
        np_seq:

    Returns:
        aa_count [list]: Frequency of amino acids at the most divergent position
    """
    norm_list_len = np.arange(len(norm_list))

    # Now plot the data
    plt.figure(figsize=(15, 10))
    ax = sns.lineplot(norm_list_len, norm_list, color="purple")

    plt.xlabel("Amino acid position", fontsize=14)
    plt.title("Variation in the CoV spike protein", fontsize=14)
    plt.ylabel("Variation score", fontsize=14)

    # i can change the 0-based index to 1-based position, if needed
    for index, variation_score, aminoacids in zip(norm_list_len, norm_list, np_seq):
        if variation_score > 5:
            sorted_aa = sorted((aminoacids.astype(str)))
            aa_s = set(sorted_aa)
            ax.text(index, variation_score, len(aa_s), color='b')
    plt.show()


def make_countplot(aa_count, aa_count2):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    sns.countplot(aa_count, ax=ax1)
    ax1.set_ylabel("Amino acid frequency")
    ax1.set_title("Amino acid ferquency in different CoV patients")

    sns.countplot(aa_count2, ax=ax2)
    ax2.set_ylabel("Amino acid frequency")
    return fig


if __name__ == '__main__':
    main()
