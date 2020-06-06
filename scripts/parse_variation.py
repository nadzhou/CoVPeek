from __future__ import annotations
from typing import List, Dict

import matplotlib.pyplot as plt
from matplotlib import figure
import numpy as np
import seaborn as sns
from Bio import SeqIO
from collections import Counter

from scipy.stats import chi2_contingency

sns.set()


class DivergenceParser:
    """Parse the globally aligned sequences to look for divergence
    """

    def __init__(self, sequences: List[List[str]] = None):
        """
        Args:
            sequences: List of sequences represented by list of aminoacids
        """
        if sequences is None:
            sequences = []
        self.npseqs = self._seq2np(sequences)
        self._conservation_scores = None

    @classmethod
    def retrieve_sequence(cls, aligned_path: str) -> DivergenceParser:
        """ Extract the protein sequence from the FASTA file

            Args
                aligned_path: Path to a fasta file with aligned sequences
                
            Returns
                seqs: A list of sequences represented by lists of aminoacids
        """
        records = list(SeqIO.parse(aligned_path, "fasta"))
        seqs = [[aa for aa in record] for record in records]
        new_parser = cls(seqs)
        return new_parser

    @staticmethod
    def _seq2np(seq: List[List[str]]) -> np.ndarray:
        """ Turn the sequence into numpy S1 array for calculations later.

            Args:
                seq [2d list]: List of lists that contain sequences

            Returns:
                np array [2d np array]: Np array that turns the chars into bytes
        """

        return np.asarray(seq, dtype='S1')

    @staticmethod
    def _normalize_data(ent_list: np.ndarray) -> np.ndarray:
        """ Takes the entropy array and normalizes the data.

            Args:
                ent_list [Nd array]: Entropy float array

            Returns:
                Normalized list [nd array]: Values between -1 and 1
        """
        return -(ent_list - np.mean(ent_list, axis=0)) / np.std(ent_list, axis=0)


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
        total_aminoacids = sum(aa_count.values())
        for aminoacid, count in aa_count.items():
            # and +1 for deletions
            pA *= (count / total_aminoacids)
        return np.sum(pA * np.log2(pA))

    def conservation_scores(self) -> np.ndarray:
        """ Calculate the normalized Shannon Entropy vertically
            for each position in the amino acid msa sequence.

            Returns:
                np apply array [nd float array]: Calculate conservation
                scores vertically into a float nd array
        """
        # caching
        if self._conservation_scores is None:
            conservation_scores = np.apply_along_axis(self._shannon, 0, self.npseqs)

            conservation_scores = self._normalize_data(conservation_scores)
            self._conservation_scores = conservation_scores
        return self._conservation_scores


    def plot_variation(self) -> figure:
        """
            Returns:
                fig: Variation score lineplot for a given sequence
        """
        scores = self.conservation_scores()
        aa_positions = np.arange(1, len(scores) + 1)
        # Now plot the data
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot()
        sns.lineplot(x=aa_positions, y=scores, color="purple", ax=ax)
        ax.set_xlabel("Amino acid position", fontsize=14)
        ax.set_ylabel("Variation score", fontsize=14)
        ax.set_title("Variation in the CoV spike protein", fontsize=14)

        minimum_labeled_variance = 5

        for position, score, aminoacids in zip(aa_positions, scores, self.npseqs.T):
            if score > minimum_labeled_variance:
                aminoacids = self.npseqs[:, position-1].astype(str)
                different_aas = set(aminoacids)
                ax.text(position, score, different_aas, color='b')
        return fig

    def aminoacids_in_variable_positions(self, minimum_var=5) -> Dict[int, List[str]]:
        variation_scores = self.conservation_scores()
        positions = {}
        assert len(variation_scores) == self.npseqs.shape[1]
        for i, score in enumerate(variation_scores):
            if score > minimum_var:
                aminoacids = self.npseqs[:, i].astype(str)
                positions[i + 1] = sorted(aminoacids)
        return positions


def main():
    """
        Loads the aligned FASTA file, plots the variation score for each AA position,
        prints aminoacids in positions with score > 5, and makes a countplot for position 614
    """
    aligned_path = "../operations/gisaid_results//aligned.fasta"
    # Call an instance of the class that converts it to ndarray then run through the functions, calculate entropy etc
    parser = DivergenceParser.retrieve_sequence(aligned_path)
    lineplot = parser.plot_variation()
    #plt.show()
    variable_positions = parser.aminoacids_in_variable_positions()
    print_variable(variable_positions)
    # countplots = make_countplots(variable_positions[614])
    plt.show()


def print_variable(aa_positions: Dict[int, List[str]], min_score: int = 5):
    """ Print amino acid and their counts at the specific position 
        Args:
            aa_positions: Dictionary of positions and aminoacids present in them
            min_score: Minimum variance score, default 5
    """
    print(f"Positions with variation score > {min_score}")
    for position, aminoacids in aa_positions.items():
        print(f"Position {position}: {dict(Counter(aminoacids))}")


def make_countplots(*aa_in_positions: List[str]) -> figure:
    """ Makes a countplot for each list of aminoacids

        Args:
            *aa_in_positions: Variable number of aminoacid 
            lists, plots 1 subplot for each list

        Returns:
            axes: A figure with several countplots
    """
    numplots = len(aa_in_positions)
    fig, axes = plt.subplots(numplots, 1, figsize=(12, 5*numplots))
    fig.suptitle("Amino acid frequency in different CoV patients")
    if aa_in_positions:
        # the function returns either a single ax or 2d ndarray of axes
        if numplots > 1:
            axes = axes.flat
        else:
            axes = [axes]
        for aa_in_position, ax in zip(aa_in_positions, axes):
            sns.countplot(aa_in_position, ax=ax)
            ax.set_ylabel("Amino acid frequency")
        return fig


if __name__ == '__main__':
    main()
