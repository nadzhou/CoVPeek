import numpy as np
import pandas as pd
from pathlib import Path

from Bio import AlignIO
from Bio import SeqIO

import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

from collections import Counter

import sys

class DivergenceParser: 
    """Parse the globally aligned sequences to look for divergence"""

    def retrieve_dna(self, aligned_path): 
        records = list(SeqIO.parse(aligned_path, "fasta"))


        seqs = [[x for x in y] for y in records]

        return seqs

    def seq2np(self, seq): 
        """"Turn the sequence into numpy S1 array for calculations later. 

        Args: 
            seq [2d list]: List of lists that contain sequences 

        Returns: 
            np array [2d np array]: Np array that turns the chars into bytes

        """

        return np.asarray(seq, dtype='S1')



    def normalize_data(self, ent_list): 
        """Takes the entropy array and normalizes the data. 

        Args: 
            ent_list [Nd array]: Entropy float array

        Returns: 
            Normalized list [nd array]: Values between -1 and 1

        """

        return -(ent_list - np.mean(ent_list, axis=0)) / np.std(ent_list, axis=0)  





    def _shannon(self, array): 
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

        return -np.sum(pA*np.log2(pA))


    def conservation_score(self, np_seq): 
        """Calculate the Shannon Entropy vertically
        for each position in the amino acid msa sequence.

        Args: 
            np_seq [Numpy nd array]: Np array of sequences

        Returns: 
            np apply array [nd float array]: Calculate conservation 
            scores vertically into a float nd array   
        """

        return np.apply_along_axis(self._shannon, 0, np_seq)  



    def label_plot(self, norm_list, norm_list_len, val, ax): 
        a = np.concatenate({'x': norm_list_len, 'y': norm_list, 'val': val}, axis=1)
        
        for i, point in a.iteritems(): 
            ax.text(point['x']+.02, point['y'], str(point['val']))




def main(): 
    # Load the aligned FASTA file
    aligned_path = "gisaid_results/aligned.fasta"

    # Call an instance of the class
    inst = DivergenceParser()
    seqs = inst.retrieve_dna(aligned_path)

    # Convert to nd array then run through 
    # the functions, calculate entropy etc.
    np_seq = inst.seq2np(seqs)

    ent_list = inst.conservation_score(np_seq)
    norm_list = inst.normalize_data(ent_list)

    norm_list_len = np.arange(len(norm_list))


    # Now plot the data
    plt.figure(figsize=(15, 10))
    ax=sns.lineplot (norm_list_len, norm_list, color="purple")

    plt.xlabel("Amino acid position", fontsize=15)
    plt.title("Variation in the CoV spike protein", fontsize=14)
    plt.ylabel("Variation score", fontsize=14)

    aa_count = []

    for x, y, name in zip(norm_list_len, norm_list, np_seq): 
        if y > 5: 
            aa_s = sorted(name.astype(str))
            aa_count = aa_s
            aa_s = set(aa_s)
            print(aa_s)
            ax.text(x,y, len(aa_s), color='b')
    plt.show()


if __name__ == '__main__': 
    main()