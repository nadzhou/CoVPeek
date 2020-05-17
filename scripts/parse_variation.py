import numpy as np
import pandas as pd
from pathlib import Path

from Bio import AlignIO
from Bio import SeqIO

import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

from collections import Counter

import sys
import numpy
numpy.set_printoptions(threshold=sys.maxsize)


def retrieve_dna(aligned_path): 
    records = list(SeqIO.parse(aligned_path, "fasta"))


    seqs = [[x for x in y] for y in records]



    return seqs

def seq2np(seq): 
    """"Turn the sequence into numpy S1 array for calculations later. 
    
    Args: 
        seq [2d list]: List of lists that contain sequences 
    
    Returns: 
        np array [2d np array]: Np array that turns the chars into bytes
        
    """
    
    return np.asarray(seq, dtype='S1')



def normalize_data(ent_list): 
    """Takes the entropy array and normalizes the data. 
    
    Args: 
        ent_list [Nd array]: Entropy float array
        
    Returns: 
        Normalized list [nd array]: Values between -1 and 1
        
    """
    
    return -(ent_list - np.mean(ent_list, axis=0)) / np.std(ent_list, axis=0)  





def _shannon( array): 
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


    
    
    
def conservation_score(np_seq): 
    """Calculate the Shannon Entropy vertically
    for each position in the amino acid msa sequence.
    
    Args: 
        np_seq [Numpy nd array]: Np array of sequences
        
    Returns: 
        np apply array [nd float array]: Calculate conservation 
        scores vertically into a float nd array   
    """
    
    return np.apply_along_axis(_shannon, 0, np_seq)  



def moving_average(data, n=3) :
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


def label_plot(norm_list, norm_list_len, val, ax): 
    a = np.concatenate({'x': norm_list_len, 'y': norm_list, 'val': val}, axis=1)
    for i, point in a.iteritems(): 
        ax.text(point['x']+.02, point['y'], str(point['val']))




def main(): 

    aligned_path = Path("/home/nadzhou/SEQs/5-6.fasta")


    seqs = retrieve_dna(aligned_path)

    np_seq = seq2np(seqs)

    ent_list = conservation_score(np_seq)
    norm_list = normalize_data(ent_list)
    #norm_list = moving_average(norm_list)

    norm_list_len = np.arange(len(norm_list))


    ax=sns.lineplot (norm_list_len, norm_list, color="purple")

    plt.xlabel("Amino acid position")
    plt.title("Variation in the CoV spike protein")
    plt.ylabel("Variation score")

    for x, y, name in zip(norm_list_len, norm_list, np_seq): 
        if y > 5: 
            print(len(set(name)))
            ax.text(x,y, ", ".join(str(name[:3], encoding='utf-8')), color='b')
    plt.show()




if __name__ == '__main__': 
    main()