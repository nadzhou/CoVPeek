import numpy as np
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

from Bio.Seq import Seq

from scipy import chi2_contingency

class IdenticalSequencesParser: 
    """ Match reference sequence to target and calculate identity then write to file
    """

    def __init__(self, ref_seq, tar_seq): 
        self.ref_seq = ref_seq
        self.tar_seq = tar_seq
        self.ref_seq_len = len(self.ref_seq)
        self.np_seq = self.seqs2np()
        self._highly_identical_seqs = None


    def seqs2np(self) -> np.ndarray:
        """ Turn the sequence into numpy S1 array for calculations later.

            Returns:
                np array [2d np array]: Np array that turns the chars into bytes
        """
        return np.asarray((self.ref_seq.seq, self.tar_seq.seq))


    def aa_hit_calc(self): 
        """ Isolate the target sequence only where it is identical 
            to reference sequence
        """
        trimd_tar_seq = np.apply_along_axis(
                        self._trim_target_seq, 0, self.np_seq)

        self.trimd_tar_seq = trimd_tar_seq.astype(np.str)


    def identical_freq_calc(self) -> np.ndarray:
        """ Calculate percent identity of given two sequences
        """
        identical_aa_freq = self.trimd_tar_seq.view(np.uint8)

        self.identical_aa_freq = np.where(identical_aa_freq > 1, 1, 0)


    @staticmethod
    def _trim_target_seq(array) -> np.ndarray: 
        """ Trim target sequence wherever identical with reference sequence.

            Args:
                array [nd array]: Pair of protein sequence characters
        """
        hit = ""

        if array[0] == array[1] and array[0] != "-":
            hit = array[1]

        return hit


    def count_checker(self):
        """ Check whether the trimmed sequence has significant idenity
            and length, if it does, write the sequences to file. 
        """    
        if not np.all(self.identical_aa_freq == 0):
            self.trimd_tar_seq = "".join(item for item in self.trimd_tar_seq)

            self.identity_score = np.true_divide(sum(self.identical_aa_freq), 
                                                        self.ref_seq_len)

            return self.seq2record()


    def seq2record(self) -> SeqRecord:
        """ Write the sequences that have identity scores greater than 80 percent.
        """

        if self.identity_score > 0.7 and len(self.trimd_tar_seq) > 50:
            print(f"hit hit {self.identity_score}")

            target_seq_record = SeqRecord(Seq(self.trimd_tar_seq),
                                            id=self.tar_seq.id,
                                            name=self.tar_seq.name,
                                            description=self.tar_seq.description)

            return target_seq_record

    def highly_identical_seqs(self): 
        """ Calculate identity and trim before giving out SeqRecord as output

            Returns: 
                _highly_identical_seqs [SeqRecord]
        """
        if self._highly_identical_seqs == None: 
            self.aa_hit_calc()      
            self.identical_freq_calc()
            self._highly_identical_seqs = self.count_checker()

        return self._highly_identical_seqs