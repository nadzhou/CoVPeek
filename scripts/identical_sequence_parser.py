import numpy as np
from pathlib import Path

from Bio.SeqRecord import SeqRecord


class IdenticalSequencesParser: 
    """ Match reference sequence to target and calculate identity then write to file
    """

    def __init__(self, ref_seq, tar_seq): 
        self.ref_seq = ref_seq
        self.tar_seq = tar_seq
        self.ref_seq_len = len(self.ref_seq)
        self.np_seq = self.seqs2np()


    @staticmethod
    def seqs2np() -> np.ndarray:
        """ Turn the sequence into numpy S1 array for calculations later.

            Returns:
                np array [2d np array]: Np array that turns the chars into bytes
        """
        return np.asarray((self.ref_seq.seq, self.tar_seq.seq))


    def identical_freq_calc(self) -> np.ndarray:
        """ Calculate percent identity of given two sequences
        """

        identical_aa_freq = np.apply_along_axis(
                                    _identity_calc, 0, self.np_seq)
        self.identical_aa_freq = identical_aa_freq.astype(np.int64)


    def aa_hit_calc(self): 
        """ Isolate the target sequence only where it is identical 
            to reference sequence
        """
        trimd_tar_seq = np.apply_along_axis(
                        _trim_target_seq, 0, self.np_seq)
        self.trimd_tar_seq = trimd_tar_seq.astype(np.str)


    @staticmethod
    def _identity_calc(array) -> np.ndarray:
        """ Calculate idenity if both sequences have a match identical_aa_freq it up.

            Args:
                array [nd array]: Pair of protein sequence characters
        """
        identical_aa_freq = 0

        if array[0] == array[1] and array[0] != "-":
            identical_aa_freq += 1

        return identical_aa_freq


    @staticmethod
    def _trim_target_seq(array) -> np.ndarray: 
        """ Trim target sequence wherever identical with reference sequence.

            Args:
                array [nd array]: Pair of protein sequence characters
        """
        hit = ""

        if array[0] == array[1] and array[0] != "-":
            hit = array[1]

        return identical_aa_freq, hit


    def count_checker(self):
        """ Check whether the trimmed sequence has significant idenity
            and length, if it does, write the sequences to file. 
        """    
        if not np.all(self.identical_aa_freq == 0):
            self.trimd_tar_seq = "".join(item for item in trimd_tar_seq)

            self.identity_score = np.true_divide(sum(identical_aa_freq), 
                                                        self.ref_seq_len)

            return self.seq2record()


    def seq2record(self) -> SeqRecord:
        """ Write the sequences that have identity scores greater than 80 percent.
        """
        try: 
            if self.identity_score > 0.7 and len(self.trimd_tar_seq) > 50:
                print(f"hit hit {identity_score}")

                target_seq_record = SeqRecord(Seq(self.trimd_tar_seq),
                                        id=self.tar_seq.id,
                                        name=self.tar_seq.name,
                                        description=self.tar_seq.description)

                return target_seq_record

