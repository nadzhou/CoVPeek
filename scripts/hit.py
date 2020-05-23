import numpy as np


def _trim_target_seq(array) -> np.ndarray: 
    """ Trim target sequence wherever identical with reference sequence.

        Args:
            array [nd array]: Pair of protein sequence characters
    """
    hit = ""

    if array[0] == array[1] and array[0] != "-":
        hit = array[1]

    return hit


def aa_hit_calc(np_seq): 
    """ Isolate the target sequence only where it is identical 
        to reference sequence
    """
    trimd_tar_seq = np.apply_along_axis(
                    _trim_target_seq, 0, np_seq)

    print(trimd_tar_seq)
    trimd_tar_seq = trimd_tar_seq.view(np.uint8)
    print(trimd_tar_seq)
    trimd_tar_seq = np.where(trimd_tar_seq > 1, 1, 0)
    print(trimd_tar_seq)


seq = np.array((['A', "C", "G", "T"], ["A", "G", "G", "T"]))


aa_hit_calc(seq)