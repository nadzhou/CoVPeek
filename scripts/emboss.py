from Bio.Emboss.Applications import NeedleCommandline
import subprocess
import sys

from io import StringIO
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleallCommandline

from Bio import AlignIO




def emboss_needle(seq_a, seq_b, out_file): 
    """Do a global pairwise alignment using EMBOSS

    Args: 
        seq_a [str]: First sequence 
        seq_b [str]: second sequence
        out_file [str]: Output file 

    Returns: 
        r [subprocess object]: Execute the commandline command for EMBOSS
    
    """
    needle_cline = NeedleallCommandline(asequence=seq_a,
                            bsequence=seq_b, gapopen=10,
                            outfile=out_file, gapextend=1, verbose=True)


    cmd = str(needle_cline)
    cmd = cmd.split(" ")
    #cmd.append("-ossingle")
    #cmd.append("-sformat1=phylip")
    cmd.append("-aformat=msf")
    #cmd.append("-ossingle")
    print(cmd)
    r = subprocess.Popen(cmd)

    return r.communicate()






def main(): 
    seq_a = "/home/nadzhou/SEQs/spike_uniprot.fasta"
    seq_b = "/home/nadzhou/SEQs/CoV/polished_seqs/gisaid5-6_polished.fasta"
    out_file = "/home/nadzhou/SEQs/needle.fasta"

    emboss_needle(seq_a, seq_b, out_file)


if __name__ == '__main__': 
    main()