from Bio.Emboss.Applications import NeedleallCommandline
import subprocess


def emboss_needle(seq_a_file: str, seq_b_file: str, out_file: str):
    """ Do a global pairwise alignment using EMBOSS

        Args: 
            seq_a_file: First sequence
            seq_b_file: second sequence
            out_file: Output file

        Returns: 
            r [subprocess object]: Execute the commandline command for EMBOSS
        
    """
    needle_cline = NeedleallCommandline(asequence=seq_a_file,
                                        bsequence=seq_b_file,
                                        outfile=out_file,
                                        verbose=True,
                                        gapextend=1,
                                        gapopen=10)

    cmd = str(needle_cline)
    cmd = cmd.split(" ")
    cmd.append("-aformat=msf")
    print(cmd)

    r = subprocess.Popen(cmd)
    if r.communicate():
        print("Global alignment done.")


def main():
    seq_a = "/home/nadzhou/SEQs/spike_uniprot.fasta"
    seq_b = "/home/nadzhou/SEQs/CoV/polished_seqs/gisaid5-6_polished.fasta"
    out_file = "/home/nadzhou/SEQs/needle.fasta"

    emboss_needle(seq_a, seq_b, out_file)


if __name__ == '__main__':
    main()
