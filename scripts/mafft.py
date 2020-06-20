from Bio.Align.Applications import MafftCommandline


def mafft(in_file: str, out_file: str):
    """ MAFFT command line for MSA. 

        Args: 
            in_file [str]: Input file 
    """

    mafft_cline = MafftCommandline(input=in_file)

    stdout, stderr = mafft_cline()
    print(mafft_cline)

    with open(out_file, "w") as file:
        file.write(stdout)


def main():
    in_file = "/home/nadzhou/Desktop/results.fasta"
    out_file = "../localdata/mafft_results.fasta"
    mafft(in_file, out_file)


if __name__ == '__main__':
    main()
