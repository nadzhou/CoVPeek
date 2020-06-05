from Bio.Align.Applications import MafftCommandline


def mafft(in_file: str):
    """ MAFFT command line for MSA. 

        Args: 
            in_file [str]: Input file 
    """

    mafft_cline = MafftCommandline(input=in_file)

    stdout, stderr = mafft_cline()
    print(mafft_cline)

    with open("../operations/gisaid_results/aligned.fasta", "w") as file:
        file.write(stdout)


def main():
    in_file = "/home/nadzhou/Desktop/results.fasta"
    mafft(in_file)


if __name__ == '__main__':
    main()
