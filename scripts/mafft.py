from Bio.Align.Applications import MafftCommandline





def mafft(in_file): 
    mafft_cline = MafftCommandline(input=in_file)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    stdout, stderr = mafft_cline()
    print(mafft_cline)

    with open("/home/nadzhou/Desktop/aligned.fasta", "w") as file: 
        file.write(stdout)





def main(): 
    in_file = "/home/nadzhou/Desktop/results.fasta"

    mafft(in_file)








if __name__ == '__main__': 
    main()