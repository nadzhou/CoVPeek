import argparse as ap

def parse_arguments(parser=None): 
    """
        Parser object for genome path, Uniprot, and output directory input. 
        
        Args: 
            parser [argparse]: Input argument. 

        Returns: \
            args [argparse]: Return the args.
    """
    if not parser: 
        parser = ap.ArgumentParser()

    parser.add_argument("cov_genome_path", 
                        help="Path to the CoV reference genome")
    
    parser.add_argument("uniprot_refseq_path", 
                        help="Path to the UniProt canonical protein")

    parser.add_argument("output_directory", 
                        help="An output directory where the code results go.")
                                         
    args = parser.parse_args()

    return args