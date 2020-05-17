

def parse_percent_id_file(text_file_path):
    """Parse the MSA matrix file to get the sequence score and ID
    
    Args: 
        text_file_path [str]: File path 

    Returns: 
        hit_records [dict]: Sequence ID as key and position in the matrix
                            with maximum score as value
    """

    prot_data = ""

    ids_hits_dict = {}

    for file in text_file_path.iterdir():
        if "txt" in file.name:
            with open(file) as file:
                uniprot_identity_score = next(file)
                prot_data = next(file)
                prot_data = prot_data.split(" ")
                identity_score = prot_data[0]

                prot_data = [x for x in prot_data[1: ] if x]
                local_hit = [index for index, score in 
                                enumerate(prot_data) if float(score) >= 80]

                if local_hit: ids_hits_dict[id] = local_hit

    return ids_hits_dict





def align_high_score_seqs(origin_file_path, id_dict):
    """Find the sequences that have the maximum identity in a directory
    
    Args: 
        id_dict [dict]: Identity dictionary with Uniprot id as key, 
                        position of the best scoring sequence as value

        origin_file_path [str]: File path for the aligned MSA sequences

    Returns: 
        hit_records [list of lists]: SeqRecord for both the Uniprot and best scoring 
                                    translated protein per row, and different files
                                    per column
    """
    hit_records = []

    for file in origin_file_path.iterdir():
        for k, v in id_dict.items():
            if k in file.name:
                records = list(SeqIO.parse(file, "fasta"))

                v = list(v)

                if len(v) > 1: 
                    print(v)

                for index in v: 
                    if int(index) > 0: 
                        records_list = [records[x] for x in v]
                        hit_records.append(records_list)


    for i in hit_records:
        print(i)
        print(len(i))


    return hit_records





def uniprot_genome_extract():
    """Extract sequences from the whole directory as Uniprot record"""

    uniprot_path = Path("/home/nadzhou/ml/CoV/polished_seqs/")

    print("Retrieving CoV genome...")

    dna_path = Path("/home/nadzhou/ml/CoV/raw_seqs/GISAID_DNA.fasta")
    print("Finding ORFS...")

    genome_record = genome_retrieve(dna_path)
    orf_proteins = genome_orf_find(genome_record)
    out_file_path = Path("/home/nadzhou/DEVELOPMENT/CoV/translated_seqs/")


    for file in uniprot_path.iterdir():
        extract_dna(genome_record, orf_proteins, file, out_file_path)



def bulk_align():
    """Run Clustal Omega MSA on a whole directory"""

    in_file_path = Path("/home/nadzhou/DEVELOPMENT/CoV/translated_seqs/")
    out_file_path = Path("/home/nadzhou/DEVELOPMENT/CoV/msa_seqs")
    out_file_path.mkdir(parents=True, exist_ok=True)
    print("Starting Bulk MSA...")

    count = 0

    for file in in_file_path.iterdir():
        count += 1
        msa(file, f"{out_file_path}/{count}_msa_result.fasta", mat_file=f"glob_scores{count}.txt")









# def aligned_seqs(genome_file_path, uniprot_file_path):
#     aligner = Align.PairwiseAligner()
#     aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
#     aligner.mode = "global"
#     table = 11
#     min_pro_len = 100

#     uniprot_record = SeqIO.read(uniprot_file_path, "fasta")
#     genome_record = SeqIO.read(genome_file_path, "fasta")



#     for protein, num in zip(orf_proteins, range(len(orf_proteins) + 1)):

#         aln = aligner.align(Seq(protein), uniprot_record.seq)
#         top_aln = list(aln)
#         top_aln = top_aln[0]

#         print(top_aln)







    # records = []

    # for protein, num in zip(orf_proteins, range(len(orf_proteins) + 1)):
    #     translated_pt_record = SeqRecord(Seq(protein, IUPAC.protein),
    #                                     id=f"{num}_{genome_record.id}", name=genome_record.name,
    #                                     description=genome_record.description)

    #     records.append(translated_pt_record)
        
    # SeqIO.write(records, f"{translated_seq_path}/gisaid.fasta", "fasta")







#     prot_data = prot_data.split(" ")
#     prot_data = [x for x in prot_data if x]

#     for i in prot_data[1:]:
#         i = float(i)
#         print(i)




    # # record = AlignIO.read(msa_out_file, "fasta")

    # # trimmed_seq1 = record[0]

    # # for i in range(len(record)):
    # #     seq2 = record[i]
    # #     percent_id_calc(trimmed_seq1, seq2, i)

    # msa(msa_in_file, msa_out_file)
    # # # print(MatrixInfo.blosum62)
    # #aligned_seqs(dna_path, uniprot_path)

    # #aligner_encodes(dna_path, uniprot_path)


    # # seqs = seqs.split("\n")


    # # text_file = "/home/nadzhou/Desktop/globin.txt"
    # # parse_percent_id_file(text_file)

    # # file = tables.openFile('/home/nadzhou/Desktop/dump/globin.mat')

    # # lon = file.root.lon[:]
    # # lat = file.root.lat[:]
    # # # Alternate syntax if the variable name is in a string
    # # varname = 'lon'
    # # lon = file.getNode('/' + varname)[:]

    # # print(file)



        # glob_file_path = Path("/home/nadzhou/Desktop/globin_dump")
    
    
    # results = parse_percent_id_file(glob_file_path)

    # origin_pre_msa_path = Path("/home/nadzhou/DEVELOPMENT/cov_seqs/msa_seqs")

    # align_high_score_seqs(origin_pre_msa_path, results)
    # # # record = extract_dna(dna_path, uniprot_path)
