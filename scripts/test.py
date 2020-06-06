from Bio import SeqIO



record = list(SeqIO.parse("gisaid19.fasta", "fasta"))

print(len(record))