from cov import trim_sequences


def test_trim_sequences():
    test_file = '../notebooks/gisaid_results/needle.fasta'
    test_result = '../notebooks/gisaid_results/trimmed_seqs.fasta'
    with open(test_result, 'r') as file:
        expected = file.read()
    results = trim_sequences(test_file)
    assert results == expected


