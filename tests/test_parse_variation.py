import scripts.parse_variation as pv


def test_get_conservation_scores():
    parser = pv.DivergenceParser.retrieve_sequence('test.fasta')
    assert parser.npseqs.shape == (10, 60)
    scores = parser.conservation_scores()
    assert len(scores) == 60
    variable = parser.aminoacids_in_variable_positions()
    assert 3 in variable
    assert len(variable[3]) == 10
