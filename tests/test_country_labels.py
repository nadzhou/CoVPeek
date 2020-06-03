import os

import scripts.country_labels as labels


def test_country_labels():
    data = labels.load_data('test.fasta')
    country_dict = labels.seq_records_for_country(data)
    assert 'England' in country_dict
    assert len(country_dict.keys()) == 1
    assert len(country_dict['England']) == 10
    var_counts = labels.variation_counts_at_position(country_dict, 3)
    assert len(var_counts['England'].items()) == 10
