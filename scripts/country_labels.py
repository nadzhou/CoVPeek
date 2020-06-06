from typing import List, Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter
from collections import defaultdict
import re


def load_data(path: str):
    # loads records from fasta
    seqs: List = list(SeqIO.parse(path, "fasta"))
    return seqs


def seq_records_for_country(seqs: List[SeqRecord]) -> Dict[str, List[SeqRecord]]:
    # sorts the sequence records based on their IDs
    records_by_country = defaultdict(list)
    skipped = []
    for seq in seqs:
        try:
            country = country_from_id(seq.id)
        except ValueError:
            skipped.append(seq.id)
            continue
        records_by_country[country].append(seq)
    print(f"{len(skipped)} records skipped")
    return records_by_country


def country_from_id(record_id: str) -> str:
    # seaches for a word, number-number-word id format
    regex_pattern = r"([a-zA-Z]+)"
    match = re.search(pattern=regex_pattern, string=record_id)
    if not match:
        raise ValueError(f"Invalid id format: {record_id}")
    countryname = match.group(1)
    return countryname


def variation_counts_at_position(records_by_country: Dict[str, List[SeqRecord]], position: int) -> Dict[str, Counter]:
    counts_by_country = {}
    for country, records in records_by_country.items():
        aa_variants = [record.seq[position-1] for record in records]
        counts = Counter(aa_variants)
        counts_by_country[country] = counts
    return counts_by_country


def main():
    # we store aligned data in non-shared folder
    seq_path = "../operations/gisaid_results/aligned.fasta"
    records = load_data(seq_path)
    country_dict = seq_records_for_country(records)

    print('Position 705')
    print("Country, Variation")
    var_counts = variation_counts_at_position(country_dict, 705)
    for country, counts in var_counts.items():
        print(f"{country}, {dict(counts)}")

    # print('Position 614')
    # var_counts = variation_counts_at_position(country_dict, 614)
    # for country, counts in var_counts.items():
    #     print(f"country: {country}, \t variations: {dict(counts)}")

    # print('Position 615')
    # var_counts = variation_counts_at_position(country_dict, 615)
    # for country, counts in var_counts.items():
    #     print(f"country: {country}, \t variations: {dict(counts)}")


if __name__ == '__main__':
    main()
