#!/usr/bin/env python

"""
Generate the background codon usage with shifted reading frame.
"""

import codon_frequency
import sys

def main():
    # Format: $0 transcripts.fa strand frame
    transcripts = list(codon_frequency.read_sequences(sys.argv[1]).values())
    strand = int(sys.argv[2])
    frame = int(sys.argv[3])
    totals = histogram(transcripts, strand, frame)
    codon_frequency.print_results(totals)


def histogram(transcripts, strand, frame):
    """
    Calculate the background frequencies of all occurring in the transcripts.

    Returns a dictionary of codons and their frequency.
    """
    def add_to_total(hist):
        for key, value in hist.items():
            total[key] += value

    total = codon_frequency.empty_codon_dict()
    failures = 0

    for transcript in transcripts:
        if transcript is None or transcript == '':
            failures += 1
            continue

        add_to_total(codon_frequency.codon_usage_for(transcript, strand, frame))

    return total


if __name__ == '__main__':
    sys.exit(main())
