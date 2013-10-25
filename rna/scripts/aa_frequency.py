#!/usr/bin/env python

import sys
from goo.genetic_code import genetic_code, empty_codon_dict
import codon_frequency as cf
from collections import defaultdict


def main():
    # Format: $0 transcripts.fa < gene_ids/expressions
    transcripts, genes = cf.read_input(sys.stdin, sys.argv[1])
    totals = histogram(transcripts, genes)
    print_results(totals)


def histogram(transcripts, genes):
    """
    Equivalent to codon_frequency.histogram, just for amino acids.
    """
    return freq_codon_to_aas(cf.histogram(transcripts, genes))


def freq_codon_to_aas(codon_total):
    aa_total = defaultdict(int)
    for codon, weight in codon_total.items():
        aa_total[genetic_code.translate(codon)] += weight

    return aa_total


def print_results(totals):
    for aa in sorted(genetic_code.amino_acids):
        print '%s\t%0.3f' % (aa, totals[aa])


if __name__ == '__main__':
    sys.exit(main())
