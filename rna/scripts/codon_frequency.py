#!/usr/bin/env python

"""
Generate codon usage bias histograms from expression data.
"""

import sys
import itertools
from Bio import SeqIO
import csv
from collections import defaultdict
from operator import itemgetter
from goo.genetic_code import genetic_code, empty_codon_dict


def main():
    # Format: $0 transcripts.fa < gene_ids/expressions
    transcripts, genes = read_input(sys.stdin, sys.argv[1])
    totals = histogram(transcripts, genes)
    print_results(totals)


def read_input(expression_stream, transcript_filename):
    """
    Read the input files and output statistics.

    The input files are parsed into appropriate structures and returned
    as a tuple.
    """
    genes = read_genes(expression_stream)
    transcripts = read_sequences(transcript_filename)
    return transcripts, genes


def read_genes(csvfile):
    """
    Return a list of (gene name, expression) tuples read from a CSV file.
    """
    reader = csv.reader(csvfile, delimiter = '\t')
    return [(row[0], float(row[1])) for row in reader]


def read_sequences(filename):
    """
    Read sequences from a FASTA file containing protein-coding gene regions.

    Return a dictionary of gene IDs and their sequences.
    """
    # Format:
    #
    #   >ENSGENE_ID|ENSTRANSCRIPT_ID
    #   SEQUENCE
    #
    records = [(record.id.split('|')[0], str(record.seq))
            for record in SeqIO.parse(filename, 'fasta')]
    valid = filter(lambda r: is_valid_transcript(r[1]), records)
    by_gene = group_by(valid, key = itemgetter(0), value = itemgetter(1))
    # For each gene, select the longest transcript as "representative".
    best_for_gene = map(lambda r: (r[0], max(r[1], key = len)), by_gene)
    return defaultdict(generate_const(None), best_for_gene)


def is_valid_transcript(seq):
    return (len(seq) > 2 and
            len(seq) % 3 == 0 and
            seq[: 3] in genetic_code.start and
            seq[-3 : ] in genetic_code.stop)


def group_by(iterable, key, value = lambda x: x):
    """
    Groups elements in an iterable by some key criterion.

    The `key` keyword argument needs to supply a function which maps entries to
    a key. The return value of this method is an iterator of tuples where the
    first item is the key and the second item is a list of values.
    An optionally provided ``value`` argument provides a list for mapping items
    to values (default: identity).
    """
    ret = defaultdict(list)
    for item in iterable:
        ret[key(item)].append(value(item))

    return ret.iteritems()


def generate_const(const):
    return lambda: const


def histogram(transcripts, genes):
    """
    Calculate the weighted frequencies of all occurring in the transcripts.

    The transcripts are weighed by their transcription level.
    Returns a dictionary of codons and their frequency.
    """
    def add_to_total(total, hist):
        for key, value in hist.items():
            total[key] += value

#    if len(genes) == 0:
#        # Calculate background distribution
#        genes = [(gene, 1.0) for gene in transcripts.keys()]
#        sys.stderr.write('{}\n'.format(len(genes)))

    total = empty_codon_dict()
    failures = 0
    for expression in genes:
        gene = expression[0]
        seq = transcripts[gene]
        if seq is None or seq == '':
            failures += 1
            continue

        hist = weighted_codon_usage_for(gene, seq, expression[1])
        add_to_total(total, hist)

    return total


def codon_usage_for(gene, seq):
    """
    Return a count list of codon usage for a given gene.
    """
    hist = empty_codon_dict()

    for pos in range(len(seq) // 3):
        codon = seq[pos * 3 : pos * 3 + 3]
        if 'N' not in codon:
            hist[codon] += 1

    return hist


def weighted_codon_usage_for(gene, seq, expression):
    """
    Return a count list of codon usage for a given gene, weighted by the
    expression base mean of the gene.
    """
    if expression == 0:
        return empty_codon_dict()

    # TODO Normalise by transcript length -- or not?! -- Nah. For now.
    hist = codon_usage_for(gene, seq)
    return dict((codon, freq * expression)
                for codon, freq in hist.items())


def print_results(totals):
    for codon in sorted(genetic_code.codons):
        print '%s\t%0.3f' % (codon, totals[codon])


if __name__ == '__main__':
    sys.exit(main())
