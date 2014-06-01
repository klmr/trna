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
    transcripts = read_sequences(transcript_filename)
    genes = read_genes(expression_stream, transcripts)
    return transcripts, genes


def read_genes(csvfile, transcripts):
    """
    Return a list of (gene name, expression) tuples read from a CSV file.
    """
    reader = csv.reader(csvfile, delimiter = '\t')
    return normalize([(row[0], float(row[1])) for row in reader], transcripts)


def normalize(genes, transcripts):
    """
    Calculate the FPKM for unnormalized gene counts

    The counts we get are already library size normalized so we ignore the
    total number of mapped reads.
    """
    def fpk(gene):
        transcript = transcripts[gene[0]]
        if transcript is None or transcript == '':
            return 0
        return 1000 * gene[1] / len(transcript)

    if len(filter(lambda g: g[1] != 1.0, genes)) == 0:
        return genes
    return [(gene[0], fpk(gene)) for gene in genes]


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


def codon_usage_for(seq, strand = 1, frame = 0):
    """
    Return a count list of codon usage for a given gene.
    """
    hist = empty_codon_dict()

    for codon in codons_for_frame(seq, strand, frame):
        if 'N' not in codon:
            hist[codon] += 1

    return hist


def codons_for_frame(seq, strand = 1, frame = 0):
    assert strand == 1 or strand == -1, 'strand must be -1 or 1'
    if strand == -1:
        seq = seq[::-1]

    upto = (len(seq) - frame) // 3
    for pos in range(upto):
        realpos = pos * 3 + frame
        yield seq[realpos : realpos + 3]


def weighted_codon_usage_for(gene, seq, expression):
    """
    Return a count list of codon usage for a given gene, weighted by the
    expression base mean of the gene.
    """
    if expression == 0:
        return empty_codon_dict()

    hist = codon_usage_for(seq)
    return dict((codon, freq * expression)
                for codon, freq in hist.items())


def print_results(totals):
    for codon in sorted(genetic_code.codons):
        print '%s\t%0.3f' % (codon, totals[codon])


if __name__ == '__main__':
    sys.exit(main())
