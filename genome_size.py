#!/usr/bin/env python
import os


def get_peak_kmers(histo_file):
    most_kmers = 0
    total_kmers = 0
    peak_depth = 0
    with open(histo_file) as hist_data:
        for row in hist_data:
            data = row.split()
            kmer_depth = int(data[0])
            num_kmers = int(data[1])
            total_kmers += num_kmers * kmer_depth
            if kmer_depth != 1 and num_kmers > most_kmers:
                most_kmers = num_kmers
                peak_depth = kmer_depth

    return peak_depth, total_kmers


def get_genome_size(total_kmers, peak_depth):
    return total_kmers/peak_depth


def run_jellyfish_histo():
    cmd = 'jellyfish histo mer_counts.jf > histogram.txt'
    os.system(cmd)

