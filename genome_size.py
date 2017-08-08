#!/usr/bin/env python
import os


def get_peak_kmers(histo_file):
    """
    :param histo_file: Histogram file created by jellyfish histo
    :return: peak_depth: Coverage depth that had the most kmers (int)
             total_kmers: Total number of kmers (int)
    """
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
    # This is a simplistic way of getting peak depth, but it works well for my purposes.
    return peak_depth, total_kmers


def get_genome_size(total_kmers, peak_depth):
    return total_kmers/peak_depth


def run_jellyfish_histo():
    cmd = 'jellyfish histo mer_counts.jf > histogram.txt'
    os.system(cmd)

