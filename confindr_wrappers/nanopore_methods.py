#!/usr/bin/env python

import gzip
from Bio import SeqIO


def nanopore_or_illumina(fastq_file, read_length_cutoff=301):
    """
    Makes a somewhat educated guess about whether something is illumina or nanopore data
    :param fastq_file: Path to a fastq file, either gzipped or uncompressed.
    :param read_length_cutoff: Read length cutoff - if average read length is greater than this, assume reads are
    nanopore, otherwise assume illumina.
    :return: 'Nanopore' if file is determined to be a nanopore file, 'Illumina' if it's illumina.
    """
    # Parse through the first 10 records of the fastq file to find average length and quality.
    if fastq_file.endswith('.gz'):
        with gzip.open(fastq_file, 'rt') as fastq_handle:
            average_read_length = average_length_and_quality(fastq_handle,
                                                             records_to_parse=10)
    else:
        with open(fastq_file) as fastq_handle:
            average_read_length = average_length_and_quality(fastq_handle,
                                                             records_to_parse=10)

    # Anything with over 300 average read length is probably nanopore (or PacBio I suppose, but that
    # should be able to be treated same as nanopore).
    if average_read_length > read_length_cutoff:
        return 'Nanopore'
    else:
        return 'Illumina'


def average_length_and_quality(fastq_handle, records_to_parse=10):
    records_parsed = 0
    read_lengths = list()
    for record in SeqIO.parse(fastq_handle, 'fastq'):
        records_parsed += 1
        read_lengths.append(len(record))
        if records_parsed >= records_to_parse:
            break
    average_read_length = sum(read_lengths)/len(read_lengths)
    return average_read_length

