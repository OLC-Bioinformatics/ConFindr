#!/usr/bin/env python

import subprocess
import shutil
import glob
import csv
import os
from Bio import SeqIO
from biotools import bbtools


def dependency_check(dependency):
    if shutil.which(dependency) is not None:
        return True
    else:
        return False


def find_paired_reads(fastq_directory, forward_id='R1', reverse_id='R2'):
    """
    Looks at a directory to try to find paired fastq files. Should be able to find anything fastq.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default R1.
    :param reverse_id: Identifier for reverse reads. Default R2.
    :return: List containing pairs of fastq files, in format [[forward_1, reverse_1], [forward_2, reverse_2]], etc.
    """
    pair_list = list()
    fastq_files = glob.glob(fastq_directory + '/*.f*q*')
    for name in fastq_files:
        if forward_id in name and os.path.isfile(name.replace(forward_id, reverse_id)):
            pair_list.append([name, name.replace(forward_id, reverse_id)])
    return pair_list


def run_mashsippr(sequence_dir, output_dir, database_dir):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    cmd = 'python -m confindr.mashsippr -s {sequence_dir} -t {database_dir} {output_dir}'.format(sequence_dir=sequence_dir,
                                                                                                 database_dir=database_dir,
                                                                                                 output_dir=output_dir)
    subprocess.call(cmd, shell=True)
    if os.path.isfile(os.path.join(output_dir, 'reports/mash.csv')):
        return True
    else:
        return False


def read_mashsippr_output(mashsippr_result_file, sample):
    genus = 'NA'
    with open(mashsippr_result_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Strain'] == sample:
                genus = row['ReferenceGenus']
    return genus


def find_genusspecific_alleles(profiles_file, target_genus):
    genes_to_exclude = list()
    with open(profiles_file) as f:
        lines = f.readlines()
    for line in lines:
        line = line.rstrip()
        genus = line.split(':')[0]
        if genus == target_genus:
            genes = line.split(':')[1]
            genes_to_exclude = genes.split(',')
    return genes_to_exclude


def extract_rmlst_genes(pair, database, forward_out, reverse_out):
    bbtools.bbduk_bait(database, pair[0], forward_out, reverse_in=pair[1], reverse_out=reverse_out)


