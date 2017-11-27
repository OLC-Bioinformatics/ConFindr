#!/usr/bin/env python

import multiprocessing
import subprocess
import argparse
import shutil
import glob
import time
import csv
import os
import pysam
from Bio import SeqIO
from biotools import bbtools
from biotools import jellyfish
from accessoryFunctions.accessoryFunctions import printtime


def dependency_check(dependency):
    if shutil.which(dependency) is not None:
        return True
    else:
        return False


def find_paired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
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


def setup_genusspecific_database(database_folder, genus, genes_to_exclude):
    with open(os.path.join(database_folder, '{}_db.fasta'.format(genus)), 'w') as f:
        sequences = SeqIO.parse(os.path.join(database_folder, 'rMLST_combined.fasta'), 'fasta')
        for item in sequences:
            if item.id.split('_')[0] not in genes_to_exclude:
                f.write('>' + item.id + '\n')
                f.write(str(item.seq) + '\n')


def extract_rmlst_genes(pair, database, forward_out, reverse_out):
    bbtools.bbduk_bait(database, pair[0], forward_out, reverse_in=pair[1], reverse_out=reverse_out)


def subsample_reads(forward_in, reverse_in, coverage_level, genome_size, forward_out, reverse_out):
    bases_target = coverage_level * genome_size
    bbtools.subsample_reads(forward_in, forward_out, bases_target, forward_out, reverse_in=reverse_in, reverse_out=reverse_out)


def generate_kmers(forward_reads, reverse_reads, counts_file, kmer_size, tmpdir):
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    jellyfish.count(forward_reads, reverse_in=reverse_reads, count_file=os.path.join(tmpdir, 'mer_counts.jf'), kmer_size=kmer_size)
    jellyfish.dump(os.path.join(tmpdir, 'mer_counts.jf'), counts_file)
    shutil.rmtree(tmpdir)


def rename_kmers(input_kmers, output_kmers, cutoff):
    with open(input_kmers) as f:
        fastas = f.readlines()

    num_mers = 0
    sequences = list()
    for i in range(len(fastas)):
        if '>' in fastas[i]:
            if int(fastas[i].replace('>', '')) >= cutoff:
                num_mers += 1
                sequences.append(fastas[i].rstrip() + '_' + str(num_mers) + '\n' + fastas[i + 1])
    # Write out our solid kmers to file to be used later.
    with open(output_kmers, 'w') as f:
        f.write(''.join(sequences))

    return num_mers


def parse_bamfile(bamfile, kmer_size):
    mismatch_kmer_headers = list()
    bam_handle = pysam.AlignmentFile(bamfile, 'rb')
    for match in bam_handle:
        if match.cigarstring is not None:
            if '1X' in match.cigarstring and match.query_alignment_length == kmer_size:
                mismatch_kmer_headers.append(bam_handle.getrname(match.reference_id))
    return mismatch_kmer_headers


def find_contamination(pair, args, genus):
    sample_start = time.time()
    # Main method for finding contamination - works on one pair at a time.
    sample_name = os.path.split(pair[0])[-1].split(args.forward_id)[0]
    printtime('Finding contamination for {}...'.format(sample_name), sample_start)
    # Need to:
    # Setup genus-specific databases, if necessary.
    if genus != 'NA':
        sample_database = os.path.join(args.databases, '{}_db.fasta'.format(genus))
        if not os.path.isfile(os.path.join(args.databases, '{}_db.fasta'.format(genus))):
            printtime('Setting up genus-specific database for genus {}...'.format(genus), sample_start)
            genes_to_excude = find_genusspecific_alleles(os.path.join(args.databases, 'profiles.txt'), genus)
            setup_genusspecific_database(args.databases, genus, genes_to_excude)
    else:
        sample_database = os.path.join(args.databases, 'rMLST_combined.fasta')
    # Extract rMLST reads and quality trim.
    sample_tmp_dir = os.path.join(args.output_name, sample_name)
    if not os.path.isdir(sample_tmp_dir):
        os.makedirs(sample_tmp_dir)
    printtime('Extracting rMLST genes...', sample_start)
    extract_rmlst_genes(pair, sample_database,
                        forward_out=os.path.join(sample_tmp_dir, 'rmlst_R1.fastq.gz'),
                        reverse_out=os.path.join(sample_tmp_dir, 'rmlst_R2.fastq.gz'))
    printtime('Quality trimming...', sample_start)
    bbtools.bbduk_trim(forward_in=os.path.join(sample_tmp_dir, 'rmlst_R1.fastq.gz'),
                       reverse_in=os.path.join(sample_tmp_dir, 'rmlst_R2.fastq.gz'),
                       forward_out=os.path.join(sample_tmp_dir, 'trimmed_R1.fastq.gz'),
                       reverse_out=os.path.join(sample_tmp_dir, 'trimmed_R2.fastq.gz'))
    # Now do the actual contamination detection cycle the number of times specified by arguments.
    printtime('Beginning {} cycles of contamination detection'.format(str(args.number_subsamples)), sample_start)
    for i in range(args.number_subsamples):
        # Subsample
        subsample_reads(forward_in=os.path.join(sample_tmp_dir, 'trimmed_R1.fastq.gz'),
                        reverse_in=os.path.join(sample_tmp_dir, 'trimmed_R2.fastq.gz'),
                        coverage_level=args.subsample_depth,
                        genome_size=35000,  # This is the sum of the longest allele for each rMLST gene.
                        forward_out=os.path.join(sample_tmp_dir, 'subsample_{}_R1.fastq.gz'.format(str(i))),
                        reverse_out=os.path.join(sample_tmp_dir, 'subsample_{}_R2.fastq.gz'.format(str(i))))
        # Kmerize subsampled reads.
        generate_kmers(forward_reads=os.path.join(sample_tmp_dir, 'subsample_{}_R1.fastq.gz'.format(str(i))),
                       reverse_reads=os.path.join(sample_tmp_dir, 'subsample_{}_R2.fastq.gz'.format(str(i))),
                       counts_file=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                       kmer_size=args.kmer_size,
                       tmpdir=os.path.join(sample_tmp_dir, 'tmp'))
        # Rename kmers so each has a unique ID, and count the number of kmers.
        num_kmers = rename_kmers(input_kmers=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                                 output_kmers=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                                 cutoff=args.kmer_cutoff)
        # Find mismatches.
        # Step 1 of mismatch finding: Run bbmap with the kmer file.
        bbtools.bbmap(reference=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                      forward_in=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                      ambig='all',
                      overwrite='true',
                      out_bam=os.path.join(sample_tmp_dir, 'subsample_{}.bam'.format(str(i))))
        # Step 2 of mismatch finding: Parse the bamfile created by bbmap to find one mismatch kmers.
        fasta_ids = parse_bamfile(os.path.join(sample_tmp_dir, 'subsample_{}.bam'.format(str(i))), args.kmer_size)
        # Step 2.5: Create a dictionary so you know which ID goes with which sequence.
        mer_dict = dict()
        with open(os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i)))) as f:
            mers = f.readlines()
        for j in range(0, len(mers), 2):
            key = mers[j].replace('>', '')
            key = key.replace('\n', '')
            mer_dict[key] = mers[j + 1]
        # Step 3 of mismatch finding: Blast kmers found from parsing against database to make sure they're
        # not overhangs into non-RMLST regions that could cause false positives.
        # Find cross contamination.
    # Create contamination report.


if __name__ == '__main__':
    start = time.time()
    cpu_count = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_directory',
                        type=str,
                        required=True,
                        help="Folder that contains fastq files you want to check for contamination. "
                             "Will find any fastq file that contains .fq or .fastq in the filename.")
    parser.add_argument('-o', '--output_name',
                        type=str,
                        required=True,
                        help='Base name for output/temporary directories.')
    parser.add_argument('-d', '--databases',
                        type=str,
                        required=True,
                        help='Databases folder. Should contain rMLST_combined.fasta, profiles.txt, '
                             'and refseq.msh as well as RefSeqSketchesDefaults.msh')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=cpu_count,
                        help='Number of threads to run analysis with.')
    parser.add_argument('-n', '--number_subsamples',
                        type=int,
                        default=5,
                        help='Number of times to subsample.')
    parser.add_argument('-k', '--kmer-size',
                        type=int,
                        default=31,
                        help='Kmer size to use for contamination detection.')
    parser.add_argument('-s', '--subsample_depth',
                        type=int,
                        default=20,
                        help='Depth to subsample to. Higher increases sensitivity, but also false positive '
                             'rate. Default is 20.')
    parser.add_argument('-c', '--kmer_cutoff',
                        type=int,
                        default=2,
                        help='Number of times you need to see a kmer before it is considered trustworthy.'
                             ' Kmers with counts below this number will be discarded.')
    parser.add_argument('-fid', '--forward_id',
                        type=str,
                        default='_R1',
                        help='Identifier for forward reads.')
    parser.add_argument('-rid', '--reverse_id',
                        type=str,
                        default='_R2',
                        help='Identifier for reverse reads.')
    # Check for dependencies.
    dependencies = ['jellyfish', 'bbmap.sh', 'bbduk.sh']
    for dependency in dependencies:
        if dependency_check(dependency) is False:
            print('WARNING: Dependency {} not found. ConFindr will likely crash!'.format(dependency))
    # TODO: Add a dependency check here.
    args = parser.parse_args()
    # Make the output directory.
    if not os.path.isdir(args.output_name):
        os.makedirs(args.output_name)
    # Run mashsippr on all sample, and then read the file individually for each sample.
    # TODO: Get the mashsippr created folders cleaned up, and redirect the stdout of mashsippr to not the terminal
    run_mashsippr(args.input_directory, args.output_name, args.databases)
    paired_reads = find_paired_reads(args.input_directory, forward_id=args.forward_id, reverse_id=args.reverse_id)
    for pair in paired_reads:
        sample_name = os.path.split(pair[0])[-1].split(args.forward_id)[0]
        genus = read_mashsippr_output(os.path.join(args.output_name, 'reports/mash.csv'), sample_name)
        find_contamination(pair, args, genus)

