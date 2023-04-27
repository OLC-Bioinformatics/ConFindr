#!/usr/bin/env python
from confindr_src.methods import check_acceptable_xmx, check_for_databases_and_download, check_valid_base_fraction, \
    dependency_check, find_paired_reads, find_unpaired_reads, find_contamination, get_version, write_output
import multiprocessing
import subprocess
import traceback
import argparse
import logging
import shutil
import os


def confindr(args):
    # Check for dependencies.
    all_dependencies_present = True
    # Re-enable minimap2 as dependency once nanopore stuff actually works.
    if args.data_type == 'Illumina':
        dependencies = ['bbmap.sh', 'bbduk.sh', 'mash', 'kma']
    else:
        dependencies = ['bbduk.sh', 'mash', 'minimap2', 'kma']

    for dependency in dependencies:
        if dependency_check(dependency) is False:
            logging.error('Dependency {} not found. Please make sure it is installed and present'
                          ' on your $PATH.'.format(dependency))
            all_dependencies_present = False
    if not all_dependencies_present:
        logging.error('Could not find all necessary dependencies, quitting...')
        quit(code=1)

    # Check that the base fraction specified actually makes sense.
    if check_valid_base_fraction(args.base_fraction_cutoff) is False:
        logging.error('Base fraction must be between 0 and 1 if specified. Input value was: {}'
                      .format(args.base_fraction_cutoff))
        quit(code=1)

    # If user specified Xmx, make sure that they actually entered a value that will work. If not, the method will tell
    # them what they did wrong. Then quit.
    if args.Xmx:
        valid_xmx = check_acceptable_xmx(args.Xmx)
        if valid_xmx is False:
            quit(code=1)

    # Don't yet have cgmlst support with Nanopore reads - don't let user do this.
    if args.cgmlst and args.data_type == 'Nanopore':
        logging.error('ERROR: cgMLST schemes not yet supported for Nanopore reads. Quitting...')
        quit(code=1)

    if args.data_type == 'Nanopore':
        logging.warning('WARNING: Nanopore contamination detection is highly experimental. Any results should be taken '
                        'with several very large grains of salt. If you are going to try this, try setting -q to '
                        'somewhere in the range of 12-15, and only count things as contaminated that have at least '
                        '10 contaminating SNVs. Even then, results may be wonky. In particular, samples with lots of '
                        'depth will probably always show up as contaminated.')

    # Make the output directory.
    if not os.path.isdir(args.output_name):
        os.makedirs(args.output_name)
    # Remove any reports created by previous iterations of ConFindr
    try:
        os.remove(os.path.join(args.output_name, 'confindr_report.csv'))
    except FileNotFoundError:
        pass
    # Set the minimum number of matching hashes
    min_matching_hashes = args.min_matching_hashes
    # Check if databases necessary to run are present, and download them if they aren't
    check_for_databases_and_download(database_location=args.databases)

    # Figure out what pairs of reads, as well as unpaired reads, are present.
    paired_reads = find_paired_reads(args.input_directory,
                                     forward_id=args.forward_id,
                                     reverse_id=args.reverse_id)
    unpaired_reads = find_unpaired_reads(args.input_directory,
                                         forward_id=args.forward_id,
                                         reverse_id=args.reverse_id,
                                         find_fasta=args.fasta)
    # Consolidate read lists
    reads = sorted(paired_reads + unpaired_reads)
    # Process paired reads, one sample at a time.
    for fastq in reads:
        if len(fastq) == 1:
            sample_name = os.path.split(fastq[0])[-1].split('.')[0]
        else:
            sample_name = os.path.split(fastq[0])[-1].split(args.forward_id)[0]
        logging.info('Beginning analysis of sample {}...'.format(sample_name))
        try:
            find_contamination(pair=fastq,
                               forward_id=args.forward_id,
                               threads=args.threads,
                               output_folder=args.output_name,
                               databases_folder=args.databases,
                               keep_files=args.keep_files,
                               quality_cutoff=args.quality_cutoff,
                               base_cutoff=args.base_cutoff,
                               base_fraction_cutoff=args.base_fraction_cutoff,
                               cgmlst_db=args.cgmlst,
                               xmx=args.Xmx,
                               tmpdir=args.tmp,
                               data_type=args.data_type,
                               use_rmlst=args.rmlst,
                               min_matching_hashes=min_matching_hashes,
                               fasta=args.fasta,
                               debug=args.verbosity)
        except subprocess.CalledProcessError:
            # If something unforeseen goes wrong, traceback will be printed to screen.
            # We then add the sample to the report with a note that it failed.
            multi_positions = 0
            genus = 'Error processing sample'
            write_output(output_report=os.path.join(args.output_name, 'confindr_report.csv'),
                         sample_name=sample_name,
                         multi_positions=multi_positions,
                         genus=genus,
                         total_gene_length=0,
                         database_download_date='ND')
            logging.warning('Encountered error when attempting to run ConFindr on sample '
                            '{sample}. Skipping...'.format(sample=sample_name))
            logging.warning('Error encountered was:\n{}'.format(traceback.format_exc()))
            if args.keep_files is False:
                shutil.rmtree(os.path.join(args.output_name, sample_name))
    if args.keep_files is False and args.tmp is not None:
        shutil.rmtree(args.tmp)
    logging.info('Contamination detection complete!')


def main():
    version = get_version()
    cpu_count = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_directory',
                        type=str,
                        required=True,
                        help='Folder that contains fastq files you want to check for contamination. '
                             'Will find any file that contains .fq or .fastq in the filename.')
    parser.add_argument('-o', '--output_name',
                        type=str,
                        required=True,
                        help='Base name for output/temporary directories.')
    parser.add_argument('-d', '--databases',
                        type=str,
                        default=os.environ.get('CONFINDR_DB', os.path.expanduser('~/.confindr_db')),
                        help='Databases folder. To download these, you will need to get access to the rMLST databases. '
                             'For complete instructions on how to do this, please see '
                             'https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases')
    parser.add_argument('--rmlst',
                        default=False,
                        action='store_true',
                        help='Activate to prefer using rMLST databases over core-gene derived databases. By default,'
                             'ConFindr will use core-gene derived databases where available.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=cpu_count,
                        help='Number of threads to run analysis with.')
    parser.add_argument('-tmp', '--tmp',
                        type=str,
                        help='If your ConFindr databases are in a location you don\'t have write access to, '
                             'you can enter this option to specify a temporary directory to put genus-specific '
                             'databases to.')
    parser.add_argument('-k', '--keep_files',
                        default=False,
                        action='store_true',
                        help='By default, intermediate files are deleted. Activate this flag to keep intermediate '
                             'files.')
    parser.add_argument('-q', '--quality_cutoff',
                        type=int,
                        default=20,
                        help='Base quality needed to support a multiple allele call. Defaults to 20.')
    parser.add_argument('-b', '--base_cutoff',
                        type=int,
                        default=3,
                        help='Number of bases necessary to support a multiple allele call, and automatically '
                        'increments based upon gene-specific quality score, length and depth of coverage. '
                        'Default is 3.')
    parser.add_argument('-bf', '--base_fraction_cutoff',
                        type=float,
                        default=0.05,
                        help='Fraction of bases necessary to support a multiple allele call. Particularly useful when '
                             'dealing with very high coverage samples. Default is 0.05.')
    parser.add_argument('-e', '--error_cutoff',
                        type=float,
                        default=1.0,
                        help='Value to use for the calculated error cutoff when setting the base cutoff value. '
                             'Default is 1.0%%.')
    parser.add_argument('-fid', '--forward_id',
                        type=str,
                        default='_R1',
                        help='Identifier for forward reads.')
    parser.add_argument('-rid', '--reverse_id',
                        type=str,
                        default='_R2',
                        help='Identifier for reverse reads.')
    parser.add_argument('-v', '--version',
                        action='version',
                        version=version)
    parser.add_argument('-dt', '--data_type',
                        choices=['Illumina', 'Nanopore'],
                        default='Illumina',
                        help='Type of input data. Default is Illumina, but can be used for Nanopore too. No PacBio '
                             'support (yet).')
    parser.add_argument('-Xmx', '--Xmx',
                        type=str,
                        help='Very occasionally, parts of the pipeline that use the BBMap suite will have their memory '
                             'reservation fail and request not enough, or sometimes negative, memory. If this happens '
                             'to you, you can use this flag to override automatic memory reservation and use an amount '
                             'of memory requested by you. -Xmx 20g will specify 20 gigs of RAM, and -Xmx 800m '
                             'will specify 800 megs.')
    parser.add_argument('-cgmlst', '--cgmlst',
                        type=str,
                        help='Path to a cgMLST database to use for contamination detection instead of using the default'
                             ' rMLST database. Sequences in this file should have headers in format '
                             '>genename_allelenumber. To speed up ConFindr runs, clustering the cgMLST database with '
                             'CD-HIT before running ConFindr is recommended. This is highly experimental, results '
                             'should be interpreted with great care.')
    parser.add_argument('--fasta',
                        default=False,
                        action='store_true',
                        help='If activated, will look for FASTA files instead of FASTQ for unpaired reads.')
    parser.add_argument('-verbosity', '--verbosity',
                        choices=['debug', 'info', 'warning'],
                        default='info',
                        help='Amount of output you want printed to the screen. Defaults to info, which should be good '
                             'for most users.')
    parser.add_argument('-m', '--min_matching_hashes',
                        default=150,
                        type=int,
                        help='Minimum number of matching hashes in a MASH screen in order for a genus to be considered '
                             'present in a sample. Default is 150')
    args = parser.parse_args()
    # Setup the logger. TODO: Different colors for different levels.
    if args.verbosity == 'info':
        logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                            level=logging.INFO,
                            datefmt='%Y-%m-%d %H:%M:%S')
    elif args.verbosity == 'debug':
        logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                            level=logging.DEBUG,
                            datefmt='%Y-%m-%d %H:%M:%S')
    elif args.verbosity == 'warning':
        logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                            level=logging.WARNING,
                            datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('Welcome to {version}! Beginning analysis of your samples...'.format(version=version))
    confindr(args)


if __name__ == '__main__':
    main()
