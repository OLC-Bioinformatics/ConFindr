#!/usr/bin/env python3
"""
confindr.py

Command-line entry point for ConFindr contamination detection.

This module implements the command-line interface used to run the
ConFindr contamination detection pipeline. It parses user-supplied
arguments, validates dependencies and inputs, and dispatches per-sample
work to the pipeline functions in :mod:`confindr_src.methods`.

Features:
- Supports probabilistic scoring via ``--use-prob-scoring``
- Configurable thresholds, memory settings, and alternate databases

See the project README for full usage and caveats.
"""

# Standard library imports
from glob import glob
import argparse
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import traceback

# Ensure the repository root is on sys.path when running the script directly
_SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_SCRIPT_DIR, '..'))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# Third party imports
try:
    import coloredlogs
except ImportError:  # pragma: no cover
    coloredlogs = None

# Local imports
from confindr_src.methods import (
    _valid_downsample_depth,
    check_acceptable_xmx,
    check_for_databases_and_download,
    check_valid_base_fraction,
    dependency_check,
    find_contamination,
    find_paired_reads,
    find_unpaired_reads,
    recommend_xmx,
    write_output,
)

from confindr_src import methods  # module import for runtime config
from confindr_src.version import __version__


def confindr(
    *,  # Enforce keyword arguments
    args: argparse.Namespace
) -> None:
    """
    Execute ConFindr analysis for parsed CLI arguments.

    This function performs argument validation, dependency checks, and
    iterates samples invoking the core pipeline for each sample.

    Args:
        args: Parsed argparse Namespace from the CLI. Expected attributes
            include at least:
            - input_directory: path containing read files
            - output_name: base output/temporary directory name
            - databases: path to ConFindr databases
            - threads: number of worker threads
            - quality_cutoff: base quality cutoff (int)
            - base_cutoff: integer base cutoff for SNV calls
            - base_fraction_cutoff: fractional cutoff for SNV calls
            - use_prob_scoring: enable probabilistic scoring (bool)
            - score_threshold: threshold for probabilistic scoring (float)

    Returns:
        None
    """
    # Check for dependencies.
    all_dependencies_present = True

    # Re-enable minimap2 as dependency once nanopore stuff actually works.
    if args.data_type == 'Illumina':
        dependencies = ['bbmap.sh', 'bbduk.sh', 'mash', 'kma']
    else:
        dependencies = ['bbduk.sh', 'mash', 'minimap2', 'kma']

    for dependency in dependencies:
        if dependency_check(
            dependency=dependency
        ) is False:
            logging.error(
                'Dependency %s not found. Please make sure it is '
                'installed and present on your $PATH.',
                dependency
            )
            all_dependencies_present = False
    if not all_dependencies_present:
        logging.error(
            'Could not find all necessary dependencies, quitting...'
        )
        sys.exit(1)

    # Check that the base fraction specified actually makes sense.
    if check_valid_base_fraction(
        base_fraction=args.base_fraction_cutoff
    ) is False:
        logging.error(
            'Base fraction must be between 0 and 1 if specified. Input value '
            'was: %s',
            args.base_fraction_cutoff
        )
        sys.exit(1)

    # If a user specified -Xmx, ensure they supplied a valid memory string.
    # The helper function will report errors and we should quit on failure.
    if args.Xmx:
        valid_xmx = check_acceptable_xmx(xmx_string=args.Xmx)
        if valid_xmx is False:
            sys.exit(1)

    # cgMLST schemes are not yet supported for Nanopore reads; prevent use.
    if args.cgmlst and args.data_type == 'Nanopore':
        logging.error(
            'ERROR: cgMLST schemes not yet supported for Nanopore reads. '
            'Quitting...'
        )
        sys.exit(1)

    # Warn users about Nanopore contamination detection caveats.
    if args.data_type == 'Nanopore':
        logging.warning(
            'WARNING: Nanopore contamination detection is highly '
            'experimental. Results should be interpreted with caution. '
            'If you try this, set -q to around 12-15 and only consider '
            'samples with at least 10 contaminating SNVs as contaminated. '
            'High-depth samples may appear contaminated and results can be '
            'unreliable.'
        )

    # Make the output directory.
    os.makedirs(args.output_name, exist_ok=True)

    # Remove any reports created by previous iterations of ConFindr
    try:
        os.remove(os.path.join(args.output_name, 'confindr_report.tsv'))
    except FileNotFoundError:
        pass

    # Set the minimum number of matching hashes
    min_matching_hashes = args.min_matching_hashes

    # Check whether the necessary databases are present, and download
    # them if they are missing.
    check_for_databases_and_download(database_location=args.databases)

    # Figure out what pairs of reads, as well as unpaired reads, are present.
    paired_reads = find_paired_reads(
        fastq_directory=args.input_directory,
        forward_id=args.forward_id,
        reverse_id=args.reverse_id
    )
    unpaired_reads = find_unpaired_reads(
        fastq_directory=args.input_directory,
        forward_id=args.forward_id,
        reverse_id=args.reverse_id,
        find_fasta=args.fasta
    )

    # Consolidate read lists
    reads = sorted(paired_reads + unpaired_reads)

    # Process paired reads, one sample at a time.
    for fastq in reads:
        if len(fastq) == 1:
            sample_name = os.path.split(fastq[0])[-1].split('.')[0]
        else:
            sample_name = os.path.split(fastq[0])[-1].split(args.forward_id)[0]
        logging.info('Beginning analysis of sample %s...', sample_name)

        # Run the contamination finding pipeline for this sample.
        try:
            find_contamination(
                pair=fastq,
                forward_id=args.forward_id,
                threads=args.threads,
                output_folder=args.output_name,
                databases_folder=args.databases,
                keep_files=args.keep_files,
                quality_cutoff=args.quality_cutoff,
                min_quality=args.min_quality,
                base_cutoff=args.base_cutoff,
                base_fraction_cutoff=args.base_fraction_cutoff,
                cgmlst_db=args.cgmlst,
                xmx=args.Xmx,
                tmpdir=args.tmp,
                data_type=args.data_type,
                use_rmlst=args.rmlst,
                min_matching_hashes=min_matching_hashes,
                fasta=args.fasta,
                use_prob_scoring=args.use_prob_scoring,
                score_threshold=args.score_threshold,
                max_expected_positions=args.max_expected_positions,
                downsample_depth=args.downsample_depth,
                subreplicates=args.subreplicates,
                subreplicate_seed=args.subreplicate_seed,
                subreplicate_consensus=args.subreplicate_consensus,
            )

            # Debug: scan per-sample position-level TSV files to check for
            # duplicated or excessive calls. This helps identify whether
            # multiple files contain the same coordinates.
            if args.verbosity == 'debug':
                sample_dir = os.path.join(args.output_name, sample_name)
                try:
                    pos_files = glob(
                        os.path.join(
                            sample_dir, '*positions*.tsv'
                        )
                    )
                    total_lines = 0
                    coords = set()
                    for pf in pos_files:
                        with open(pf, encoding='utf-8') as fh:
                            for line in fh:
                                if not line.strip():
                                    continue
                                if line.startswith(
                                    'contig'
                                ) or line.startswith('#'):
                                    # skip header lines
                                    continue
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    coords.add(f'{parts[0]}:{parts[1]}')
                                    total_lines += 1
                    # Check for duplicated coordinates
                    if len(coords) > 0:
                        sample_coords_list = list(sorted(coords))[:10]
                        logging.debug(
                            'Sample %s: first unique coords: %s',
                            sample_name, sample_coords_list
                        )
                except IndexError:
                    logging.debug(
                        'Error while scanning sample position files: %s',
                        traceback.format_exc()
                    )
        except subprocess.CalledProcessError:
            # If something unforeseen goes wrong, the traceback will be
            # printed to screen. We then add the sample to the report with a
            # note that it failed.
            multi_positions = 0
            genus = 'Error processing sample'
            write_output(
                output_report=os.path.join(
                    args.output_name, 'confindr_report.tsv'
                ),
                sample_name=sample_name,
                multi_positions=multi_positions,
                genus=genus,
                total_gene_length=0,
                database_download_date='ND',
            )
            logging.warning(
                'Encountered error when attempting to run ConFindr on sample '
                '%s. Skipping...',
                sample_name,
            )
            logging.warning(
                'Error encountered was:\n%s', traceback.format_exc()
            )
            if args.keep_files is False:
                shutil.rmtree(os.path.join(args.output_name, sample_name))

    # Clean up temporary directories if requested.
    if args.keep_files is False and args.tmp is not None:
        if os.path.isdir(args.tmp):
            shutil.rmtree(args.tmp)

    # Finished all samples.
    logging.info('Contamination detection complete!')


def main() -> None:
    """
    Command-line entrypoint: parse arguments and start analysis.

    This builds the CLI parser, validates user options, configures logging,
    and delegates to :func:`confindr` to run the analysis for each sample.

    Returns:
        None
    """
    # Get CPU count for defaults
    cpu_count = multiprocessing.cpu_count()

    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Define arguments
    parser.add_argument(
        '-i', '--input_directory',
        type=str,
        required=True,
        help='Folder that contains fastq files you want to check for '
        'contamination. Will find any file that contains .fq or .fastq in '
        'the filename.'
    )
    parser.add_argument(
        '-o', '--output_name',
        type=str,
        required=True,
        help='Base name for output/temporary directories.'
    )
    parser.add_argument(
        '-d', '--databases',
        type=str,
        default=os.environ.get(
            'CONFINDR_DB',
            os.path.expanduser('~/.confindr_db')
        ),
        help=(
            'Databases folder. To download these, you will need to get '
            'access to the rMLST databases. For complete instructions see '
            'https://olc-bioinformatics.github.io/ConFindr/install/'
            '#downloading-confindr-databases'
        ),
    )
    parser.add_argument(
        '--rmlst',
        default=False,
        action='store_true',
        help=(
            'Activate to prefer using rMLST databases over core-gene '
            'derived databases. By default, ConFindr will use '
            'core-gene derived databases where available.'
        ),
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=cpu_count,
        help=(
            'Number of threads to use for parallel analysis. Larger values '
            'can reduce runtime but increase CPU and memory usage.'
        )
    )
    parser.add_argument(
        '-tmp', '--tmp',
        type=str,
        help=(
            'Temporary directory for writing genus-specific database files '
            'when the primary database folder is not writable. Useful for '
            'read-only database mounts or restricted environments.'
        ),
    )
    parser.add_argument(
        '-k', '--keep_files',
        default=False,
        action='store_true',
        help=(
            'Keep intermediate files produced during analysis for debugging '
            'or inspection. By default, these files are removed when the '
            'run completes.'
        ),
    )
    parser.add_argument(
        '-q', '--quality_cutoff',
        type=int,
        default=20,
        help=(
            'Base quality threshold used for read trimming. Defaults to 20. '
            'High-quality SNV support is controlled by --min_quality.'
        ),
    )
    parser.add_argument(
        '--min_quality',
        type=int,
        default=15,
        help=(
            'Minimum base quality required to count a base as SNV support. '
            'This filters low-quality bases at individual positions without '
            'discarding the read. Default is 15.'
        )
    )
    parser.add_argument(
        '--use-prob-scoring',
        action='store_true',
        default=False,
        help=(
            'Use probabilistic scoring by requiring statistically supported '
            'positions rather than raw SNP counts. This mode is more '
            'conservative and is intended for higher-confidence contamination '
            'detection.'
        ),
    )
    parser.add_argument(
        '--score-threshold',
        type=float,
        default=2.0,
        help=(
            'Minimum number of statistically supported positions required '
            'to call a sample contaminated in probabilistic mode. '
            'A higher value makes contamination calls more stringent. '
            'Default is 2.0.'
        ),
    )
    parser.add_argument(
        '-b', '--base_cutoff',
        type=int,
        default=3,
        help=(
            'Minimum number of supporting bases required to consider an '
            'alternate allele for SNV calling. This value is adjusted by '
            'the pipeline based on gene-specific data quality and coverage.'
            ' Default is 3.'
        ),
    )
    parser.add_argument(
        '-bf', '--base_fraction_cutoff',
        type=float,
        default=0.05,
        help=(
            'Minimum fraction of the usable pileup depth that must support '
            'an alternate allele for it to be considered a candidate SNV. '
            'This helps avoid low-frequency noise in high-depth samples. '
            'Default is 0.05.'
        ),
    )
    parser.add_argument(
        '-e', '--error_cutoff',
        type=float,
        default=1.0,
        help=(
            'Error cutoff used when computing dynamic base support thresholds.'
            ' Lower values make SNV calls more conservative. Default is 1.0%%.'
        ),
    )
    parser.add_argument(
        '--max-expected-positions',
        type=float,
        default=0.001,
        help=(
            'Maximum expected number of false-positive positions allowed '
            'per gene when computing a dynamic base support cutoff. '
            'Smaller values make the cutoff stricter; set to 0 to disable '
            'dynamic tightening entirely.'
        ),
    )
    parser.add_argument(
        '--downsample_depth',
        type=_valid_downsample_depth,
        default=None,
        metavar='DEPTH',
        help=(
            'Approximate target coverage depth for downsampling reads before '
            'analysis. Useful to reduce runtime or make comparisons more '
            'consistent across samples with very high coverage.'
        ),
    )
    parser.add_argument(
        '--subreplicates',
        type=int,
        default=1,
        help=(
            'Number of independent downsample replicates to run. Default is 1 '
            '(no replicates). If >1, each replicate is downsampled separately '
            'and candidate positions are reported only when enough replicates '
            'agree.'
        ),
    )
    parser.add_argument(
        '--subreplicate-seed',
        type=int,
        default=None,
        help=(
            'Optional seed to make downsample replicates deterministic. '
            'When provided, replicate i uses seed+i so results are repeatable.'
        ),
    )
    parser.add_argument(
        '--subreplicate-consensus',
        type=float,
        default=0.5,
        help=(
            'Consensus fraction of replicates that must support a multibase '
            'position before it is reported. Higher values require more '
            'agreement across replicates. Default is 0.5 (majority).'
        ),
    )
    parser.add_argument(
        '-fid', '--forward_id',
        type=str,
        default='_R1',
        help='Identifier for forward reads.'
    )
    parser.add_argument(
        '-rid', '--reverse_id',
        type=str,
        default='_R2',
        help='Identifier for reverse reads.'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=__version__
    )
    parser.add_argument(
        '-dt', '--data_type',
        choices=['Illumina', 'Nanopore'],
        default='Illumina',
        help=(
            'Type of input data. Default is Illumina, but Nanopore is '
            'also supported (experimental).'
        ),
    )
    parser.add_argument(
        '-Xmx', '--Xmx',
        type=str,
        default=recommend_xmx(),
        help=(
            'Very occasionally, parts of the pipeline that use BBMap may '
            'fail to reserve memory correctly. Use this option to override '
            'automatic memory reservation, e.g. -Xmx 20g or -Xmx 800m. '
            'Default is 80%% of the available virtual memory'
        ),
    )
    parser.add_argument(
        '-cgmlst', '--cgmlst',
        type=str,
        help=(
            'Path to a cgMLST database to use instead of the default rMLST '
            'database. Sequences should have headers like '
            '>genename_allelenumber. Clustering with CD-HIT is recommended '
            'for speed. This is experimental; interpret results with care.'
        ),
    )
    parser.add_argument(
        '--fasta',
        default=False,
        action='store_true',
        help='If activated, will look for FASTA files instead of FASTQ for '
             'unpaired reads.',
    )
    parser.add_argument(
        '--contig-chunk-multiplier',
        type=int,
        default=3,
        help=(
            'Controls how many chunks are created per thread. The total number'
            ' of chunks is threads * multiplier. Larger values give finer '
            'workload distribution but increase scheduling overhead. '
            'Default is 3.'
        )
    )
    parser.add_argument(
        '--contig-chunk-bases',
        type=int,
        default=200000,
        help=(
            'Maximum number of bases assigned to a single chunk. Contigs '
            'larger than this are split into smaller subranges, improving '
            'parallelism for large targets. Smaller values increase overhead. '
            'Default is 200000.'
        )
    )
    parser.add_argument(
        '-verbosity', '--verbosity',
        choices=['debug', 'info', 'warning'],
        default='info',
        help=(
            'Amount of output you want printed to the screen. Defaults '
            'to info (sensible for most users).'
        ),
    )
    parser.add_argument(
        '-m', '--min_matching_hashes',
        default=150,
        type=int,
        help=(
            'Minimum number of matching hashes in a MASH screen required '
            'for a genus to be considered present in a sample. Default '
            'is 150'
        ),
    )
    parser.add_argument(
        '--max_expected_positions',
        default=None,
        type=float,
        help=(
            'Optional maximum expected positions per gene for probabilistic '
            'scoring. If specified, samples with expected positions above '
            'this threshold will be rejected.'
        ),
    )
    args = parser.parse_args()

    # Setup the logger
    fmt = '%(asctime)s %(message)s'
    datefmt = '%Y-%m-%d %H:%M:%S'
    level = {
        'info': 'INFO',
        'debug': 'DEBUG',
        'warning': 'WARNING'
    }.get(args.verbosity, 'INFO')

    if coloredlogs is not None:
        coloredlogs.install(level=level, fmt=fmt, datefmt=datefmt)
    else:
        logging.basicConfig(level=level, format=fmt, datefmt=datefmt)

    logging.info(
        'Welcome to %s! Beginning analysis of your samples...',
        __version__,
    )
    logging.debug(
        'Parsed command-line arguments: %s',
        args
    )

    # Propagate chunking options into methods module
    try:
        methods.CONTIG_CHUNK_MULTIPLIER = int(args.contig_chunk_multiplier)
        methods.CONTIG_CHUNK_MAX_BASES = int(args.contig_chunk_bases)
    except (ValueError, TypeError) as exc:
        logging.debug(
            'Could not configure contig chunking parameters from CLI; using '
            'defaults: %s',
            exc
        )

    # Run ConFindr with the parsed arguments
    confindr(
        args=args,
    )


if __name__ == '__main__':
    main()
