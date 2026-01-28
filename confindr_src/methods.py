#!/usr/bin/env python3
"""
Core ConFindr methods.

This module implements the low-level helpers used by ConFindr's
contamination-detection pipeline (I/O, BAM parsing, statistics and
reporting helpers).
"""

# Standard library imports
from glob import glob
from itertools import chain
from statistics import mean
import traceback
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
)
import csv
import gzip
import logging
import math
import multiprocessing
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request

# Third-party imports
from Bio import SeqIO
from pysam.utils import SamtoolsError
from scipy.stats import (
    betabinom,
    chi2
)
import numpy as np
import pkg_resources
import pysam

# Local imports
from confindr_src.wrappers import (
    bbtools,
    mash,
)


def download_mash_sketch(
    *,  # Enforce keyword arguments
    output_folder: str
) -> None:
    """
    Download the prebuilt RefSeq mash sketch.

    Args:
        output_folder: Directory where the sketch will be saved.

    Returns:
        None
    """
    logging.info('Downloading mash refseq sketch...')
    urllib.request.urlretrieve(
        'https://github.com/OLC-Bioinformatics/ConFindr/raw/master/'
        'refseq_sketch/refseq.msh',
        os.path.join(output_folder, 'refseq.msh')
    )


def download_cgmlst_derived_data(
    *,  # Enforce keyword arguments
    output_folder: str
) -> None:
    """
    Download and extract cgMLST-derived databases.

    The function downloads a tarball containing precomputed databases and
    extracts it into ``output_folder``. The temporary tarball is removed
    after successful extraction.

    Args:
        output_folder: Directory to download and extract files into.

    Returns:
        None
    """
    logging.info(
        'Downloading cgMLST-derived data for Escherichia, Salmonella, '
        'and Listeria...'
    )
    dest = os.path.join(output_folder, 'confindr_db.tar.gz')
    urllib.request.urlretrieve(
        'https://ndownloader.figshare.com/files/14771267',
        dest
    )
    try:
        with tarfile.open(dest) as tar:
            tar.extractall(path=output_folder)
    finally:
        try:
            os.remove(dest)
        except OSError:
            # Best-effort cleanup; ignore failures to remove
            pass
    index(
        output_folder=output_folder,
        genera=['Escherichia', 'Listeria', 'Salmonella'],
        cgderived=True
    )


def index(
    *,  # Enforce keyword arguments
    output_folder: str,
    genera: List[str],
    cgderived: bool = False
) -> None:
    """Ensure genus-specific databases exist and index them.

    Args:
        output_folder: Directory where genus databases live or will be
            created.
        genera: Sequence of genus names to process.
        cgderived: If True, prefer *_db_cgderived.fasta filenames.

    Returns:
        None
    """
    for predominant_genus in genera:
        if cgderived:
            sample_database = os.path.join(
                output_folder,
                f'{predominant_genus}_db_cgderived.fasta'
            )
        else:
            sample_database = os.path.join(
                output_folder,
                f'{predominant_genus}_db.fasta'
            )

        if not os.path.isfile(sample_database):
            if (
                os.path.isfile(
                    os.path.join(
                        output_folder,
                        'gene_allele.txt'
                    )
                ) and os.path.isfile(
                    os.path.join(
                        output_folder,
                        'rMLST_combined.fasta'
                    )
                )
            ):
                logging.info(
                    'Setting up rMLST genus-specific database for genus %s...',
                    predominant_genus
                )
                allele_list = find_genus_specific_allele_list(
                    profiles_file=os.path.join(
                        output_folder, 'gene_allele.txt'
                    ),
                    target_genus=predominant_genus
                )
                # Create the allele-specific database
                setup_allelespecific_database(
                    fasta_file=sample_database,
                    database_folder=output_folder,
                    allele_list=allele_list
                )
        # Perform the necessary samtools and KMA indexing
        index_databases(sample_database=sample_database)


def run_cmd(
    *,  # Enforce keyword arguments
    cmd: str
) -> Tuple[str, str]:
    """
    Run a shell command and return its stdout and stderr.

    Args:
        cmd: The command to run (as would be typed in the shell).

    Returns:
        A tuple (stdout, stderr) as decoded UTF-8 strings.

    Raises:
        subprocess.CalledProcessError: If the command exits with a non-zero
            return code. `output` and `stderr` attributes will contain the
            captured stdout/stderr.
    """
    p = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    out, err = p.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if p.returncode != 0:
        raise subprocess.CalledProcessError(
            p.returncode,
            cmd,
            output=out,
            stderr=err
        )

    return out, err


def write_to_logfile(
    *,  # Enforce keyword arguments
    logfile: str,
    out: str,
    err: str,
    cmd: str
) -> None:
    """
    Append command, stdout and stderr to a logfile.

    Args:
        logfile: Path to file to write output to.
        out: Stdout of program called, as a string.
        err: Stderr of program called, as a string.
        cmd: Command that was used.
    """
    with open(logfile, 'a+', encoding='utf-8') as outfile:
        outfile.write(f'Command used: {cmd}\n\n')
        outfile.write(f'STDOUT: {out}\n\n')
        outfile.write(f'STDERR: {err}\n\n')


def dependency_check(
    *,  # Enforce keyword arguments
    dependency: str
) -> bool:
    """
    Check whether a command-line dependency is available on PATH.

    Args:
        dependency: The executable name to check for (e.g. 'blastn').

    Returns:
        True if the executable is found on PATH, False otherwise.
    """
    return shutil.which(dependency) is not None


def find_paired_reads(
    *,  # Enforce keyword arguments
    fastq_directory: str,
    forward_id: str = '_R1',
    reverse_id: str = '_R2'
) -> List[List[str]]:
    """
    Find paired FASTQ files in a directory.

    Args:
        fastq_directory: Path containing FASTQ files.
        forward_id: Identifier substring for forward reads (default: '_R1').
        reverse_id: Identifier substring for reverse reads (default: '_R2').

    Returns:
        A list of [forward, reverse] pairs (paths as strings).
    """
    # Define list to hold pairs
    pair_list: List[List[str]] = []

    # Find all FASTQ files in the directory
    fastq_files = glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in sorted(fastq_files):
        # Check to see if the file is a forward read and that the corresponding
        # reverse read exists
        if forward_id in name and os.path.isfile(
            name.replace(
                forward_id,
                reverse_id
            )
        ):
            pair_list.append([name, name.replace(forward_id, reverse_id)])

    return pair_list


def find_unpaired_reads(
    *,  # Enforce keyword arguments
    fastq_directory: str,
    forward_id: str = '_R1',
    reverse_id: str = '_R2',
    find_fasta: bool = False
) -> List[List[str]]:
    """
    Find unpaired FASTQ (or FASTA) files in a directory.

    Args:
        fastq_directory: Path to directory containing reads.
        forward_id: Identifier substring for forward reads (default: '_R1').
        reverse_id: Identifier substring for reverse reads (default: '_R2').
        find_fasta: If True, search for FASTA (*.f*a*). Otherwise search for
        FASTQ.

    Returns:
        A list of single-element lists containing paths to unpaired read files.
    """
    # Define list to hold unpaired reads
    read_list: List[List[str]] = []

    # Find all FASTQ or FASTA files in the directory
    if find_fasta is False:
        sequence_files = glob(os.path.join(fastq_directory, '*.f*q*'))
    else:
        # Find FASTA files
        sequence_files = glob(os.path.join(fastq_directory, '*.f*a*'))

    # Iterate through files, adding them to our list of unpaired reads if:
    for name in sorted(sequence_files):
        # 1) They don't have the forward identifier or the reverse identifier
        # in their name.
        if forward_id not in name and reverse_id not in name:
            read_list.append([name])

        # 2) They have forward but the reverse isn't there.
        elif forward_id in name and not os.path.isfile(
            name.replace(
                forward_id,
                reverse_id
            )
        ):
            read_list.append([name])

        # 3) They have reverse but the forward isn't there.
        elif reverse_id in name and not os.path.isfile(
            name.replace(
                reverse_id,
                forward_id
            )
        ):
            read_list.append([name])

    return read_list


def find_genus_specific_allele_list(
    *,  # Enforce keyword arguments
    profiles_file: str,
    target_genus: str
) -> List[str]:
    """
    Return allele list for a target genus from a profiles file.

    Args:
        profiles_file: Path to profiles file containing genus:allele1,
        allele2,... lines.
        target_genus: Genus for which alleles should be returned.

    Returns:
        A list of allele identifiers for the target genus.
    """
    # Initialize empty allele list
    alleles: List[str] = []

    # Read through the profiles file to find alleles for the target genus
    with open(profiles_file, encoding='utf-8') as f:
        lines = f.readlines()

    # Parse the lines to find the target genus
    for line in lines:
        line = line.rstrip()
        genus = line.split(':')[0]

        # If this line corresponds to the target genus, extract alleles
        if genus == target_genus:
            alleles = line.split(':')[1].split(',')[:-1]

    return alleles


def setup_allelespecific_database(
    *,  # Enforce keyword arguments
    fasta_file: str,
    database_folder: str,
    allele_list: List[str]
) -> None:
    """
    Create a genus-specific FASTA file containing only the alleles listed.

    Args:
        fasta_file: Path to FASTA file to write allele-specific database to.
        database_folder: Path where rMLST_combined.fasta can be found.
        allele_list: List of allele identifiers to include.
    """
    # Create an index of the rMLST combined FASTA file
    rmlst_index = SeqIO.index(
        os.path.join(
            database_folder,
            'rMLST_combined.fasta'
        ),
        'fasta'
    )

    # Create a list to hold the sequences to write
    seqs = []

    # Iterate through the allele list, adding sequences to the seqs list
    for s in allele_list:
        try:
            seqs.append(rmlst_index[s])
        except KeyError:
            logging.warning(
                'Tried to add %s to allele-specific database, but could not '
                'find it.', s
            )
    try:
        SeqIO.write(seqs, fasta_file, 'fasta')
    except FileNotFoundError:
        pass


def find_cross_contamination(
    *,  # Enforce keyword arguments
    databases: str,
    reads: Any,
    sample_name: str,
    tmpdir: str = 'tmp',
    log: str = 'log.txt',
    threads: int = 1,
    min_matching_hashes: int = 40
) -> str:
    """
    Uses mash to find out whether or not a sample has more than one genus
    present, indicating cross-contamination.

    Args:
        databases: Path to folder containing mash sketch database (refseq.msh).
        reads: Either a string path to single-end reads, or a list of two
        string paths for paired-end reads.
        sample_name: Sample name to use for temporary files.
        tmpdir: Path to temporary directory to use.
        log: Path to logfile to write mash output to.
        threads: Number of threads to use.
        min_matching_hashes: Minimum number of matching hashes to consider a
        genus present.

    Returns:
        A string representing the genera present. If only one genus is found,
        the string is that genus. If no genera are found, the string is 'ND'.
        If more than one genus is found, the string is a list of genera
        present, separated by colons.
    """
    # Initialize empty list to hold genera present
    genera_present = []

    # Only run the MASH analyses if the screen.tab output file does not
    # already exist
    screen_file = os.path.join(
        tmpdir,
        f'{sample_name}_screen.tab'
    )

    # Run mash screen if output file doesn't already exist
    if not os.path.isfile(screen_file):
        # If reads is a string, it's unpaired. If it's a list, it's paired
        if isinstance(reads, str):
            out, err, cmd = mash.screen(
                f'{databases}/refseq.msh',
                reads,
                threads=threads,
                w='',
                i='0.85',
                output_file=screen_file,
                returncmd=True
            )
        else:
            out, err, cmd = mash.screen(
                f'{databases}/refseq.msh',
                reads[0],
                reads[1],
                threads=threads,
                w='',
                i='0.85',
                output_file=screen_file,
                returncmd=True
            )

        # Write mash output to logfile
        write_to_logfile(
            logfile=log,
            out=out,
            err=err,
            cmd=cmd
        )

    # Read mash screen output
    screen_output = mash.read_mash_screen(
        os.path.join(
            tmpdir,
            f'{sample_name}_screen.tab'
        )
    )

    # Parse mash screen output to find genera present
    for item in screen_output:
        # Extract genus from the query ID
        mash_genus = item.query_id.split('/')[-3]

        # Convert Shigella to Escherichia
        if 'Shigella' in mash_genus:
            mash_genus = 'Escherichia'

        # Extract number of matching hashes
        matching_hashes = int(item.shared_hashes.split('/')[0])

        # Only add the genus to the genera_present list of the number of
        # matching hashes exceeds the cutoff
        if matching_hashes >= min_matching_hashes:
            if mash_genus not in genera_present:
                genera_present.append(mash_genus)

    # Format the output string depending on how many genera were found
    if len(genera_present) == 1:
        genera_present = genera_present[0]
    elif len(genera_present) == 0:
        genera_present = 'ND'
    else:
        # Concatenate multiple genera with colons
        tmpstr = ''
        for mash_genus in genera_present:
            tmpstr += mash_genus + ':'
        genera_present = tmpstr[:-1]

    return genera_present


def number_of_bases_above_threshold(
    *,  # Enforce keyword arguments
    high_quality_base_count: Dict[str, int],
    base_count_cutoff: int = 2,
    base_fraction_cutoff: Optional[float] = None
) -> int:
    """
    Return number of bases that meet configured high-quality thresholds.

    Args:
        high_quality_base_count: Mapping of base->count at a position.
        base_count_cutoff: Absolute count cutoff for a base to be considered
        high-quality.
        base_fraction_cutoff: Optional fraction of total coverage required for
        base to count.

    Returns:
        The number of bases that meet the threshold (integer).
    """
    # Dictionary comprehension where values are True or False for each base
    # depending on whether the count meets the threshold.
    # Method differs depending on whether absolute or fraction cutoff is
    # specified
    if base_fraction_cutoff:
        # Calculate total high-quality base count for fraction comparison
        total_hq_base_count = sum(high_quality_base_count.values())

        # Calculate which bases meet both the fraction and count cutoffs
        bases_above_threshold = {
            base:
                float(count) / total_hq_base_count >= base_fraction_cutoff
                and count >= base_count_cutoff for (base, count)
                in high_quality_base_count.items()
        }
    else:
        bases_above_threshold = {
            base:
                count >= base_count_cutoff for (base, count)
                in high_quality_base_count.items()
            }

    # True is equal to 1 so sum of the number of Trues in the
    # bases_above_threshold dict is the number of bases passing threshold
    return sum(bases_above_threshold.values())


def parse_bam(
    *,  # Enforce keyword arguments
    bamfile_name: str,
    contig_name: str,
    pysam_fasta: Any
) -> Tuple[Any, Any]:
    """
    Open a BAM file with pysam and produce a pileup iterator for a contig.

    Args:
        bamfile_name: Path to sorted, indexed BAM file.
        contig_name: Contig/sequence name to produce pileup for.
        pysam_fasta: A pysam.FastaFile (or compatible) used for reference base.

    Returns:
        A tuple (bamfile, pileup) where bamfile is the opened AlignmentFile and
        pileup is the iterator produced by pysam's pileup() for the contig.
    """
    # Load the sorted BAM-formatted file using pysam
    bamfile = pysam.AlignmentFile(bamfile_name, 'rb')

    # These parameters seem to be fairly undocumented with pysam, but I think
    # that they should make the output that I'm getting to match up with
    # what I'm seeing in Tablet.
    pileup = bamfile.pileup(
        contig_name,
        stepper='samtools',
        ignore_orphans=False,
        fastafile=pysam_fasta,
        min_base_quality=0
    )

    return bamfile, pileup


def characterise_read(
    *,  # Enforce keyword arguments
    column: Any,
    reference_sequence: str,
    fastq_records: Dict[str, Any],
    quality_cutoff: int,
    fasta: bool = False,
    nanopore: bool = False
) -> Dict[str, Any]:
    """
    Extract read-level characteristics from a pileup column.

    Args:
        column: Pysam pileup column object for the position.
        reference_sequence: The full reference sequence string for the contig.
        fastq_records: Mapping of read name -> SeqRecord with quality info.
        quality_cutoff: Minimum base quality to be considered high-quality.

    Keyword Args:
        fasta: If True, operate in FASTA-only mode (no quality checks).
        nanopore: If True, use nanopore-specific heuristics.

    Returns:
        A dict containing read characterisations (bases, qualities, mapping
        qualities, strand information, etc.)."
    """
    # Initialise a dictionary to store the base details
    filtered_read_dict = {
        'congruent_SNV': {},
        'congruent_ref': {},
        'forward_SNV_reverse_SNV1': {},
        'reverse_SNV_forward_SNV1': {},
        'forward_SNV_reverse_ref': {},
        'reverse_SNV_forward_ref': {},
        'forward_SNV_reverse_UM_QF': {},
        'forward_ref_reverse_UM_QF': {},
        'forward_quality_filtered': {},
        'reverse_SNV_forward_UM_QF': {},
        'reverse_ref_forward_UM_QF': {},
        'reverse_quality_filtered': {}
    }

    # Initialise a dictionary to store the details parsed from the pileup
    unfiltered_read_details = {}

    # Initialise a list to store all the phred scores for the bases passing
    # filter in the column
    qualities = []

    # Initialise a dict to store per-base support metadata for statistical
    # tests
    base_support = {}

    # Iterate through every read present in the column of the pileup
    for read in column.pileups:
        # Not entirely sure why this is sometimes None, but it causes bad stuff
        if read.query_position is not None:
            #  Initialise the read name in the dictionary as required
            if read.alignment.qname not in unfiltered_read_details:
                unfiltered_read_details[read.alignment.qname] = {}

            # Extract the sequence of the base in the read
            query_base = read.alignment.query_sequence[read.query_position]

            # Extract the sequence of the base in the reference gene
            ref_base = reference_sequence[column.pos]

            # Create a boolean of whether the query base matches the reference
            # base (is it a SNV?)
            match = query_base == ref_base

            # Read names in the pileup have the direction removed - add it
            # back for future parsing. Not an issue for FASTA files
            if not fasta and not nanopore:
                if read.alignment.is_read1:
                    if read.alignment.qname.split(' ')[0].endswith('/1'):
                        read_name = read.alignment.qname.split(' ')[0]
                    else:
                        read_name = read.alignment.qname.split(' ')[0] + '/1'
                else:
                    if read.alignment.qname.split(' ')[0].endswith('/2'):
                        read_name = read.alignment.qname.split(' ')[0]
                    else:
                        read_name = read.alignment.qname.split(' ')[0] + '/2'
            else:
                read_name = read.alignment.qname

            # Extract the phred quality score from the FASTQ records
            quality = fastq_records[read_name].letter_annotations[
                "phred_quality"
            ][read.query_position]

            # Initialise a dictionary to store the range
            range_dict = {}

            # Determine whether there are SNVs clustered together - they will
            # be discarded from the analysis. Iterate through a range of the
            # five positions preceding, and five subsequent positions
            for iterator, contig_pos in enumerate(
                chain(
                    range(
                        column.pos - 5, column.pos
                    ),
                    range(
                        column.pos + 1, column.pos + 6
                    )
                )
            ):
                # Ensure that the contig position being examined isn't beyond
                # the length of the gene
                if 0 <= contig_pos < len(reference_sequence) - 1:
                    # Calculate the read position corresponding to the current
                    # column position
                    read_pos = list(
                        chain(
                            range(
                                read.query_position - 5, read.query_position
                            ),
                            range(
                                read.query_position + 1,
                                read.query_position + 6
                            )
                        )
                    )[iterator]
                    # Ensure that read position being examined isn't beyond
                    # the length of the read
                    if 0 <= read_pos < read.alignment.query_alignment_end - 1:
                        # Populate the dictionary with the calculated positions
                        range_dict[contig_pos] = read_pos

            # Initialise a boolean of whether the current base passes filters
            # and should be added to the dictionary
            add_base = True

            # Iterate through the range_dict to extract the sequence of the
            # bases in the range for both the read and the reference sequence
            for contig_pos, read_pos in range_dict.items():
                try:
                    # Extract the sequence of the base
                    reference_base = reference_sequence[contig_pos]
                    base = read.alignment.query_sequence[read_pos]
                    # If any of the downstream or upstream bases do not match,
                    # set the boolean to False
                    if not match and reference_base != base:
                        add_base = False
                except KeyError:
                    pass

            # Populate the dictionary only if there are no other SNVs within
            # five downstream and five upstream bases
            if add_base:
                unfiltered_read_details[
                    read.alignment.qname
                ][read.alignment.is_read1] = {
                    'mate unmapped': read.alignment.mate_is_unmapped,
                    'paired': read.alignment.is_paired,
                    'forward': read.alignment.is_read1,
                    'reverse': read.alignment.is_read2,
                    'match': match,
                    'qbase': query_base,
                    'rbase': ref_base,
                    'qual': quality,
                    'pos': column.pos,
                    'gene': column.reference_name
                }

                # Add the quality of the base
                if quality >= quality_cutoff:
                    qualities.append(quality)
                    # record per-base support qualities, mapping qualities and
                    # strand counts for later statistical tests
                    if query_base not in base_support:
                        base_support[query_base] = {
                            'quals': [],
                            'forward_quals': [],
                            'reverse_quals': [],
                            'forward': 0,
                            'reverse': 0,
                            'mapqs': [],
                            'forward_mapqs': [],
                            'reverse_mapqs': [],
                            'positions': [],
                            'forward_positions': [],
                            'reverse_positions': []
                        }
                    base_support[query_base]['quals'].append(quality)

                    # Mapping quality and read position
                    try:
                        mapq = int(read.alignment.mapping_quality)
                    except TypeError:
                        mapq = None

                    # Add mapping quality if it exists
                    if mapq is not None:
                        base_support[query_base]['mapqs'].append(mapq)

                    # Add strand and position information
                    base_support[query_base]['positions'].append(
                        read.query_position
                    )

                    # Strand-specific information
                    if read.alignment.is_read1:
                        base_support[query_base]['forward_quals'].append(
                            quality
                        )

                        # Increment forward strand count
                        base_support[query_base]['forward'] += 1

                        # Add mapping quality if it exists
                        if mapq is not None:
                            base_support[query_base]['forward_mapqs'].append(
                                mapq
                            )

                        # Add read position
                        base_support[query_base]['forward_positions'].append(
                            read.query_position
                        )
                    else:
                        # Reverse strand
                        base_support[query_base]['reverse_quals'].append(
                            quality
                        )

                        # Increment reverse strand count
                        base_support[query_base]['reverse'] += 1

                        # Add mapping quality if it exists
                        if mapq is not None:
                            base_support[query_base]['reverse_mapqs'].append(
                                mapq
                            )

                        # Add read position
                        base_support[query_base]['reverse_positions'].append(
                            read.query_position
                        )

    # Parse the unfiltered reads to characterise the bases
    for _, dir_dict in unfiltered_read_details.items():
        # Check to see if paired reads are present at this position
        if len(dir_dict) > 1:
            # SNV in both forward and reverse reads
            if not dir_dict[True]['match'] and not dir_dict[False]['match']:
                # Same SNV sequence - don't quality filter as the reads
                # agreeing acts as a quality check
                if dir_dict[True]['qbase'] == dir_dict[False]['qbase']:
                    # Forward and reverse SNV - add two matches (for forward
                    # and reverse reads)
                    if (
                        dir_dict[True]['qbase']
                        not in filtered_read_dict['congruent_SNV']
                    ):
                        filtered_read_dict[
                            'congruent_SNV'
                        ][dir_dict[True]['qbase']] = 2
                    else:
                        filtered_read_dict[
                            'congruent_SNV'
                        ][dir_dict[True]['qbase']] += 2

                # Different SNV sequences
                else:
                    # Both SNV sequences pass quality
                    if (
                        dir_dict[True]['qual'] >= quality_cutoff
                        and dir_dict[False]['qual'] >= quality_cutoff
                    ):
                        # Forward SNV1
                        if (
                            dir_dict[True]['qbase']not in
                            filtered_read_dict['forward_SNV_reverse_SNV1']
                        ):
                            filtered_read_dict[
                                'forward_SNV_reverse_SNV1'
                            ][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_SNV_reverse_SNV1'
                            ][dir_dict[True]['qbase']] += 1

                        # Reverse SNV2
                        if (
                            dir_dict[False]['qbase'] not in
                            filtered_read_dict['reverse_SNV_forward_SNV1']
                        ):
                            filtered_read_dict[
                                'reverse_SNV_forward_SNV1'
                            ][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_SNV_forward_SNV1'
                            ][dir_dict[False]['qbase']] += 1

                    # Only the forward reads pass quality
                    elif dir_dict[True]['qual'] >= quality_cutoff:
                        # Forward SNV reverse QF
                        if (
                            dir_dict[True]['qbase'] not in
                            filtered_read_dict['forward_SNV_reverse_UM_QF']
                        ):
                            filtered_read_dict[
                                'forward_SNV_reverse_UM_QF'
                            ][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_SNV_reverse_UM_QF'
                            ][dir_dict[True]['qbase']] += 1

                        # Reverse QF
                        if (
                            dir_dict[False]['qbase'] not in
                            filtered_read_dict['reverse_quality_filtered']
                        ):
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[False]['qbase']] += 1

                    # Only the reverse reads pass quality
                    elif dir_dict[False]['qual'] >= quality_cutoff:
                        # Reverse SNV forward QF
                        if (
                            dir_dict[False]['qbase'] not in
                            filtered_read_dict['reverse_SNV_forward_UM_QF']
                        ):
                            filtered_read_dict[
                                'reverse_SNV_forward_UM_QF'
                            ][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_SNV_forward_UM_QF'
                            ][dir_dict[False]['qbase']] += 1

                        # Forward QF
                        if (
                            dir_dict[True]['qbase'] not in
                            filtered_read_dict['forward_quality_filtered']
                        ):
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[True]['qbase']] += 1

                    # Neither forward nor reverse reads pass quality
                    else:
                        # Forward QF
                        if (
                            dir_dict[True]['qbase'] not in
                            filtered_read_dict['forward_quality_filtered']
                        ):
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[True]['qbase']] += 1

                        # Reverse QF
                        if (
                            dir_dict[False]['qbase'] not in
                            filtered_read_dict['reverse_quality_filtered']
                        ):
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[False]['qbase']] += 1

            # SNV in forward read only
            elif not dir_dict[True]['match'] and dir_dict[False]['match']:
                # Since only the forward read supports the SNV, quality filter
                if dir_dict[True]['qual'] >= quality_cutoff:
                    if (
                        dir_dict[True]['qbase'] not in
                        filtered_read_dict['forward_SNV_reverse_ref']
                    ):
                        filtered_read_dict[
                            'forward_SNV_reverse_ref'
                        ][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'forward_SNV_reverse_ref'
                        ][dir_dict[True]['qbase']] += 1
                    # Since the reverse base matches the reference,
                    # don't quality filter
                    if (
                        dir_dict[False]['qbase'] not in
                        filtered_read_dict['forward_SNV_reverse_ref']
                    ):
                        filtered_read_dict[
                            'forward_SNV_reverse_ref'
                        ][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'forward_SNV_reverse_ref'
                        ][dir_dict[False]['qbase']] += 1

                # Reverse ref forward QF
                else:
                    if (
                        dir_dict[False]['qbase'] not in
                        filtered_read_dict['reverse_ref_forward_UM_QF']
                    ):
                        filtered_read_dict[
                            'reverse_ref_forward_UM_QF'
                        ][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'reverse_ref_forward_UM_QF'
                        ][dir_dict[False]['qbase']] += 1

            # SNV in reverse read only
            elif dir_dict[True]['match'] and not dir_dict[False]['match']:
                # Quality filter
                if dir_dict[False]['qual'] >= quality_cutoff:
                    if (
                        dir_dict[False]['qbase'] not in
                        filtered_read_dict['reverse_SNV_forward_ref']
                    ):
                        filtered_read_dict[
                            'reverse_SNV_forward_ref'
                        ][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'reverse_SNV_forward_ref'
                        ][dir_dict[False]['qbase']] += 1
                    # Since the forward base matches the reference, don't
                    # quality filter
                    if (
                        dir_dict[True]['qbase'] not in
                        filtered_read_dict['reverse_SNV_forward_ref']
                    ):
                        filtered_read_dict[
                            'reverse_SNV_forward_ref'
                        ][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'reverse_SNV_forward_ref'
                        ][dir_dict[True]['qbase']] += 1

                # Forward ref reverse QF
                else:
                    if (
                        dir_dict[True]['qbase'] not in
                        filtered_read_dict['forward_ref_reverse_UM_QF']
                    ):
                        filtered_read_dict[
                            'forward_ref_reverse_UM_QF'
                        ][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict[
                            'forward_ref_reverse_UM_QF'
                        ][dir_dict[True]['qbase']] += 1

            # Both reads match the reference sequence - don't quality filter,
            # and add two matches (for forward and reverse reads)
            else:
                if (
                    dir_dict[True]['qbase'] not in
                    filtered_read_dict['congruent_ref']
                ):
                    filtered_read_dict[
                        'congruent_ref'
                    ][dir_dict[True]['qbase']] = 2
                else:
                    filtered_read_dict[
                        'congruent_ref'
                    ][dir_dict[True]['qbase']] += 2

        # Either the reads are unpaired, or only a single read aligns to this
        # position on the gene (no overlap)
        else:
            for direction in dir_dict:
                # SNV supported by a single read
                if not dir_dict[direction]['match']:
                    if dir_dict[direction]['qual'] >= quality_cutoff:
                        # Forward
                        if direction:
                            if (
                                dir_dict[direction]['qbase'] not in
                                filtered_read_dict['forward_SNV_reverse_UM_QF']
                            ):
                                filtered_read_dict[
                                    'forward_SNV_reverse_UM_QF'
                                ][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict[
                                    'forward_SNV_reverse_UM_QF'
                                ][dir_dict[direction]['qbase']] += 1

                        # Reverse
                        else:
                            if (
                                dir_dict[direction]['qbase'] not in
                                filtered_read_dict['reverse_SNV_forward_UM_QF']
                            ):
                                filtered_read_dict[
                                    'reverse_SNV_forward_UM_QF'
                                ][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict[
                                    'reverse_SNV_forward_UM_QF'
                                ][dir_dict[direction]['qbase']] += 1
                    else:
                        # Forward QF
                        if (
                            dir_dict[direction]['qbase'] not in
                            filtered_read_dict['forward_quality_filtered']
                        ):
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[direction]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'forward_quality_filtered'
                            ][dir_dict[direction]['qbase']] += 1

                # Match to the reference sequence supported by a single read
                else:
                    if dir_dict[direction]['qual'] >= quality_cutoff:
                        # Forward
                        if direction:
                            if (
                                dir_dict[direction]['qbase'] not in
                                filtered_read_dict['forward_ref_reverse_UM_QF']
                            ):
                                filtered_read_dict[
                                    'forward_ref_reverse_UM_QF'
                                ][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict[
                                    'forward_ref_reverse_UM_QF'
                                ][dir_dict[direction]['qbase']] += 1

                        # Reverse
                        else:
                            if (
                                dir_dict[direction]['qbase'] not in
                                filtered_read_dict['reverse_ref_forward_UM_QF']
                            ):
                                filtered_read_dict[
                                    'reverse_ref_forward_UM_QF'
                                ][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict[
                                    'reverse_ref_forward_UM_QF'
                                ][dir_dict[direction]['qbase']] += 1
                    else:
                        # Reverse QF
                        if (
                            dir_dict[direction]['qbase'] not in
                            filtered_read_dict['reverse_quality_filtered']
                        ):
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[direction]['qbase']] = 1
                        else:
                            filtered_read_dict[
                                'reverse_quality_filtered'
                            ][dir_dict[direction]['qbase']] += 1

    return filtered_read_dict, qualities, base_support


def determine_cutoff(
    *,  # Enforce keyword arguments
    qualities: List[int],
    reference_sequence: str,
    base_cutoff: int,
    error_cutoff: float = 1.0
) -> int:
    """
    Determine the smallest integer k >= base_cutoff such that the probability
    of observing >= k sequencing errors by chance (given per-base qualities)
    is <= alpha, where alpha = (error_cutoff% / 100) / len(reference_sequence).

    Uses exact Poisson-binomial pmf by iterative convolution when depth is
    reasonable; otherwise uses a normal approximation with continuity
    correction.

    Args:
        qualities: Per-base quality scores at the position.
        reference_sequence: Reference sequence string (used for length/GC
        heuristics).
        base_cutoff: Base cutoff parameter (minimum bases required before SNV
        considered).
        error_cutoff: Expected error percentage (as 1.0 for 1%).

    Returns:
        The computed integer cutoff (>= 0).
    """
    # Determine the maximum read length for error percentage calculation
    max_len = max(1, len(reference_sequence))

    # Calculate alpha for error percentage calculation
    alpha = (error_cutoff / 100.0) / max_len

    # Handle edge cases
    if not qualities or base_cutoff < 1:
        return max(1, base_cutoff), 0.0

    # Per-observation error probabilities from Phred Q
    p_list = [10 ** (-q / 10.0) for q in qualities]
    n = len(p_list)

    # Try exact poisson-binomial via iterative convolution if possible and not
    # too large
    if n <= 1000:
        pmf = np.array([1.0], dtype=float)
        for p in p_list:
            pmf = np.convolve(pmf, np.array([1.0 - p, p], dtype=float))

        # Tail probability function: cumulative sum from end:
        # cumulative_sum[k] = P(X>=k)
        cumulative_sum = pmf[::-1].cumsum()[::-1]
        k = max(1, base_cutoff)
        cap = n
        while k <= cap and cumulative_sum[k] > alpha:
            k += 1

        # Calculate error percentage
        tail_prob = cumulative_sum[k] if k <= cap else 0.0
        error_perc = tail_prob * 100.0 * max_len
        return k, error_perc

    # Fallback: normal approximation to Poisson-binomial with continuity
    # correction
    mu = sum(p_list)
    var = sum(p * (1.0 - p) for p in p_list)
    sigma = math.sqrt(var) if var > 0 else 0.0

    def normal_tail(k):
        if sigma == 0.0:
            # degenerate: if mu rounds to < k then probability is 0, else 1
            return 1.0 if mu + 1e-12 >= k else 0.0
        z = (k - 0.5 - mu) / (sigma * math.sqrt(2.0))
        return 0.5 * math.erfc(z)

    k = max(1, base_cutoff)
    # cap to avoid long loops; in practice k won't exceed n
    cap = max(n, base_cutoff + 1000)
    while k <= cap and normal_tail(k) > alpha:
        k += 1
    tail_prob = normal_tail(k) if k <= cap else 0.0
    error_perc = tail_prob * 100.0 * max_len

    return k, error_perc


def _poisson_binomial_pmf_fft(
    *,  # Enforce keyword arguments
    p_list: List[float]
) -> 'np.ndarray':
    """
    Compute Poisson-binomial PMF via divide-and-conquer FFT convolution.
    Returns an array of length n+1 where pmf[j] = P(X==j).

    Args:
        p_list: List of success probabilities for independent Bernoulli trials.

    Returns:
        A numpy array containing the PMF for 0..n successes.
    """
    # convert to numpy array
    p_arr = np.asarray(p_list, dtype=float)
    n = p_arr.size
    if n == 0:
        return np.array([1.0], dtype=float)
    if n == 1:
        return np.array([1.0 - p_arr[0], p_arr[0]], dtype=float)

    # recursive divide-and-conquer
    def _rec(pv):
        m = len(pv)
        if m == 0:
            return np.array([1.0], dtype=float)
        if m == 1:
            return np.array([1.0 - pv[0], pv[0]], dtype=float)
        mid = m // 2
        left = _rec(pv[:mid])
        right = _rec(pv[mid:])
        size = left.size + right.size - 1
        nfft = 1 << ((size - 1).bit_length())
        fa = np.fft.rfft(left, nfft)
        fb = np.fft.rfft(right, nfft)
        fr = fa * fb
        conv = np.fft.irfft(fr, nfft)[:size]

        # numerical rounding cleanup
        conv[conv < 0] = 0.0
        return conv
    pmf = _rec(p_arr.tolist())

    # renormalize to sum to 1
    s = pmf.sum()
    if s > 0:
        pmf = pmf / s
    return pmf


def _beta_binomial_tail(
    *,  # Enforce keyword arguments
    k: int,
    n: int,
    a: float,
    b: float
) -> float:
    """
    Compute tail probability P(X>=k) for Beta-Binomial(n, a, b).
    Prefer SciPy's betabinom.sf when available for numerical stability;
    fall back to log-gamma summation

    Args:
        k: Minimum number of successes to include in tail.
        n: Number of trials.
        a: Alpha parameter of Beta prior.
        b: Beta parameter of Beta prior.

    Returns:
        The tail probability as a float.
    """
    # handle edge cases
    if n <= 0:
        return 0.0 if k > 0 else 1.0

    # use SciPy's betabinom.sf to calculate tail probability
    return float(betabinom.sf(k - 1, n, a, b))


def _estimate_beta_params(
    *,  # Enforce keyword arguments
    p_list: List[float]
) -> Optional[Tuple[float, float]]:
    """
    Estimate Beta distribution (alpha, beta) from probabilities list.

    Args:
        p_list: List of observed probabilities.

    Returns:
        (alpha, beta) tuple, or (None, None) if estimation is not possible
        (e.g. zero variance).
    """
    if not p_list:
        return None, None
    m = float(np.mean(p_list))
    v = float(np.var(p_list, ddof=0))
    # ensure positive variance less than m*(1-m)
    denom = m * (1.0 - m)
    if v <= 0 or v >= denom:
        return None, None
    common = (denom / v) - 1.0
    a = m * common
    b = (1.0 - m) * common
    # ensure >0
    if a <= 0 or b <= 0:
        return None, None
    return a, b


def poisson_binomial_tail(
    *,  # Enforce keyword arguments
    k: int,
    p_list: List[float]
) -> float:
    """
    Tail probability for Poisson-Binomial distribution P(X >= k). Uses
    FFT-based convolution for exact pmf when feasible, otherwise normal
    approximation with continuity correction.

    Args:
        k: Threshold number of successes.
        p_list: List of success probabilities.

    Returns:
        Tail probability (float).
    """
    # Set the value of n from the length of p_list
    n = len(p_list)
    if n == 0:
        return 0.0 if k > 0 else 1.0
    # If n is small-ish, compute exact pmf via FFT divide-and-conquer (fast
    # and accurate)
    pmf = _poisson_binomial_pmf_fft(
        p_list=p_list
    )
    cumsum = pmf[::-1].cumsum()[::-1]

    return float(cumsum[k]) if k <= n else 0.0


def mann_whitney_u_p(
    *,  # Enforce keyword arguments
    x: List[float],
    y: List[float]
) -> float:
    """
    Two-sample Mann-Whitney U test p-value (two-sided) using normal
    approximation with continuity correction. Assigns average ranks for ties
    and applies tie correction to variance.

    Args:
        x, y: Lists of samples.

    Returns:
        Two-sided p-value as float.
    """
    # Initialize sample sizes
    n1 = len(x)
    n2 = len(y)

    # Handle edge cases
    if n1 == 0 or n2 == 0:
        return None

    # Calculate merged list of (value, group) tuples
    merged = [(v, 0) for v in x] + [(v, 1) for v in y]
    merged.sort(key=lambda t: t[0])
    total = n1 + n2

    # assign average ranks for ties
    ranks = [0.0] * total
    i = 0
    while i < total:
        j = i
        while j + 1 < total and merged[j + 1][0] == merged[i][0]:
            j += 1
        avg_rank = sum(range(i + 1, j + 2)) / (j - i + 1)
        for k in range(i, j + 1):
            ranks[k] = avg_rank
        i = j + 1

    # sum ranks for group x
    summed = sum(ranks[idx] for idx, (_, grp) in enumerate(merged) if grp == 0)
    u1 = summed - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0

    # tie correction
    tie_sum = 0.0
    i = 0
    while i < total:
        j = i
        while j + 1 < total and merged[j + 1][0] == merged[i][0]:
            j += 1
        t = j - i + 1
        if t > 1:
            tie_sum += t * (t * t - 1)
        i = j + 1
    denom = 12.0 * total * (total - 1)
    var = n1 * n2 * ((total + 1) - (tie_sum / denom)) / 12.0
    if var <= 0:
        return 1.0
    sigma = math.sqrt(var)
    z = (u1 - mu) / sigma

    return 2.0 * 0.5 * math.erfc(-abs(z) / math.sqrt(2.0))


def fisher_two_sided_p(
    *,  # Enforce keyword arguments
    a: int,
    b: int,
    c: int,
    d: int
) -> float:
    """Two-sided Fisher exact test p-value for a 2x2 contingency table.

    Args:
        a, b, c, d: Table counts arranged as [[a, b], [c, d]].

    Returns:
        Two-sided p-value as float.
    """
    # Hypergeometric probabilities for fixed margins
    def hyper_p(x, m1, n1, k1):
        # prob of x successes in sample size k1 from population with m1
        # successes and n1 failures
        return (
            math.comb(m1, x) * math.comb(n1, k1 - x)
        ) / math.comb(m1 + n1, k1)

    m1 = a + c  # successes in population (col1)
    n1 = b + d  # failures in population (col2)
    k1 = a + b  # sample size (row1)
    if m1 < 0 or n1 < 0 or k1 < 0:
        return 1.0
    # compute observed probability
    try:
        p_obs = hyper_p(a, m1, n1, k1)
    except ValueError:
        return 1.0
    # iterate all possible x and sum probabilities <= p_obs
    lo = max(0, k1 - n1)
    hi = min(k1, m1)
    p_sum = 0.0
    for x in range(lo, hi + 1):
        try:
            p_x = hyper_p(x, m1, n1, k1)
        except ValueError:
            p_x = 0.0
        if p_x <= p_obs + 1e-20:
            p_sum += p_x

    return min(1.0, p_sum)


def benjamini_hochberg(
    *,  # Enforce keyword arguments
    pvals: List[float]
) -> List[float]:
    """
    Compute Benjamini-Hochberg adjusted q-values.

    Args:
        pvals: List of p-values.

    Returns:
        List of q-values (same order as input pvals).
    """
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    qvals = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        unadj = pvals[i]
        adj = min(prev, unadj * n / rank)
        qvals[i] = adj
        prev = adj
    return qvals


def combine_pvalues_fisher(
    *,  # Enforce keyword arguments
    pvals: List[float]
) -> float:
    """
    Combine a list of p-values using Fisher's method. Returns combined p-value

    Args:
        pvals: List of p-values.

    Returns:
        Combined p-value as float.
    """
    if not pvals:
        return None

    # clip p-values to avoid log(0)
    eps = 1e-300
    adj = [max(min(p, 1.0), eps) for p in pvals]
    stat = -2.0 * sum(math.log(p) for p in adj)
    k = len(adj)

    return float(chi2.sf(stat, 2 * k))


def find_multibase_positions(
    *,  # Enforce keyword arguments
    ref_base: str,
    filtered_read_dict: dict,
    base_cutoff: int,
    base_fraction_cutoff: float,
    base_support=None
) -> Tuple[dict, dict, int]:
    """
    Determine the characterised bases at the current position in the pileup,
    and whether any bases pass the specified cutoffs for SNV calling.

    Args:
        ref_base: String of the sequence of the reference gene at the current
        position
        filtered_read_dict: Dictionary with of base types: read direction:
        count
        base_cutoff: Integer of the number of identical mismatches in a column
        required for a SNV call
        base_fraction_cutoff: Float fraction of bases necessary to support a
        SNV call
        base_support: Optional dictionary of per-base support information

    Returns:
        snv_dict: Dictionary of characterised bases of the current column e.g.
        total, number matching reference
        sequence, number of SNVs in both forward and reverse reads, etc.
        passing_snv_dict: Dictionary summarising counts of categories of bases
        e.g. congruent, forward, reverse,
        paired
        total_coverage: Integer of the total number of bases passing filter in
        the pileup
    """
    # Initialise a dictionary to store the counts of the characterised bases
    snv_dict = {
        'total': 0,
        'total_congruent': 0,
        'total_congruent_SNV': 0,
        'total_forward': 0,
        'total_forward_SNV': 0,
        'total_reverse': 0,
        'total_reverse_SNV': 0,
        'total_SNV': 0
    }
    # Initialise the total depth of the column to zero
    total_coverage = 0

    # Initialise a dictionary to store the count for each individual
    # nucleotide at this position e.g. A:24, G:2
    base_count = {}

    # Iterate through the categories e.g. congruent_ref in the dictionary
    for category, base_dict in filtered_read_dict.items():
        # Iterate through the sequence of each query base, and the
        # corresponding count
        for base, count in base_dict.items():
            # Do not process the base if it has been flagged as 'filtered'
            if 'filtered' not in category:
                # Populate the base counting dictionary with the base sequence
                # and increment the count
                if base not in base_count:
                    base_count[base] = count
                else:
                    base_count[base] += count

                # Update the total coverage
                total_coverage += count
                snv_dict['total'] += count

            # Forward and reverse reads agree on base
            if 'congruent' in category:
                # Congruent reference sequence
                snv_dict['total_congruent'] += count
                snv_dict['total_forward'] += int(count / 2)
                snv_dict['total_reverse'] += int(count / 2)

                # Congruent SNVs
                if 'SNV' in category and base != ref_base:
                    snv_dict['total_congruent_SNV'] += count
                    snv_dict['total_forward_SNV'] += int(count / 2)
                    snv_dict['total_reverse_SNV'] += int(count / 2)
                    snv_dict['total_SNV'] += count

            # SNV in forward read
            elif category.startswith('forward_SNV'):
                snv_dict['total_forward'] += count
                if base != ref_base:
                    snv_dict['total_forward_SNV'] += count
                    snv_dict['total_SNV'] += count

            # Forward read matches reference
            elif category.startswith('forward_ref'):
                snv_dict['total_forward'] += count

            # SNV in reverse read
            elif category.startswith('reverse_SNV'):
                snv_dict['total_reverse'] += count
                if base != ref_base:
                    snv_dict['total_reverse_SNV'] += count
                    snv_dict['total_SNV'] += count

            # Reverse read match reference
            elif category.startswith('reverse_ref'):
                snv_dict['total_reverse'] += count

    # Initialise a dictionary to store the summary of characterised base types
    passing_snv_dict = {
        'congruent': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
        'forward': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
        'reverse': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
        'paired': {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    }

    # Boolean of whether there are bases passing filter, and the
    # passing_snv_dict should be used
    return_dict = False

    # Iterate through the categories in the dictionary
    for category, base_dict in filtered_read_dict.items():
        # Iterate through each query base and its corresponding count
        for base, count in base_dict.items():
            # Congruent SNVs
            if 'congruent' in category and base != ref_base:
                # If the base_cutoff_fraction has been provided, use it in
                # making SNV calls
                if base_fraction_cutoff:
                    # Ensure that the number of SNVs in the pileup is greater
                    # than the base_cutoff and the fraction of SNVs in the
                    # pileup is greater than the base_cutoff_fraction
                    if (
                        float(
                            base_count[base] / snv_dict['total']
                        ) >= base_fraction_cutoff and base_count[base]
                        >= base_cutoff
                    ):
                        # Update the summary dictionary with the base sequence
                        # and count
                        if base not in passing_snv_dict['forward']:
                            passing_snv_dict['congruent'][base] = count
                        else:
                            passing_snv_dict['congruent'][base] += count
                        return_dict = True

                # If base_cutoff_fraction is not supplied, only the number of
                # SNVs in the pileup must be greater than the base_cutoff
                # value in order for a SNV call
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['congruent']:
                            passing_snv_dict['congruent'][base] = count
                        else:
                            passing_snv_dict['congruent'][base] += count
                    return_dict = True

            if category.startswith('forward_SNV') and base != ref_base:
                if base_fraction_cutoff:
                    if (
                        float(
                            base_count[base] / snv_dict['total']
                        ) >= base_fraction_cutoff
                        and base_count[base] >= base_cutoff
                    ):
                        if base not in passing_snv_dict['forward']:
                            passing_snv_dict['forward'][base] = count
                        else:
                            passing_snv_dict['forward'][base] += count
                        return_dict = True
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['forward']:
                            passing_snv_dict['forward'][base] = count
                        else:
                            passing_snv_dict['forward'][base] += count
                        return_dict = True

            # SNVs in the reverse read
            if category.startswith('reverse_SNV') and base != ref_base:
                if base_fraction_cutoff:
                    if (
                        float(
                            base_count[base] / snv_dict['total']
                        ) >= base_fraction_cutoff
                        and base_count[base] >= base_cutoff
                    ):
                        if base not in passing_snv_dict['reverse']:
                            passing_snv_dict['reverse'][base] = count
                        else:
                            passing_snv_dict['reverse'][base] += count
                        return_dict = True
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['reverse']:
                            passing_snv_dict['reverse'][base] = count
                        else:
                            passing_snv_dict['reverse'][base] += count
                        return_dict = True

            # All unfiltered bases
            if base != ref_base and 'filtered' not in category:
                if base_fraction_cutoff:
                    if (
                        float(
                            base_count[base] / snv_dict['total']
                        ) >= base_fraction_cutoff
                        and base_count[base] >= base_cutoff
                    ):
                        if base not in passing_snv_dict['paired']:
                            passing_snv_dict['paired'][base] = count
                        else:
                            passing_snv_dict['paired'][base] += count
                        return_dict = True
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['paired']:
                            passing_snv_dict['paired'][base] = count
                        else:
                            passing_snv_dict['paired'][base] += count
                        return_dict = True

    # Compute per-base statistical tests if base_support was provided
    position_stats = {}

    if base_support:
        # flatten all qualities to make p_list for poisson-binomial
        # (prob of any error -> specific alt base ~ p/3)
        all_quals = []
        for b in base_support:
            all_quals.extend(base_support[b].get('quals', []))

        # Calculate p-values for each base
        p_list = [
            10 ** (-q / 10.0) / 3.0 for q in all_quals
        ] if all_quals else []

        # Iterate through each base in the base_count dictionary
        for base in base_count:
            if base == ref_base:
                continue
            k = base_count.get(base, 0)
            # p-value from poisson-binomial (exact) and mean-binomial
            # (conservative) - take the more conservative (larger) p
            if p_list:
                p_exact = poisson_binomial_tail(
                    k=k,
                    p_list=p_list
                )

                # Calculate mean probability for mean-binomial
                mean_p = sum(p_list) / len(p_list) if len(p_list) > 0 else 0.0
                p_binom = poisson_binomial_tail(
                    k=k,
                    p_list=[mean_p] * len(p_list)
                )

                # try beta-binomial if overdispersion present
                a, b = _estimate_beta_params(p_list=p_list)

                # Combine p-values
                if a is not None and b is not None:
                    p_bb = _beta_binomial_tail(
                        k=k,
                        n=len(p_list),
                        a=a,
                        b=b
                    )
                    p_value = max(p_exact, p_binom, p_bb)
                else:
                    p_value = max(p_exact, p_binom)
            else:
                p_value = None

            # strand bias
            fcount = base_support.get(base, {}).get('forward', 0)
            rcount = base_support.get(base, {}).get('reverse', 0)
            total_forward = snv_dict.get('total_forward', 0)
            total_reverse = snv_dict.get('total_reverse', 0)

            # Run Fisher's exact test
            strand_p = fisher_two_sided_p(
                a=fcount,
                b=total_forward - fcount,
                c=rcount,
                d=total_reverse - rcount
            )

            # Position bias
            base_pos = base_support.get(base, {}).get('positions', [])
            ref_pos = base_support.get(ref_base, {}).get('positions', [])

            # Mann-Whitney U test
            pos_p = mann_whitney_u_p(
                x=base_pos,
                y=ref_pos
            ) if base_pos and ref_pos else None

            # Mean qualities and mapping qualities
            mean_q_val = mean(
                base_support.get(
                    base, {}
                ).get(
                    'quals', []
                )
            ) if base_support.get(base, {}).get('quals') else None

            # Mean mapping quality
            mapqs = base_support.get(base, {}).get('mapqs', [])
            mean_mapq = mean(mapqs) if mapqs else None

            # Store the statistics for this base
            position_stats[base] = {
                'p_value': p_value,
                'strand_p': strand_p,
                'pos_p': pos_p,
                'mean_qual': mean_q_val,
                'mean_mapq': mean_mapq
            }
    if return_dict:
        return snv_dict, passing_snv_dict, total_coverage, position_stats

    return snv_dict, {}, total_coverage, position_stats


def position_details(
    *,  # Enforce keyword arguments
    actual_position: int,
    passing_snv_dict: Dict[str, int],
    contig_name: str,
    ref_base: str,
    total_coverage: int,
    base_cutoff: int,
    error_perc: float,
    p_value: Optional[float] = None,
    adj_p_value: Optional[float] = None,
    strand_p: Optional[float] = None,
    pos_p: Optional[float] = None,
    mean_q: Optional[float] = None,
    mean_mapq: Optional[float] = None
) -> Dict[str, Any]:
    """
    Format per-position statistics for reporting and return a dict
    representing the row.

    Args:
        actual_position: 1-based position in reference sequence.
        passing_snv_dict: Mapping base->coverage for bases that passed filters.
        contig_name: Contig name.
        ref_base: Reference base at this position.
        total_coverage: Total depth of coverage at the position.
        base_cutoff: Base cutoff used for calling.
        error_perc: Error percentage used in calculations.
        p_value, adj_p_value, strand_p, pos_p: Optional p-values calculated
        elsewhere.
        mean_q, mean_mapq: Optional mean base and mapping qualities.

    Returns:
        A dictionary with keys matching the TSV columns used by the pipeline.
    """
    # List of the base categories present in passing_snv_dict
    read_types = ['congruent', 'paired', 'forward', 'reverse']

    # Update the to_write string with basic information
    to_write = f'{contig_name}\t{actual_position}\t{ref_base}\t'

    # Initialise the coverage value to zero
    snv_coverage = 0

    # Iterate through all the read types
    for read_type in read_types:
        # Boolean of whether a semi-colon needs to be added to the string to
        # separate all the base sequences
        semi_colon = False
        for base, coverage in passing_snv_dict[read_type].items():
            # Ensure that the base isn't empty
            if base:
                if semi_colon:
                    to_write += ';'
                to_write += f'{base}:{coverage}'
                semi_colon = True

                # Increment the coverage only if the read type isn't 'paired'
                if read_type != 'paired':
                    snv_coverage += coverage
        to_write += ','

    # Format the error_perc to be more readable
    error_perc = f'{error_perc:0.2f}' if error_perc else 'ND'

    # Format the statistical fields
    p_str = f'{p_value:0.3e}' if p_value is not None else 'ND'
    q_str = f'{adj_p_value:0.3e}' if adj_p_value is not None else 'ND'
    strand_str = f'{strand_p:0.3e}' if strand_p is not None else 'ND'
    pos_str = f'{pos_p:0.3e}' if pos_p is not None else 'ND'
    mean_q_str = f'{mean_q:0.2f}' if mean_q is not None else 'ND'
    mean_mapq_str = f'{mean_mapq:0.2f}' if mean_mapq is not None else 'ND'

    # Finalise the to_write string with the remaining statistics
    to_write += (
        f'{snv_coverage}\t{total_coverage}\t{base_cutoff}\t{error_perc}\t'
        f'{p_str}\t{q_str}\t{strand_str}\t{pos_str}\t{mean_q_str}\t'
        f'{mean_mapq_str}\n'
    )
    return to_write


def read_contig(
    *,  # Enforce keyword arguments
    contig_name: str,
    bamfile_name: str,
    reference_fasta: str,
    allele_records: Any,
    fastq_records: Dict[str, Any],
    quality_cutoff: int = 20,
    base_cutoff: Optional[int] = None,
    base_fraction_cutoff: Optional[float] = None,
    fasta: bool = False,
    error_cutoff: float = 1.0,
    nanopore: bool = False
) -> Tuple[Dict[str, Any], str]:
    """
    Analyse a contig in a BAM file and return multibase dict and TSV text.

    Args:
        contig_name: Name of contig to analyse.
        bamfile_name: Path to BAM file.
        reference_fasta: Path to reference FASTA.
        allele_records: SeqIO index or similar mapping of alleles.
        fastq_records: Mapping of read_name -> SeqRecord with quality info.

    Keyword Args:
        quality_cutoff: Minimum base quality to be considered (default: 20).
        base_cutoff: Absolute base-count cutoff (optional).
        base_fraction_cutoff: Fractional cutoff of coverage (optional).
        fasta: If True, operate in FASTA mode (no base qualities).
        error_cutoff: Error cutoff used when computing base_cutoff adjustments.
        nanopore: If True, enable nanopore-specific handling.

    Returns:
        A tuple (multibase_dict, tsv_output) where multibase_dict is a dict of
        detected multibase positions and gene summary stats, and tsv_output is
        the position-level TSV text for the contig.
    """
    # Initialize pysam FastaFile object and other variables
    pysam_fasta = pysam.FastaFile(reference_fasta)
    multibase_position_dict = {}
    to_write = str()

    # If analysing FASTA files, a single base difference is all that
    # is expected
    if fasta:
        base_cutoff = 1

    # Extract the reference sequence for the contig being analysed
    reference_sequence = str(allele_records[contig_name].seq)

    # Parse the BAM file with pysam to create AlignmentFile, and
    # AlignmentFile.pileup objects
    bamfile, pileup = parse_bam(
        bamfile_name=bamfile_name,
        contig_name=contig_name,
        pysam_fasta=pysam_fasta
    )

    # Define dictionaries and lists to store filtered reads, base support and
    # quality scores
    filtered_read_dict = {}
    base_support_dict = {}
    quality_list = []

    # Containers for report aggregation and multiple-testing correction
    report_entries = []
    all_tests = []

    # Iterate through each column in the pileup
    for i, column in enumerate(pileup):
        filtered_reads, qualities, base_support = characterise_read(
            column=column,
            reference_sequence=reference_sequence,
            fastq_records=fastq_records,
            quality_cutoff=quality_cutoff,
            fasta=fasta,
            nanopore=nanopore
        )

        # Populate the dictionaries and lists with the results from
        # characterise_read
        filtered_read_dict[i] = filtered_reads
        base_support_dict[i] = base_support
        quality_list += qualities

    # Initialise the calculated error percentage to zero
    error_perc = None

    # If the base_cutoff is set to zero, determine the appropriate cutoff value
    if base_cutoff == 0:
        try:
            computed = determine_cutoff(
                qualities=quality_list,
                reference_sequence=reference_sequence,
                base_cutoff=base_cutoff,
                error_cutoff=error_cutoff
            )
            # defensive: determine_cutoff historically returns an int
            if isinstance(computed, tuple):
                computed_cutoff = int(computed[0])
                if len(computed) > 1:
                    error_perc = computed[1]
            else:
                computed_cutoff = int(computed)
            base_cutoff = computed_cutoff
            logging.debug(
                'Contig %s: computed dynamic base_cutoff=%s, error_perc=%s',
                contig_name, base_cutoff, error_perc
            )
        except Exception:
            logging.debug(
                'Contig %s: error computing dynamic base cutoff: %s',
                contig_name, traceback.format_exc()
            )
    bamfile.close()

    # It seems that the pileup (generator?) is used up above, so it must be
    # recreated
    bamfile, pileup = parse_bam(
        bamfile_name=bamfile_name,
        contig_name=contig_name,
        pysam_fasta=pysam_fasta
    )

    # Iterate through each column in the pileup
    for i, column in enumerate(pileup):
        # Extract the sequence of the reference gene at the current position
        ref_base = reference_sequence[column.pos]

        # Summarise the pileup (now collecting statistics per position)
        _, passing_snv_dict, total_coverage, position_stats = \
            find_multibase_positions(
                ref_base=ref_base,
                filtered_read_dict=filtered_read_dict[i],
                base_cutoff=base_cutoff,
                base_fraction_cutoff=base_fraction_cutoff,
                base_support=base_support_dict.get(i, {})
            )

        # If there are any SNVs called for the gene, update the
        # multibase_position_dict and aggregate stats
        if passing_snv_dict:
            # Pysam starts counting at 0, whereas we actually want to start
            # counting at 1.
            actual_position = column.pos + 1

            # Initialise the gene name in the dictionary as required
            if column.reference_name not in multibase_position_dict:
                multibase_position_dict[column.reference_name] = {}

            # Update the dictionary with the actual position:
            multibase_position_dict[column.reference_name].update(
                {actual_position: passing_snv_dict}
            )

            # Store entry for per-gene reporting and FDR correction
            report_entries.append({
                'gene': column.reference_name,
                'position': actual_position,
                'ref_base': ref_base,
                'passing_snv_dict': passing_snv_dict,
                'total_coverage': total_coverage,
                'base_cutoff': base_cutoff,
                'error_perc': error_perc,
                'position_stats': position_stats
            })

            # Collect p-values for BH adjustment (per-base tests)
            for base, stats in (position_stats or {}).items():
                pval = stats.get('p_value')
                if pval is not None:
                    all_tests.append(
                        {
                            'entry_idx': len(report_entries) - 1,
                            'base': base, 'p': pval
                        }
                    )

    # Multiple-testing correction across bases in this gene
    if all_tests:
        pvals = [t['p'] for t in all_tests]
        qvals = benjamini_hochberg(pvals=pvals)
        for idx, t in enumerate(all_tests):
            t['q'] = qvals[idx]
            entry = report_entries[t['entry_idx']]
            base = t['base']

            # Attach test summary to entry
            if 'tests' not in entry:
                entry['tests'] = []
            entry['tests'].append({'base': base, 'p': t['p'], 'q': t['q']})

        # Annotate q-values back into position_stats
        for entry in report_entries:
            pos_stats = entry.get('position_stats') or {}
            for test in entry.get('tests', []):
                b = test['base']
                if b in pos_stats:
                    pos_stats[b]['q_value'] = test['q']

    # Compute per-gene combined p-value (Fisher) and a gene-level score
    # (sum -log10(q)) p-values used for Fisher should be the unadjusted
    # position p-values; q-values are used for scoring.
    pvals_for_gene = [t['p'] for t in all_tests] if all_tests else []
    qvals_for_gene = [
        t.get('q') for t in all_tests if 'q' in t
    ] if all_tests else []

    # Calculate combined p-value and gene score
    combined_p = combine_pvalues_fisher(
        pvals=pvals_for_gene
    ) if pvals_for_gene else None

    # Initialize gene_score and num_sig_positions
    gene_score = None

    # Calculate gene_score and num_sig_positions if q-values are available
    if qvals_for_gene:
        gene_score = sum([-math.log10(max(q, 1e-300)) for q in qvals_for_gene])

    # Calculate number of significant positions based on q-value threshold
    num_sig_positions = sum(
        1 for q in qvals_for_gene if q is not None and q <= 0.05
    ) if qvals_for_gene else 0

    # Store gene-level stats in the multibase_position_dict under a
    # special key for this gene
    if contig_name not in multibase_position_dict:
        multibase_position_dict[contig_name] = {}
    multibase_position_dict[contig_name]['_gene_stats'] = {
        'combined_p': combined_p,
        'gene_score': gene_score,
        'num_sig_positions': num_sig_positions
    }

    # Build a per-gene summary TSV file for this sample. Aggregate across all
    # per-gene dicts returned later. (Actual file will be written once in
    # find_contamination after collecting all genes.) build output lines
    to_write = ''
    for entry in report_entries:
        pos_stats = entry.get('position_stats') or {}
        if pos_stats:
            # Choose the base with smallest p-value for summary metrics
            min_base = min(
                pos_stats.items(), key=lambda kv: kv[1].get('p_value', 1.0)
            )[0]
            stats = pos_stats[min_base]
            pval = stats.get('p_value')
            qval = stats.get('q_value')
            strand_p = stats.get('strand_p')
            pos_p = stats.get('pos_p')
            mean_q = stats.get('mean_qual')
            mean_mapq = stats.get('mean_mapq')
        else:
            pval = qval = strand_p = pos_p = mean_q = mean_mapq = None
        to_write += position_details(
            actual_position=entry['position'],
            passing_snv_dict=entry['passing_snv_dict'],
            contig_name=entry['gene'],
            ref_base=entry['ref_base'],
            total_coverage=entry['total_coverage'],
            base_cutoff=entry['base_cutoff'],
            error_perc=entry['error_perc'],
            p_value=pval,
            adj_p_value=qval,
            strand_p=strand_p,
            pos_p=pos_p,
            mean_q=mean_q,
            mean_mapq=mean_mapq
        )
    bamfile.close()

    return multibase_position_dict, to_write


def _read_contig_dispatch(kwargs: Dict[str, Any]) -> Tuple[Dict[str, Any], str]:
    """Helper for multiprocessing that calls `read_contig` with kwargs.

    ``multiprocessing.Pool.map`` passes a single argument to the worker
    function. Since ``read_contig`` enforces keyword-only arguments, we
    build a dict of kwargs and dispatch via this helper.
    """
    return read_contig(**kwargs)


def count_multibase_positions(
    *,  # Enforce keyword arguments
    multibase_dict_list: List[Dict[str, Any]]
) -> int:
    """
    Count multibase positions across the list of per-gene multibase dicts.

    The per-gene dicts have the shape {gene: {position: {...}, '_gene_stats':
    {...}}}.

    This function excludes meta keys that start with '_' (e.g. '_gene_stats').

    Args:
        multibase_dict_list: List of dicts returned by `read_contig` for each
        gene.

    Returns:
        Total number of positions (int).
    """
    # Initialize total to zero
    total = 0

    # Iterate through each gene's multibase dict
    for multibase_position_dict in multibase_dict_list:
        # Iterate through each gene and its SNP positions
        for _, snp_positions in multibase_position_dict.items():
            # Increment total by the number of positions excluding meta keys
            total += sum(
                1 for k in snp_positions.keys() if not str(k).startswith('_')
            )

    return total


def find_rmlst_type(
    *,  # Enforce keyword arguments
    kma_report: str, rmlst_report: str
) -> List[str]:
    """
    Uses a report generated by KMA to determine what allele is present for
    each rMLST gene.

    Args:
        kma_report: Path to KMA report file.
        rmlst_report: Path to output rMLST report file.

    Returns:
        List of strings representing the gene and allele called e.g.
        ['abcZ_1', 'adk_4', ...]
    """
    # Initialize a dictionary to store the best scoring allele for each gene
    genes_to_use = {}

    # Initialize a dictionary to store the highest score for each gene
    score_dict = {}

    # Initialize a list to store the final gene_allele strings
    gene_alleles = []

    # Read through the KMA report file
    with open(kma_report, encoding='utf-8') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        # Iterate through each row in the KMA report
        for row in reader:
            # Extract the gene_allele and score from the row
            gene_allele = row['#Template']

            # Extract the score from the row
            score = int(row['Score'])

            # Split the gene_allele into gene and allele components
            gene = gene_allele.split('_')[0]
            allele = gene_allele.split('_')[1]

            # Determine if this allele has the highest score for the gene
            if gene not in score_dict:
                score_dict[gene] = score
                genes_to_use[gene] = allele
            else:
                if score > score_dict[gene]:
                    score_dict[gene] = score
                    genes_to_use[gene] = allele

    # Create the gene_allele strings and write to the rMLST report file
    for gene, allele in genes_to_use.items():
        gene_alleles.append(gene + '_' + allele.replace(' ', ''))

    # Sort the gene_alleles list before writing to file
    gene_alleles = sorted(gene_alleles)

    # Write the rMLST report file
    with open(rmlst_report, 'w', encoding='utf-8') as f:
        # rMLST report is a TSV: Gene\tAllele
        f.write('Gene\tAllele\n')

        # Write each gene and allele to the report file
        for gene_allele in gene_alleles:
            gene = gene_allele.split('_')[0]
            allele = gene_allele.split('_')[1]
            f.write(f'{gene}\t{allele}\n')
    return gene_alleles


def base_dict_to_string(
    *,  # Enforce keyword arguments
    base_dict: Dict[str, int]
) -> str:
    """
    Converts a dictionary to a string. {'C': 12, 'A':4} gets converted to
    C:12;A:4

    Args:
        base_dict: Dictionary of bases and counts created by find_if_multibase

    Returns:
        A string representation of the base_dict
    """
    # Initialize the output string
    outstr = ''
    # First, sort base_dict so that major allele always comes first -
    # makes output report nicer to look at.
    base_list = sorted(base_dict.items(), key=lambda kv: kv[1], reverse=True)
    for base in base_list:
        outstr += f'{base[0]}:{base[1]};'

    return outstr[:-1]


def find_total_sequence_length(
    *,  # Enforce keyword arguments
    fasta_file: str
) -> int:
    """
    Totals up number of bases in a fasta file.

    Args:
        fasta_file: Path to FASTA file.

    Returns:
        Total number of bases in the FASTA file.
    """
    # Initialize total_length to zero
    total_length = 0

    # Iterate through each sequence in the FASTA file and sum lengths
    for sequence in SeqIO.parse(fasta_file, 'fasta'):
        total_length += len(sequence.seq)

    return total_length


def load_fastq_records(
    *,  # Enforce keyword arguments
    gz: str,
    paired: bool,
    forward: bool
) -> Dict[str, Any]:
    """
    Use SeqIO to load FASTQ records from file

    Args:
        gz: Path to FASTQ file (can be gzipped).
        paired: Boolean of whether reads are paired.
        forward: Boolean of whether reads are forward reads.

    Returns:
        Dictionary of SeqIO records keyed by read ID.
    """
    # Initialise a dictionary to store the FASTQ records
    records = {}
    # Iterate through the reads
    for record in SeqIO.parse(gz, 'fastq'):
        # Only update the naming scheme for paired reads
        if paired:
            if forward:
                # Don't worry if the record.id already has a /1
                if not record.id.endswith('/1'):
                    record.id = record.id + '/1'
            # Process reverse reads in a similar fashion to forward reads
            else:
                if not record.id.endswith('/2'):
                    record.id = record.id + '/2'
        else:
            pass
        records.update(SeqIO.to_dict([record]))
    return records


def index_databases(
    *,  # Enforce keyword arguments
    sample_database: str
) -> None:
    """
    Index the database file with pysam and kma

    Args:
        sample_database: Path to FASTA database file.
    """
    # Don't bother re-indexing, this only needs to happen once.
    if not os.path.isfile(sample_database + '.fai'):
        try:
            pysam.faidx(sample_database)
        except pysam.utils.SamtoolsError:
            pass

    # Set the KMA database name
    kma_database = sample_database.replace('.fasta', '') + '_kma'

    # The .name is one of the files KMA creates when making a database.
    if not os.path.isfile(kma_database + '.name'):
        logging.info(
            'Since this is the first time you are using this database, it '
            'needs to be indexed by KMA. This might take a while'
        )

        # Set the KMA indexing command
        cmd = f'kma index -i {sample_database} -o {kma_database}'

        # Set output variable
        out = str()

        # Run the KMA indexing command
        try:
            out, err = run_cmd(cmd=cmd)
        except subprocess.CalledProcessError as exc:
            err = exc

        # Write to logfile
        log = sample_database + '_log.txt'
        write_to_logfile(
            logfile=log,
            out=out,
            err=err,
            cmd=cmd
        )


def find_contamination(
    *,  # Enforce keyword arguments
    pair: List[str],
    output_folder: str,
    databases_folder: str,
    base_cutoff: int,
    forward_id: str = '_R1',
    threads: int = 1,
    keep_files: bool = False,
    quality_cutoff: int = 20,
    base_fraction_cutoff: float = 0.05,
    cgmlst_db: Optional[str] = None,
    xmx: Optional[str] = None,
    tmpdir: Optional[str] = None,
    data_type: str = 'Illumina',
    use_rmlst: bool = False,
    min_matching_hashes: int = 40,
    fasta: bool = False,
    error_cutoff: float = 1.0,
    debug: bool = False,
    use_prob_scoring: bool = False,
    score_threshold: float = 10.0
) -> Optional[Tuple[str, bool]]:
    """
    Run contamination detection for a sample (paired or single reads).

    Args:
        pair: Pair or single read list for the sample (e.g. [R1, R2] or
        [single]).
        output_folder: Output folder for results.
        databases_folder: Path to ConFindr databases.
        base_cutoff: Number of bases required to support a variant call.
        forward_id: Identifier for forward reads.
        threads: Number of threads to use.
        keep_files: Keep temporary files if True.
        quality_cutoff: Base quality cutoff (default 20).
        base_fraction_cutoff: Fractional cutoff for variant support.
        cgmlst_db: Optional cgMLST DB name.
        xmx: Optional memory string for BBMap tools.
        tmpdir: Temporary directory to use.
        data_type: 'Illumina' or 'Nanopore'.
        use_rmlst: Prefer rMLST DB when True.
        min_matching_hashes: Mash matching threshold.
        fasta: If True, operate in FASTA-only mode.
        error_cutoff: Error cutoff percentage (default 1.0).
        debug: Enable debug logging.
        use_prob_scoring: Use probabilistic scoring (sum -log10(q)).
        score_threshold: Threshold for probabilistic scoring.

    Returns:
        Optional tuple (sample_report_tsv, contamination_boolean), or None on
        failure.
    """
    # Set the name and path for the database download date file
    download_date_file = os.path.join(databases_folder, 'download_date.txt')
    if os.path.isfile(download_date_file):
        with open(download_date_file, 'r', encoding='utf-8') as f:
            database_download_date = f.readline().rstrip()
    else:
        database_download_date = 'ND'

    # Define the log file for this sample
    log = os.path.join(output_folder, 'confindr_log.txt')
    if len(pair) == 2:
        sample_name = os.path.split(pair[0])[-1].split(forward_id)[0]
        paired = True
        logging.debug('Sample is paired. Sample name is %s', sample_name)
    else:
        sample_name = os.path.split(pair[0])[-1].split('.')[0]
        paired = False
        logging.debug('Sample is unpaired. Sample name is %s', sample_name)
    sample_tmp_dir = os.path.join(output_folder, sample_name)
    if not os.path.isdir(sample_tmp_dir):
        os.makedirs(sample_tmp_dir)

    logging.info('Checking for cross-species contamination...')
    if paired:
        genus = find_cross_contamination(
            databases=databases_folder,
            reads=pair,
            sample_name=sample_name,
            tmpdir=sample_tmp_dir,
            log=log,
            threads=threads,
            min_matching_hashes=min_matching_hashes
        )
    else:
        genus = find_cross_contamination(
            databases=databases_folder,
            reads=pair[0],
            sample_name=sample_name,
            tmpdir=sample_tmp_dir,
            log=log,
            threads=threads,
            min_matching_hashes=min_matching_hashes
        )

    # Setup genus-specific databases, if necessary.
    if cgmlst_db is not None:
        # Sanity check that the DB specified is actually a file, otherwise,
        # quit with appropriate error message.
        if not os.path.isfile(cgmlst_db):
            logging.error(
                'ERROR: Specified cgMLST file (%s) does not exist. Please '
                'check the path and try again.', cgmlst_db
            )
            sys.exit(1)
        sample_database = cgmlst_db
    else:
        db_folder = databases_folder
        if genus != 'ND':
            # Logic here is as follows: users can either have both rMLST
            # databases, which cover all of bacteria, cgmlst-derived databases
            # which cover only Escherichia, Salmonella, and Listeria
            # (may add more at some point), or they can have both. They can
            # also set priority to either always use rMLST, or to use my
            # core-genome derived stuff and fall back on rMLST if they're
            # trying to look at a genus I haven't created a scheme for.
            if len(genus.split(':')) > 1:
                predominant_genus = genus.split(':')[0]
            else:
                predominant_genus = genus
            # In the event rmlst databases have priority, always use them.
            if use_rmlst is True:
                sample_database = os.path.join(
                    db_folder,
                    f'{predominant_genus}_db.fasta'
                )

                # Create genus specific database if it doesn't already exist
                # and we have the necessary rMLST files.
                if not os.path.isfile(sample_database):
                    if os.path.isfile(
                        os.path.join(
                            db_folder,
                            'gene_allele.txt'
                            )
                    ) and os.path.isfile(
                        os.path.join(
                            db_folder,
                            'rMLST_combined.fasta'
                        )
                    ):
                        logging.info(
                            'Setting up rMLST genus-specific database for '
                            'genus %s...', predominant_genus
                        )

                        # Calculate the list of alleles
                        allele_list = find_genus_specific_allele_list(
                            profiles_file=os.path.join(
                                db_folder,
                                'gene_allele.txt'
                            ),
                            target_genus=predominant_genus
                        )

                        # Create the database in a temporary directory if
                        # specified
                        if tmpdir:
                            logging.info(
                                'Using temporary directory for database '
                                'creation: %s', tmpdir
                            )

                            # Create the temporary directory if necessary
                            os.makedirs(tmpdir, exist_ok=True)

                            # Set the sample database path
                            sample_database = os.path.join(
                                tmpdir,
                                f'{predominant_genus}_db.fasta'
                            )

                            # Create the allele-specific database
                            setup_allelespecific_database(
                                fasta_file=sample_database,
                                database_folder=db_folder,
                                allele_list=allele_list
                            )
                        else:
                            setup_allelespecific_database(
                                fasta_file=sample_database,
                                database_folder=db_folder,
                                allele_list=allele_list
                            )
            else:
                # Check if a cgderived database is available. If not, try to
                # use rMLST database.
                sample_database = os.path.join(
                    db_folder,
                    f'{predominant_genus}_db_cgderived.fasta'
                )

                # If cgderived database doesn't exist, fall back on rMLST db.
                if not os.path.isfile(sample_database):
                    sample_database = os.path.join(
                        db_folder,
                        f'{predominant_genus}_db.fasta'
                    )

                    # Create genus specific database if it doesn't already
                    # exist and we have the necessary rMLST files.
                    if os.path.isfile(
                        os.path.join(
                            db_folder,
                            'rMLST_combined.fasta'
                        )
                    ) and os.path.isfile(
                        os.path.join(
                            db_folder,
                            'gene_allele.txt'
                        )
                    ) and not os.path.isfile(
                        sample_database
                    ):
                        logging.info(
                            'Setting up core genome genus-specific database '
                            'for genus %s...', predominant_genus
                        )

                        # Calculate the list of alleles
                        allele_list = find_genus_specific_allele_list(
                            profiles_file=os.path.join(
                                db_folder,
                                'gene_allele.txt'
                            ),
                            target_genus=predominant_genus
                        )

                        # Create the database in a temporary directory if
                        # specified
                        if tmpdir:
                            logging.info(
                                'Using temporary directory for database '
                                'creation: %s', tmpdir
                            )
                            # Create the temporary directory if necessary
                            os.makedirs(tmpdir, exist_ok=True)

                            # Set the sample database path
                            sample_database = os.path.join(
                                tmpdir,
                                f'{predominant_genus}_db.fasta'
                            )

                            # Create the allele-specific database
                            setup_allelespecific_database(
                                fasta_file=sample_database,
                                database_folder=db_folder,
                                allele_list=allele_list
                            )
                        else:
                            setup_allelespecific_database(
                                fasta_file=sample_database,
                                database_folder=db_folder,
                                allele_list=allele_list
                            )

        else:
            sample_database = os.path.join(db_folder, 'rMLST_combined.fasta')

    # If a user has gotten to this point and they don't have any database
    # available to do analysis because they don't have rMLST downloaded and
    # we don't have a cg-derived database available, boot them with a helpful
    # message.
    if not os.path.isfile(sample_database):
        write_output(
            output_report=os.path.join(output_folder, 'confindr_report.tsv'),
            sample_name=sample_name,
            multi_positions=0,
            genus=genus,
            total_gene_length=0,
            database_download_date=database_download_date
        )
        logging.info(
            'Did not find databases for genus %s. You can download the rMLST '
            'database to get access to all genera (see https://'
            'olc-bioinformatics.github.io/ConFindr/install/). Alternatively, '
            'if you have a high-quality core-genome derived database for your '
            'genome of interest, we would be happy to add it - open an issue '
            'at https://github.com/OLC-Bioinformatics/ConFindr/issues with '
            'the title "Add genus-specific database: %s"\n', genus, genus
            )
        if keep_files is False:
            shutil.rmtree(sample_tmp_dir)
        return

    # Extract rMLST reads and quality trim.
    logging.info('Extracting conserved core genes...')
    out, err, cmd = '', '', ''
    forward_bait = str()
    forward_trimmed = str()
    reverse_bait = str()
    reverse_trimmed = str()
    unpaired_bait = str()
    unpaired_trimmed = str()

    # Bait reads matching the database
    if paired:
        forward_bait = os.path.join(
            sample_tmp_dir,
            f'{sample_name}_baited_R1.fastq.gz'
        )
        reverse_bait = forward_bait.replace('_R1', '_R2')

        # Only run if the baited files don't already exist
        if not os.path.isfile(forward_bait):
            if xmx is None:
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair[0],
                    reverse_in=pair[1],
                    forward_out=forward_bait,
                    reverse_out=reverse_bait,
                    threads=threads,
                    returncmd=True
                )
            else:
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair[0],
                    reverse_in=pair[1],
                    forward_out=forward_bait,
                    reverse_out=reverse_bait,
                    threads=threads,
                    Xmx=xmx,
                    returncmd=True
                )
    else:
        # Still name the file '_baited_trimmed' even if the file won't be
        # trimmed
        if data_type == 'Nanopore' or fasta:
            unpaired_bait = os.path.join(
                sample_tmp_dir,
                f'{sample_name}_baited_trimmed.fastq.gz'
            )
        else:
            unpaired_bait = os.path.join(
                sample_tmp_dir,
                f'{sample_name}_baited.fastq.gz'
            )

        # Only run if the baited file doesn't already exist
        if not os.path.isfile(unpaired_bait):
            if xmx is None:
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair[0],
                    forward_out=unpaired_bait,
                    returncmd=True, threads=threads
                )
            else:
                out, err, cmd = bbtools.bbduk_bait(
                    reference=sample_database,
                    forward_in=pair[0],
                    forward_out=unpaired_bait,
                    Xmx=xmx,
                    returncmd=True,
                    threads=threads
                )
    if out:
        write_to_logfile(
            logfile=log,
            out=out,
            err=err,
            cmd=cmd
        )

    # Verify baiting output exists before trimming. If the expected bait
    # files are missing, write a failure line to the summary report and
    # abort this sample early.
    if paired:
        if not (os.path.isfile(forward_bait) and os.path.isfile(reverse_bait)):
            logging.error(
                'Baiting failed to create files for %s: %s, %s',
                sample_name, forward_bait, reverse_bait
            )
            write_output(
                output_report=os.path.join(
                    output_folder, 'confindr_report.tsv'
                ),
                sample_name=sample_name,
                multi_positions=0,
                genus='Error processing sample',
                total_gene_length=0,
                database_download_date=database_download_date
            )
            if keep_files is False:
                try:
                    shutil.rmtree(sample_tmp_dir)
                except OSError:
                    pass
            return
    else:
        if not os.path.isfile(unpaired_bait):
            logging.error(
                'Baiting failed to create file for %s: %s',
                sample_name, unpaired_bait
            )
            write_output(
                output_report=os.path.join(
                    output_folder, 'confindr_report.tsv'
                ),
                sample_name=sample_name,
                multi_positions=0,
                genus='Error processing sample',
                total_gene_length=0,
                database_download_date=database_download_date
            )
            if keep_files is False:
                try:
                    shutil.rmtree(sample_tmp_dir)
                except OSError:
                    pass
            return

    # Run quality trimming on the baited reads
    logging.info('Quality trimming...')
    out, err, cmd = '', '', ''

    # Handle Illumina and Nanopore data differently
    if data_type == 'Illumina':
        if paired:
            forward_trimmed = os.path.join(
                sample_tmp_dir,
                f'{sample_name}_baited_trimmed_R1.fastq.gz'
            )

            # Set the reverse trimmed filename
            reverse_trimmed = forward_trimmed.replace('_R1', '_R2')

            # Only run trimming if the trimmed files don't already exist
            if not os.path.isfile(forward_trimmed):
                # Use the appropriate bbduk trimming command based on whether
                # xmx is set
                if xmx is None:
                    out, err, cmd = bbtools.bbduk_trim(
                        forward_in=forward_bait,
                        reverse_in=reverse_bait,
                        forward_out=forward_trimmed,
                        reverse_out=reverse_trimmed,
                        threads=str(threads), returncmd=True
                    )
                else:
                    out, err, cmd = bbtools.bbduk_trim(
                        forward_in=forward_bait,
                        reverse_in=reverse_bait,
                        forward_out=forward_trimmed,
                        reverse_out=reverse_trimmed,
                        Xmx=xmx,
                        threads=str(threads),
                        returncmd=True
                    )

            # Load the trimmed FASTQ records into a dictionary
            with gzip.open(forward_trimmed, 'rt') as gz:
                fastq_records = load_fastq_records(
                    gz=gz,
                    paired=True,
                    forward=True
                )
            with gzip.open(reverse_trimmed, 'rt') as gz:
                # fastq_records.update(SeqIO.to_dict(SeqIO.parse(gz, 'fastq')))
                fastq_records.update(
                    load_fastq_records(
                        gz=gz,
                        paired=True,
                        forward=False
                    )
                )
        # Handle unpaired reads
        else:
            unpaired_trimmed = os.path.join(
                sample_tmp_dir,
                f'{sample_name}_baited_trimmed.fastq.gz'
            )

            # Process FASTA mode differently - no trimming
            if not fasta:
                # Only run trimming if the trimmed file doesn't already exist
                if not os.path.isfile(unpaired_trimmed):
                    # Use the appropriate bbduk trimming command based on
                    # whether xmx is set
                    if xmx is None:
                        out, err, cmd = bbtools.bbduk_trim(
                            forward_in=unpaired_bait,
                            forward_out=unpaired_trimmed,
                            returncmd=True,
                            threads=threads
                        )
                    else:
                        out, err, cmd = bbtools.bbduk_trim(
                            forward_in=unpaired_bait,
                            forward_out=unpaired_trimmed,
                            returncmd=True,
                            threads=threads,
                            Xmx=xmx
                        )

                # Load the trimmed FASTQ records into a dictionary
                with gzip.open(unpaired_trimmed, 'rt') as gz:
                    # Load the FASTQ records
                    fastq_records = load_fastq_records(
                        gz=gz,
                        paired=False,
                        forward=True
                    )
            else:
                # Unpaired_bait
                with gzip.open(unpaired_bait, 'rt') as gz:
                    # Load the FASTQ records
                    fastq_records = load_fastq_records(
                        gz=gz,
                        paired=False,
                        forward=True
                    )
        write_to_logfile(
            logfile=log,
            out=out,
            err=err,
            cmd=cmd
        )
    # If Nanopore data, no trimming - just load the baited reads
    else:
        if paired:
            with gzip.open(forward_bait, 'rt') as gz:
                # Load the FASTQ records
                fastq_records = load_fastq_records(
                    gz=gz,
                    paired=True,
                    forward=True
                )
            with gzip.open(reverse_bait, 'rt') as gz:
                # Load the FASTQ records
                fastq_records.update(load_fastq_records(
                    gz=gz,
                    paired=True,
                    forward=False
                ))
        else:
            with gzip.open(unpaired_bait, 'rt') as gz:
                # Load the FASTQ records
                fastq_records = load_fastq_records(
                    gz=gz,
                    paired=False,
                    forward=True
                )
    logging.info('Detecting contamination...')

    # Now do mapping in two steps - first, map reads back to database with
    # ambiguous reads matching all - this will be used to get a count of
    # number of reads aligned to each gene/allele so we can create a custom
    # rMLST file with only the most likely allele for each gene.
    kma_report = os.path.join(
        sample_tmp_dir,
        f'{sample_name}_kma'
    )

    # Set the KMA database name
    kma_database = sample_database.replace('.fasta', '') + '_kma'

    # Index the database if necessary
    index_databases(sample_database=sample_database)

    # Run KMA.
    if paired:
        if not os.path.isfile(kma_report + '.res'):
            cmd = (
                f'kma -ipe {forward_trimmed} {reverse_trimmed} '
                f'-t_db {kma_database} -o {kma_report} -t {threads}'
            )
            out, err = run_cmd(cmd=cmd)

            # Write to logfile
            write_to_logfile(
                logfile=log,
                out=out,
                err=err,
                cmd=cmd
            )
    # Unpaired reads
    else:
        if not os.path.isfile(kma_report + '.res'):
            if data_type == 'Illumina':
                # Use the FASTA file (rather than the reads) as the input
                if fasta:
                    cmd = (
                        f'kma -i {pair[0]} -t_db {kma_database} -mem_mode '
                        f'-ID 100 -ConClave 2 -ex_mode -o {kma_report} '
                        f'-t {threads}'
                    )
                else:
                    cmd = (
                        f'kma -i {unpaired_trimmed} -t_db {kma_database} '
                        f'-o {kma_report} -t {threads}'
                    )
            else:
                # Recommended Nanopore settings from KMA repo:
                # https://bitbucket.org/genomicepidemiology/kma
                cmd = (
                    f'kma -i {unpaired_bait} -t_db {kma_database} '
                    f'-o {kma_report} -mem_mode -mp 20 -mrs 0.0 -bcNano '
                    f'-t {threads}'
                )
            out, err = run_cmd(cmd=cmd)
            write_to_logfile(
                logfile=log,
                out=out,
                err=err,
                cmd=cmd
            )

    # Set the rMLST report path
    rmlst_report = os.path.join(output_folder, sample_name + '_alleles.tsv')

    # Parse the KMA report to find the best allele for each rMLST gene
    gene_alleles = find_rmlst_type(
        kma_report=kma_report + '.res',
        rmlst_report=rmlst_report
    )

    # Set the rMLST FASTA path
    rmlst_fasta = os.path.join(
        sample_tmp_dir,
        f'{sample_name}_alleles.fasta'
    )

    # Check if the rMLST FASTA file already exists
    if not os.path.isfile(rmlst_fasta):
        # Create a FASTA file that has only one allele per rMLST gene
        with open(rmlst_fasta, 'w', encoding='utf-8') as f:
            # Use SeqIO to parse the database and write out only the
            # relevant alleles
            for contig in SeqIO.parse(sample_database, 'fasta'):
                if contig.id in gene_alleles:
                    SeqIO.write(contig, f, 'fasta')

    # Get total gene length for later reporting
    rmlst_gene_length = find_total_sequence_length(
        fasta_file=rmlst_fasta
    )

    logging.debug('Total gene length is %s', rmlst_gene_length)

    # Initialize contamination boolean
    pysam_pass = True
    # Second step of mapping - Do a mapping of our baited reads against a
    # fasta file that has only one allele per rMLST gene.
    try:
        # Index the rMLST FASTA file if necessary
        if not os.path.isfile(rmlst_fasta + '.fai'):
            pysam.faidx(rmlst_fasta)

        # Set the path to the contamination BAM, sorted BAM, and SAM files
        outbam = os.path.join(
            sample_tmp_dir,
            f'{sample_name}_contamination.bam'
        )
        sorted_bam = os.path.join(
            sample_tmp_dir,
            f'{sample_name}_contamination_sorted.bam'
        )

        # Run the mapping if the sorted BAM doesn't already exist
        if not os.path.isfile(sorted_bam):
            # Perform mapping differently for paired and unpaired reads
            if paired:
                cmd = (
                    f'bbmap.sh ref={rmlst_fasta} in={forward_trimmed} '
                    f'in2={reverse_trimmed} out={outbam} threads={threads} '
                    'mdtag nodisk'
                )
                # Add subfilter for cgMLST databases
                if cgmlst_db is not None:
                    # Lots of core genes seem to have relatives within a
                    # genome that are at ~70% identity. This means that reads
                    # that shouldn't map do, and cause false positives. Adding
                    # in this sub-filter means that reads can only have one
                    # mismatch, so they have to be from the right gene for
                    # this to work.
                    cmd += ' subfilter=1'
                # Add in memory string if specified
                if xmx:
                    cmd += f' -Xmx{xmx}'

                # Run the command
                out, err = run_cmd(cmd=cmd)
                write_to_logfile(
                    logfile=log,
                    out=out,
                    err=err,
                    cmd=cmd
                )
            else:
                # If Illumina FASTQ:
                if data_type == 'Illumina' and not fasta:
                    cmd = (
                        f'bbmap.sh ref={rmlst_fasta} in={unpaired_trimmed} '
                        f'out={outbam} threads={threads} mdtag nodisk'
                    )

                    # Add subfilter for cgMLST databases
                    if cgmlst_db is not None:
                        # Core genes can have relatives within a genome that
                        # are at ~70 percent identity. This means that reads
                        # that shouldn't map do, and cause false positives.
                        # Adding in this sub-filter means reads can only have
                        # one mismatch, so they have to be from the right gene
                        # for this to work.
                        cmd += ' subfilter=1'
                    # Add in memory string if specified
                    if xmx:
                        cmd += f' -Xmx{xmx}'

                    # Run the command
                    out, err = run_cmd(cmd=cmd)
                    write_to_logfile(
                        logfile=log,
                        out=out,
                        err=err,
                        cmd=cmd
                    )
                else:
                    if fasta:
                        ax = 'asm5'
                    # If Nanopore FASTQ:
                    else:
                        ax = 'map-ont'
                    # Use minimap2 for Nanopore data or FASTA data
                    cmd = (
                        f'minimap2 --MD -t {threads} -ax {ax} {rmlst_fasta} '
                        f'{unpaired_bait}'
                    )
                    # Convert SAM to sorted BAM in one step to save time/disk
                    # space
                    cmd += (
                        f' | samtools view -@ {threads} -h -bT {rmlst_fasta} -'
                        f' | samtools sort - -@ {threads} -o {sorted_bam}'
                    )

                    # Run the command
                    out, err = run_cmd(cmd=cmd)
                    write_to_logfile(
                        logfile=log,
                        out=out,
                        err=err,
                        cmd=cmd
                    )

        # Ensure the BAM is sorted and indexed
        if not os.path.isfile(sorted_bam):
            # Use pysam to sort and index the BAM
            pysam.sort('-o', sorted_bam, outbam)
        if not os.path.isfile(sorted_bam + '.bai'):
            pysam.index(sorted_bam)

        # Now find number of multi-positions for each rMLST gene/allele
        # combination
        multi_positions = 0

        # Run the BAM parsing in parallel! Some refactoring of the code would
        # likely be a good idea so this isn't quite so ugly, but it works.
        p = multiprocessing.Pool(processes=threads)
        nanopore = True if data_type == 'Nanopore' else False
        nanopore_list = [nanopore] * len(gene_alleles)
        allele_records = SeqIO.to_dict(SeqIO.parse(rmlst_fasta, 'fasta'))
        bamfile_list = [sorted_bam] * len(gene_alleles)
        reference_fasta_list = [rmlst_fasta] * len(gene_alleles)
        fasta_list = [fasta] * len(gene_alleles)
        quality_cutoff_list = [quality_cutoff] * len(gene_alleles)
        base_cutoff_list = [base_cutoff] * len(gene_alleles)
        base_fraction_list = [base_fraction_cutoff] * len(gene_alleles)
        records_list = [allele_records] * len(gene_alleles)
        fastq_records_list = [fastq_records] * len(gene_alleles)
        error_cutoff_list = [error_cutoff] * len(gene_alleles)
        multibase_dict_list = []
        report_write_list = []
        if debug == 'debug':
            for i, gene in enumerate(gene_alleles):
                multibase_dict, report_write = read_contig(
                    contig_name=gene,
                    bamfile_name=bamfile_list[i],
                    reference_fasta=reference_fasta_list[i],
                    allele_records=records_list[i],
                    fastq_records=fastq_records_list[i],
                    quality_cutoff=quality_cutoff_list[i],
                    base_cutoff=base_cutoff_list[i],
                    base_fraction_cutoff=base_fraction_list[i],
                    fasta=fasta_list[i],
                    nanopore=nanopore_list[i],
                    error_cutoff=error_cutoff_list[i]
                )
                multibase_dict_list.append(multibase_dict)
                # Keep representation consistent with parallel branch: a
                # list of lines per gene
                report_write_list.append(report_write.splitlines(True))
        else:
            # Build kwargs list so `read_contig` is invoked with keyword
            # arguments. This is necessary because the function enforces
            # keyword-only parameters (leading ``*`` in the signature).
            kwargs_list = [
                {
                    'contig_name': gene_alleles[i],
                    'bamfile_name': bamfile_list[i],
                    'reference_fasta': reference_fasta_list[i],
                    'allele_records': records_list[i],
                    'fastq_records': fastq_records_list[i],
                    'quality_cutoff': quality_cutoff_list[i],
                    'base_cutoff': base_cutoff_list[i],
                    'base_fraction_cutoff': base_fraction_list[i],
                    'fasta': fasta_list[i],
                    'error_cutoff': error_cutoff_list[i],
                    'nanopore': nanopore_list[i],
                }
                for i in range(len(gene_alleles))
            ]

            # Use map with the dispatch helper which unpacks kwargs
            for multibase_dict, report_write in p.map(
                _read_contig_dispatch, kwargs_list, chunksize=1
            ):
                multibase_dict_list.append(multibase_dict)
                report_write_list.append(report_write.splitlines(True))

            # Close the pool
            p.close()
            p.join()
    except SamtoolsError:
        # Handle Samtools errors gracefully
        pysam_pass = False
        multi_positions = 0
        multibase_dict_list = []
        report_write_list = []

    # Write out report info.
    report_file = os.path.join(
        output_folder,
        sample_name + '_contamination.tsv'
    )

    # Write the contamination report TSV
    with open(report_file, 'w', encoding='utf-8') as r:
        # Contamination report TSV header (tab-separated columns):
        # Gene, Position, RefBase, CongruentSNVs, TotalSNVs, ForwardSNVs,
        # ReverseSNVs, SNVCoverage, TotalCoverage, BaseCutoff, ErrorPercent,
        # PValue, AdjPValue, StrandP, PosP, MeanQual, MeanMapQ
        # Note: the per-position output from `position_details` concatenates
        # the four SNV "read-type" fields (Congruent/Paired/Forward/Reverse)
        # into a single field separated by commas — this means consumers that
        # index columns by numeric position should be careful (prefer matching
        # headers by name when possible).
        r.write(
            'Gene\tPosition\tRefBase\tCongruentSNVs\tTotalSNVs\tForwardSNVs\t'
            'ReverseSNVs\tSNVCoverage\tTotalCoverage\tBaseCutoff\tErrorPercent'
            '\tPValue\tAdjPValue\tStrandP\tPosP\tMeanQual\tMeanMapQ\n'
        )

        # Iterate through the report_write_list and write each line. Items may
        # be lists (parallel branch) or strings (debug branch); normalize to
        # a list of lines before writing.
        for item in report_write_list:
            lines = item if isinstance(item, list) else item.splitlines(True)
            for contamination_info in lines:
                r.write(contamination_info)

    # Total up the number of multibase positions using helper (excludes meta
    # keys)
    multi_positions = count_multibase_positions(
        multibase_dict_list=multibase_dict_list
    )
    logging.debug('Number of contaminating SNVs found: %s', multi_positions)
    # Determine SNP cutoff based on database type
    if cgmlst_db is None:
        snp_cutoff = math.ceil(rmlst_gene_length / 10000) + 1
    elif fasta:
        snp_cutoff = 1
    else:
        snp_cutoff = 10

    # Compute a simple per-sample score: sum of -log10(q) across reported
    # positions (q = adj p-value). Used as a diagnostic.
    sample_score = None
    qs = []

    # Flatten report_write_list into lines for robust parsing. Some entries
    # are lists of lines (parallel branch) and some may be strings.
    all_lines: List[str] = []
    for item in report_write_list:
        if isinstance(item, list):
            all_lines.extend(item)
        elif isinstance(item, str):
            all_lines.extend(item.splitlines(True))

    for line in all_lines:
        # Lines are tab-separated; splitting on tabs avoids breaking fields
        # that contain commas.
        fields = line.rstrip('\n').split('\t')
        # AdjPValue is the 13th column (0-based index 12)
        if len(fields) > 12:
            q_str = fields[12]
            if q_str != 'ND':
                try:
                    qs.append(float(q_str))
                except ValueError:
                    pass
    if qs:
        sample_score = sum([-math.log10(q + 1e-300) for q in qs])

    # Write gene-summary TSV file
    gene_summary_file = os.path.join(
        output_folder,
        sample_name + '_gene_summary.tsv'
    )
    with open(gene_summary_file, 'w', encoding='utf-8') as gf:
        gf.write('Gene\tCombinedP\tGeneScore\tNumSigPositions\tNumPositions\n')
        for multibase_dict in multibase_dict_list:
            for gene, content in multibase_dict.items():
                # Skip meta keys
                if gene.startswith('_'):
                    continue

                # Get stats
                stats = content.get('_gene_stats', {})

                # Skip if no stats
                if not stats:
                    continue

                # Extract relevant stats
                combined_p = stats.get('combined_p')
                gene_score_val = stats.get('gene_score')
                num_sig_positions = stats.get('num_sig_positions', 0)
                num_positions = sum(
                    1 for k in content.keys() if not str(k).startswith('_')
                )
                combined_p_str = (
                    f'{combined_p:0.3e}' if combined_p is not None else 'ND'
                )
                gene_score_str = (
                    f'{gene_score_val:0.3f}' if gene_score_val is not None
                    else 'ND'
                )
                gf.write(
                    f'{gene}\t{combined_p_str}\t{gene_score_str}\t'
                    f'{num_sig_positions}\t{num_positions}\n'
                )

    logging.info(
        'Done! Number of contaminating SNVs found: %s\n', multi_positions
    )

    # Write summary report
    write_output(
        output_report=os.path.join(output_folder, 'confindr_report.tsv'),
        sample_name=sample_name,
        multi_positions=multi_positions,
        genus=genus,
        total_gene_length=rmlst_gene_length,
        snp_cutoff=snp_cutoff,
        database_download_date=database_download_date,
        pysam_pass=pysam_pass,
        sample_score=sample_score,
        use_probabilistic=use_prob_scoring,
        score_threshold=score_threshold
    )

    # Clean up temporary files unless keep_files is set
    if keep_files is False:
        shutil.rmtree(sample_tmp_dir)


def write_output(
    *,  # Enforce keyword-only arguments
    output_report: str,
    sample_name: str,
    multi_positions: int,
    genus: str,
    total_gene_length: int,
    database_download_date: str,
    snp_cutoff: int = 3,
    pysam_pass: bool = True,
    sample_score: Optional[float] = None,
    use_probabilistic: bool = False,
    score_threshold: Optional[float] = None
) -> None:
    """
    Write ConFindr summary report

    Args:
        output_report: Path to CSV/TSV output report file. If it does not
        exist, a header will be written automatically.
        sample_name: Name of the sample.
        multi_positions: Number of multibase positions found for the sample.
        genus: Genus string (may include ':' when multiple genera are present).
        total_gene_length: Total bases examined across genes.
        database_download_date: Date string indicating database download date.
        *  # Enforce keyword-only arguments for the following params
        snp_cutoff: Number of cSNVs required to call a sample contaminated.
        pysam_pass: Whether pysam mapping completed successfully.
        sample_score: Optional per-sample score (diagnostic).
        use_probabilistic: If True and score_threshold supplied, decide
        contamination based on sample_score >= score_threshold.
        score_threshold: Threshold used when use_probabilistic is True.

    Returns:
        None
    """
    # If the report file hasn't been created, make it, with appropriate header.
    if not os.path.isfile(output_report):
        with open(os.path.join(output_report), 'w', encoding='utf-8') as f:
            # Summary report is a TSV: Sample\tGenus\tNumContamSNVs\t
            # ContamStatus\tBasesExamined\tDatabaseDownloadDate\tScore
            f.write(
                'Sample\tGenus\tNumContamSNVs\tContamStatus\t'
                'BasesExamined\tDatabaseDownloadDate\tScore\n'
            )

    # Determine contamination status
    if pysam_pass:
        # Check contamination based on probabilistic or deterministic method
        if use_probabilistic and score_threshold is not None:
            # Probabilistic decision based on sample_score (diagnostic)
            contaminated = bool(
                sample_score is not None and sample_score >= score_threshold
            )
        else:
            contaminated = bool(
                multi_positions >= snp_cutoff or len(genus.split(':')) > 1
            )
    else:
        contaminated = 'Pysam SamtoolsError'
        multi_positions = 'ND'

    # Write out the summary line
    score_str = f'{sample_score:0.3f}' if sample_score is not None else 'ND'
    with open(output_report, 'a+', encoding='utf-8') as f:
        f.write(
            f'{sample_name}\t{genus}\t{multi_positions}\t{contaminated}\t'
            f'{total_gene_length}\t{database_download_date}\t{score_str}\n'
        )


def check_for_databases_and_download(
    *,  # Enforce keyword-only arguments
    database_location: str
) -> None:
    """
    Check for necessary ConFindr databases, download if not present.

    Args:
        database_location: Path to ConFindr database folder.

    Returns:
        None
    """
    # Check for the files necessary - should have rMLST_combined.fasta,
    # gene_allele.txt, profiles.txt, and refseq.msh
    necessary_files = [
        'Escherichia_db_cgderived.fasta',
        'Listeria_db_cgderived.fasta',
        'Salmonella_db_cgderived.fasta',
        'refseq.msh'
    ]

    # Also check for optional rMLST files
    optional_files = [
        'rMLST_combined.fasta',
        'gene_allele.txt',
        'profiles.txt'
    ]

    # Initialize flag
    all_files_present = True

    # Iterate through necessary files
    for necessary_file in necessary_files:
        if not os.path.isfile(os.path.join(database_location, necessary_file)):
            logging.warning('Could not find %s', necessary_file)
            all_files_present = False

    # Download necessary files if not present
    if not all_files_present:
        logging.warning(
            'Databases not present - downloading basic databases now...'
        )

        # Create database directory if it doesn't exist
        os.makedirs(database_location, exist_ok=True)

        # Download necessary database files
        download_mash_sketch(output_folder=database_location)
        download_cgmlst_derived_data(output_folder=database_location)

    # Check for optional files
    optional_files_present = True
    for optional_file in optional_files:
        if not os.path.isfile(os.path.join(database_location, optional_file)):
            optional_files_present = False
    if not optional_files_present:
        logging.warning(
            'Did not find rMLST databases, if you want to use ConFindr on '
            'genera other than Listeria, Salmonella, and Escherichia, '
            'you\'ll need to download them. Instructions are available at '
            'https://olc-bioinformatics.github.io/ConFindr/install/'
            '#downloading-confindr-databases\n'
        )


def check_valid_base_fraction(
    *,  # Enforce keyword-only arguments
    base_fraction: Optional[float]
) -> bool:
    """
    Validate that a base fraction is None or between 0 and 1 inclusive.

    Args:
        base_fraction: The fraction to validate (or None).

    Returns:
        True if valid, False otherwise.
    """
    if base_fraction is None:
        return True
    if 0 <= base_fraction <= 1:
        return True

    return False


def check_acceptable_xmx(
    *,  # Enforce keyword-only arguments
    xmx_string: str
) -> bool:
    """
    Validate a BBTools -Xmx memory specification string.

    Args:
        xmx_string: Memory string such as '20g', '800m', '1024K'.

    Returns:
        True if acceptable, False otherwise.
    """
    # Initialize flag
    acceptable_xmx = True

    # Set of acceptable suffixes
    acceptable_suffixes = ['K', 'M', 'G']

    # Check suffix
    if xmx_string[-1].upper() not in acceptable_suffixes:
        acceptable_xmx = False
        logging.error(
            'ERROR: Memory must be specified as K (kilobytes), M (megabytes), '
            'or G (gigabytes). Your specified suffix was %s.', xmx_string[-1]
        )

    # Check that the rest is an integer
    if '.' in xmx_string:
        acceptable_xmx = False
        logging.error(
            'ERROR: Xmx strings must be integers, floating point numbers '
            'are not accepted.'
        )

    # Check that the rest is digits
    if not str.isdigit(xmx_string[:-1]):
        acceptable_xmx = False
        logging.error(
            'ERROR: The amount of memory requested was not an integer.'
        )

    return acceptable_xmx


def get_version() -> str:
    """
    Get ConFindr version string.

    Returns:
        Version string.
    """
    try:
        version = (
            f'ConFindr {pkg_resources.get_distribution("confindr").version}'
        )
    except pkg_resources.DistributionNotFound:
        version = 'ConFindr (Unknown version)'
    return version
