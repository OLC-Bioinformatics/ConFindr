#!/usr/bin/env python
from confindr_src.wrappers import bbtools, mash
from Bio import SeqIO
from pysam.utils import SamtoolsError
from itertools import chain
from statistics import mean
import multiprocessing
import urllib.request
import pkg_resources
import numpy as np
import subprocess
import logging
import shutil
import tarfile
import pysam
import glob
import gzip
import math
import csv
import os


def download_mash_sketch(output_folder):
    logging.info('Downloading mash refseq sketch...')
    urllib.request.urlretrieve('https://github.com/OLC-Bioinformatics/ConFindr/raw/master/refseq_sketch/refseq.msh',
                               os.path.join(output_folder, 'refseq.msh'))


def download_cgmlst_derived_data(output_folder):
    logging.info('Downloading cgMLST-derived data for Escherichia, Salmonella, and Listeria...')
    urllib.request.urlretrieve('https://ndownloader.figshare.com/files/14771267',
                               os.path.join(output_folder, 'confindr_db.tar.gz'))
    confindr_tar = os.path.join(output_folder, 'confindr_db.tar.gz')
    tar = tarfile.open(confindr_tar)
    tar.extractall(path=output_folder)
    tar.close()
    os.remove(confindr_tar)
    index(output_folder=output_folder,
          genera=['Escherichia', 'Listeria', 'Salmonella'],
          cgderived=True)


def index(output_folder, genera, cgderived=False):
    """

    :param output_folder:
    :param genera:
    :param cgderived:
    :return:
    """
    for predominant_genus in genera:
        if cgderived:
            sample_database = os.path.join(output_folder, '{pg}_db_cgderived.fasta'.format(pg=predominant_genus))
        else:
            sample_database = os.path.join(output_folder, '{pg}_db.fasta'.format(pg=predominant_genus))

        if not os.path.isfile(sample_database):

            if os.path.isfile(os.path.join(output_folder, 'gene_allele.txt')) and \
                    os.path.isfile(os.path.join(output_folder, 'rMLST_combined.fasta')):
                logging.info('Setting up rMLST genus-specific database for genus {pg}...'
                             .format(pg=predominant_genus))
                allele_list = find_genus_specific_allele_list(os.path.join(output_folder, 'gene_allele.txt'),
                                                              predominant_genus)
                # Create the allele-specific database
                setup_allelespecific_database(fasta_file=sample_database,
                                              database_folder=output_folder,
                                              allele_list=allele_list)
        # Perform the necessary samtools and KMA indexing
        index_databases(sample_database=sample_database)


def run_cmd(cmd):
    """
    Runs a command using subprocess, and returns both the stdout and stderr from that command
    If exit code from command is non-zero, raises subprocess.CalledProcessError
    :param cmd: command to run as a string, as it would be called on the command line
    :return: out, err: Strings that are the stdout and stderr from the command called.
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd=cmd)
    return out, err


def write_to_logfile(logfile, out, err, cmd):
    """
    Writes stdout, stderr, and a command to a logfile
    :param logfile: Path to file to write output to.
    :param out: Stdout of program called, as a string
    :param err: Stderr of program called, as a string
    :param cmd: command that was used
    """
    with open(logfile, 'a+') as outfile:
        outfile.write('Command used: {}\n\n'.format(cmd))
        outfile.write('STDOUT: {}\n\n'.format(out))
        outfile.write('STDERR: {}\n\n'.format(err))


def dependency_check(dependency):
    """
    Uses shutil to check if a dependency is installed (won't check version of anything - just presence)
    :param dependency: The dependency as it would be called on the command line (i.e. for blastn, would be blastn)
    :return: True if dependency is present, False if it is not found.
    """
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

    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in sorted(fastq_files):
        if forward_id in name and os.path.isfile(name.replace(forward_id, reverse_id)):
            pair_list.append([name, name.replace(forward_id, reverse_id)])
    return pair_list


def find_unpaired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2', find_fasta=False):
    """
    Looks at a directory to find unpaired fastq files.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default _R1.
    :param reverse_id: Identifier for forward reads. Default _R2.
    :param find_fasta: If False, will look for fastq files. Otherwise, looks for Fasta files.
    :return: List of files that appear to be unpaired reads.
    """
    read_list = list()
    if find_fasta is False:
        fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    else:
        # Very misnamed!
        fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*a*'))
    for name in sorted(fastq_files):
        # Iterate through files, adding them to our list of unpaired reads if:
        # 1) They don't have the forward identifier or the reverse identifier in their name.
        # 2) They have forward but the reverse isn't there.
        # 3) They have reverse but the forward isn't there.
        if forward_id not in name and reverse_id not in name:
            read_list.append([name])
        elif forward_id in name and not os.path.isfile(name.replace(forward_id, reverse_id)):
            read_list.append([name])
        elif reverse_id in name and not os.path.isfile(name.replace(reverse_id, forward_id)):
            read_list.append([name])
    return read_list


def find_genus_specific_allele_list(profiles_file, target_genus):
    """
    A new way of making our specific databases: Make our profiles file have lists of every gene/allele present for
    each genus instead of just excluding a few genes for each. This way, should have much smaller databases
    while managing to make ConFindr a decent bit faster (maybe)
    :param profiles_file: Path to profiles file.
    :param target_genus: Genus you want to make a custom database for (STR)
    :return: List of gene/allele combinations that should be part of species-specific database.
    """
    alleles = list()
    with open(profiles_file) as f:
        lines = f.readlines()
    for line in lines:
        line = line.rstrip()
        genus = line.split(':')[0]
        if genus == target_genus:
            alleles = line.split(':')[1].split(',')[:-1]
    return alleles


def setup_allelespecific_database(fasta_file, database_folder, allele_list):
    """
    Since some genera have some rMLST genes missing, or two copies of some genes, genus-specific databases are needed.
    This will take only the alleles known to be part of each genus and write them to a genus-specific file.
    :param database_folder: Path to folder where rMLST_combined is stored.
    :param fasta_file: Path to fasta file to write allele-specific database to.
    :param allele_list: allele list generated by find_genus_specific_allele_list
    """
    rmlst_index = SeqIO.index(os.path.join(database_folder, 'rMLST_combined.fasta'), 'fasta')
    seqs = list()
    for s in allele_list:
        try:
            seqs.append(rmlst_index[s])
        except KeyError:
            logging.warning('Tried to add {} to allele-specific database, but could not find it.'.format(s))
    try:
        SeqIO.write(seqs, fasta_file, 'fasta')
    except FileNotFoundError:
        pass


def find_cross_contamination(databases, reads, sample_name, tmpdir='tmp', log='log.txt', threads=1,
                             min_matching_hashes=40):
    """
    Uses mash to find out whether or not a sample has more than one genus present, indicating cross-contamination.
    :param databases: A databases folder, which must contain refseq.msh, a mash sketch that has one representative
    per genus from refseq.
    :param reads: Relative path(s) to either unpaired (type STR) or paired (type LIST) FASTQ reads
    :param sample_name:
    :param tmpdir: Temporary directory to store mash result files in.
    :param log: Logfile to write to.
    :param threads: Number of threads to run mash with.
    :param min_matching_hashes: Minimum number of matching hashes in a MASH screen in order for a genus to be
    considered present in a sample. Default is 40
    :return: cross_contam: a bool that is True if more than one genus is found, and False otherwise.
    :return: genera_present: A string. If only one genus is found, string is NA. If more than one genus is found,
    the string is a list of genera present, separated by colons (i.e. for Escherichia and Salmonella found, string would
    be 'Escherichia:Salmonella'
    """
    genera_present = list()
    # Only run the MASH analyses if the screen.tab output file does not already exist
    screen_file = os.path.join(tmpdir, '{sn}_screen.tab'.format(sn=sample_name))
    if not os.path.isfile(screen_file):
        if type(reads) is str:
            out, err, cmd = mash.screen('{database}/refseq.msh'.format(database=databases), reads,
                                        threads=threads,
                                        w='',
                                        i='0.85',
                                        output_file=screen_file,
                                        returncmd=True)
        else:
            out, err, cmd = mash.screen('{database}/refseq.msh'.format(database=databases), reads[0],
                                        reads[1],
                                        threads=threads,
                                        w='',
                                        i='0.85',
                                        output_file=screen_file,
                                        returncmd=True)
        write_to_logfile(log, out, err, cmd)
    screen_output = mash.read_mash_screen(os.path.join(tmpdir, '{sn}_screen.tab'.format(sn=sample_name)))
    for item in screen_output:
        mash_genus = item.query_id.split('/')[-3]
        if 'Shigella' in mash_genus:
            mash_genus = 'Escherichia'
        matching_hashes = int(item.shared_hashes.split('/')[0])
        # Only add the genus to the genera_present list of the number of matching hashes exceeds the cutoff
        if matching_hashes >= min_matching_hashes:
            if mash_genus not in genera_present:
                genera_present.append(mash_genus)
    if len(genera_present) == 1:
        genera_present = genera_present[0]
    elif len(genera_present) == 0:
        genera_present = 'ND'
    else:
        tmpstr = ''
        for mash_genus in genera_present:
            tmpstr += mash_genus + ':'
        genera_present = tmpstr[:-1]
    return genera_present


def number_of_bases_above_threshold(high_quality_base_count, base_count_cutoff=2, base_fraction_cutoff=None):
    """
    Finds if a site has at least two bases of high quality, enough that it can be considered
    fairly safe to say that base is actually there.
    :param high_quality_base_count: Dictionary of count of HQ bases at a position where key is base and values is the
    count of that base.
    :param base_count_cutoff: Number of bases needed to support multiple allele presence.
    :param base_fraction_cutoff: Fraction of bases needed to support multiple allele presence.
    :return: True if site has at least base_count_cutoff/base_fraction_cutoff bases, False otherwise
    (changeable by user)
    """

    # make a dict by dictionary comprehension where values are True or False for each base depending on whether the
    # count meets the threshold.
    # Method differs depending on whether absolute or fraction cutoff is specified
    if base_fraction_cutoff:
        total_hq_base_count = sum(high_quality_base_count.values())
        bases_above_threshold = {base: float(count) / total_hq_base_count >= base_fraction_cutoff and count
                                 >= base_count_cutoff for (base, count) in high_quality_base_count.items()}
    else:
        bases_above_threshold = {base: count >= base_count_cutoff for (base, count) in high_quality_base_count.items()}

    # True is equal to 1 so sum of the number of Trues in the bases_above_threshold dict is the number of bases
    # passing threshold
    return sum(bases_above_threshold.values())


def parse_bam(bamfile_name, contig_name, pysam_fasta):
    """
    Use pysam to load the sorted BAM-formatted file, and subsequently the gene-specific read pileup
    :param bamfile_name: Name and path of the sorted BAM file
    :param contig_name: Name of the current gene
    :param pysam_fasta: pysam.FastaFile object created from the reference gene sequence
    :return: bamfile: pysam.AlignmentFile object of the BAM-formatted file
    :return pileup: A pysam pileup object created from the BAM file
    """
    # Load the sorted BAM-formatted file using pysam
    bamfile = pysam.AlignmentFile(bamfile_name, 'rb')
    # These parameters seem to be fairly undocumented with pysam, but I think that they should make the output
    # that I'm getting to match up with what I'm seeing in Tablet.
    pileup = bamfile.pileup(contig_name,
                            stepper='samtools',
                            ignore_orphans=False,
                            fastafile=pysam_fasta,
                            min_base_quality=0)
    return bamfile, pileup


def characterise_read(column, reference_sequence, fastq_records, quality_cutoff,
fasta=False, nanopore=False):
    """
    Parses a column to characterize all the bases present. Determines the number of bases that fit certain criteria
    :param column: A pileupColumn generated by pysam
    :param reference_sequence: String of the FASTA reference gene sequence
    :param fastq_records: Dictionary of SeqIO records parsed from filtered FASTQ reads
    :param quality_cutoff: Desired min phred quality for a base in order to be counted towards a multi-allelic column
    If specified, both the base_cutoff and base_fraction_cutoff will have to be met
    :param fasta: Boolean of whether the samples are in FASTQ or FASTA format
    :param nanopore: Boolean of whether Nanopore reads were provided. Default is False.
    :return: filtered_read_dict: Dictionary with of base types: read direction: count
    :return: qualities: List of all phred scores in the pileup
    """
    # Initialise a dictionary to store the base details
    filtered_read_dict = {
        'congruent_SNV': dict(),
        'congruent_ref': dict(),
        'forward_SNV_reverse_SNV1': dict(),
        'reverse_SNV_forward_SNV1': dict(),
        'forward_SNV_reverse_ref': dict(),
        'reverse_SNV_forward_ref': dict(),
        'forward_SNV_reverse_UM_QF': dict(),
        'forward_ref_reverse_UM_QF': dict(),
        'forward_quality_filtered': dict(),
        'reverse_SNV_forward_UM_QF': dict(),
        'reverse_ref_forward_UM_QF': dict(),
        'reverse_quality_filtered': dict()
    }
    # Initialise a dictionary to store the details parsed from the pileup
    unfiltered_read_details = dict()
    # Initialise a list to store all the phred scores for the bases passing filter in the column
    qualities = list()
    # Iterate through every read present in the column of the pileup
    for read in column.pileups:
        # Not entirely sure why this is sometimes None, but it causes bad stuff
        if read.query_position is not None:
            #  Initialise the read name in the dictionary as required
            if read.alignment.qname not in unfiltered_read_details:
                unfiltered_read_details[read.alignment.qname] = dict()
            # Extract the sequence of the base in the read
            query_base = read.alignment.query_sequence[read.query_position]
            # Extract the sequence of the base in the reference gene
            ref_base = reference_sequence[column.pos]
            # Create a boolean of whether the query base matches the reference base (is it a SNV?)
            match = query_base == ref_base
            # Read names in the pileup have the direction removed - add it back for future parsing. Not an issue for
            # FASTA files
            if not fasta and not nanopore:
                read_name = read.alignment.qname.split(' ')[0] + '/1' if read.alignment.is_read1 else \
                    read.alignment.qname.split(' ')[0] + '/2'
            else:
                read_name = read.alignment.qname
            # Extract the phred quality score from the FASTQ records
            quality = fastq_records[read_name].letter_annotations["phred_quality"][read.query_position]
            # Initialise a dictionary to store the ran
            range_dict = dict()
            # Determine whether there are SNVs clustered together - they will be discarded from the analysis
            # Iterate through a range of the five positions preceding, and five subsequent
            # positions
            for iterator, contig_pos in enumerate(
                    chain(range(column.pos - 5, column.pos), range(column.pos + 1, column.pos + 6))):
                # Ensure that the contig position being examined isn't beyond the length of the gene
                if 0 <= contig_pos < len(reference_sequence) - 1:
                    # Calculate the read position corresponding to the current column position
                    read_pos = list(chain(range(read.query_position - 5, read.query_position),
                                          range(read.query_position + 1, read.query_position + 6)))[iterator]
                    # Ensure that read position being examined isn't beyond the length of the read
                    if 0 <= read_pos < read.alignment.query_alignment_end - 1:
                        # Populate the dictionary with the calculated positions
                        range_dict[contig_pos] = read_pos
            # Initialise a boolean of whether the current base passes filters, and should be added to the dictionary
            add_base = True
            # Iterate through the range_dict to extract the sequence of the bases in the range for both the read and
            # the reference sequence
            for contig_pos, read_pos in range_dict.items():
                try:
                    # Extract the sequence of the base
                    reference_base = reference_sequence[contig_pos]
                    base = read.alignment.query_sequence[read_pos]
                    # If any of the downstream or upstream bases do not match, set the boolean to False
                    if not match and reference_base != base:
                        add_base = False
                except KeyError:
                    pass
            # Populate the dictionary only if there are no other SNVs within five downstream and five upstream bases
            if add_base:
                unfiltered_read_details[read.alignment.qname][read.alignment.is_read1] = {
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
    # Parse the unfiltered reads to characterise the bases
    for read_base, dir_dict in unfiltered_read_details.items():
        # Check to see if paired reads are present at this position
        if len(dir_dict) > 1:
            # SNV in both forward and reverse reads
            if not dir_dict[True]['match'] and not dir_dict[False]['match']:
                # Same SNV sequence - don't quality filter as the reads agreeing acts as a quality check
                if dir_dict[True]['qbase'] == dir_dict[False]['qbase']:
                    # Forward and reverse SNV - add two matches (for forward and reverse reads)
                    if dir_dict[True]['qbase'] not in filtered_read_dict['congruent_SNV']:
                        filtered_read_dict['congruent_SNV'][dir_dict[True]['qbase']] = 2
                    else:
                        filtered_read_dict['congruent_SNV'][dir_dict[True]['qbase']] += 2
                # Different SNV sequences
                else:
                    # Both SNV sequences pass quality
                    if dir_dict[True]['qual'] >= quality_cutoff and dir_dict[False]['qual'] >= quality_cutoff:
                        # Forward SNV1
                        if dir_dict[True]['qbase'] not in filtered_read_dict['forward_SNV_reverse_SNV1']:
                            filtered_read_dict['forward_SNV_reverse_SNV1'][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict['forward_SNV_reverse_SNV1'][dir_dict[True]['qbase']] += 1
                        # Reverse SNV2
                        if dir_dict[False]['qbase'] not in filtered_read_dict['reverse_SNV_forward_SNV1']:
                            filtered_read_dict['reverse_SNV_forward_SNV1'][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict['reverse_SNV_forward_SNV1'][dir_dict[False]['qbase']] += 1
                    # Only the forward reads pass quality
                    elif dir_dict[True]['qual'] >= quality_cutoff:
                        # Forward SNV reverse QF
                        if dir_dict[True]['qbase'] not in filtered_read_dict['forward_SNV_reverse_UM_QF']:
                            filtered_read_dict['forward_SNV_reverse_UM_QF'][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict['forward_SNV_reverse_UM_QF'][dir_dict[True]['qbase']] += 1
                        # Reverse QF
                        if dir_dict[False]['qbase'] not in filtered_read_dict['reverse_quality_filtered']:
                            filtered_read_dict['reverse_quality_filtered'][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict['reverse_quality_filtered'][dir_dict[False]['qbase']] += 1
                    # Only the reverse reads pass quality
                    elif dir_dict[False]['qual'] >= quality_cutoff:
                        # Reverse SNV forward QF
                        if dir_dict[False]['qbase'] not in filtered_read_dict['reverse_SNV_forward_UM_QF']:
                            filtered_read_dict['reverse_SNV_forward_UM_QF'][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict['reverse_SNV_forward_UM_QF'][dir_dict[False]['qbase']] += 1
                        # Forward QF
                        if dir_dict[True]['qbase'] not in filtered_read_dict['forward_quality_filtered']:
                            filtered_read_dict['forward_quality_filtered'][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict['forward_quality_filtered'][dir_dict[True]['qbase']] += 1
                    # Neither forward nor reverse reads pass quality
                    else:
                        # Forward QF
                        if dir_dict[True]['qbase'] not in filtered_read_dict['forward_quality_filtered']:
                            filtered_read_dict['forward_quality_filtered'][dir_dict[True]['qbase']] = 1
                        else:
                            filtered_read_dict['forward_quality_filtered'][dir_dict[True]['qbase']] += 1
                        # Reverse QF
                        if dir_dict[False]['qbase'] not in filtered_read_dict['reverse_quality_filtered']:
                            filtered_read_dict['reverse_quality_filtered'][dir_dict[False]['qbase']] = 1
                        else:
                            filtered_read_dict['reverse_quality_filtered'][dir_dict[False]['qbase']] += 1
            # SNV in forward read only
            elif not dir_dict[True]['match'] and dir_dict[False]['match']:
                # Since only the forward read supports the SNV, quality filter
                if dir_dict[True]['qual'] >= quality_cutoff:
                    if dir_dict[True]['qbase'] not in filtered_read_dict['forward_SNV_reverse_ref']:
                        filtered_read_dict['forward_SNV_reverse_ref'][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict['forward_SNV_reverse_ref'][dir_dict[True]['qbase']] += 1
                    # Since the reverse base matches the reference, don't quality filter
                    if dir_dict[False]['qbase'] not in filtered_read_dict['forward_SNV_reverse_ref']:
                        filtered_read_dict['forward_SNV_reverse_ref'][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict['forward_SNV_reverse_ref'][dir_dict[False]['qbase']] += 1
                # Reverse ref forward QF
                else:
                    if dir_dict[False]['qbase'] not in filtered_read_dict['reverse_ref_forward_UM_QF']:
                        filtered_read_dict['reverse_ref_forward_UM_QF'][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict['reverse_ref_forward_UM_QF'][dir_dict[False]['qbase']] += 1
            # SNV in reverse read only
            elif dir_dict[True]['match'] and not dir_dict[False]['match']:
                # Quality filter
                if dir_dict[False]['qual'] >= quality_cutoff:
                    if dir_dict[False]['qbase'] not in filtered_read_dict['reverse_SNV_forward_ref']:
                        filtered_read_dict['reverse_SNV_forward_ref'][dir_dict[False]['qbase']] = 1
                    else:
                        filtered_read_dict['reverse_SNV_forward_ref'][dir_dict[False]['qbase']] += 1
                    # Since the forward base matches the reference, don't quality filter
                    if dir_dict[True]['qbase'] not in filtered_read_dict['reverse_SNV_forward_ref']:
                        filtered_read_dict['reverse_SNV_forward_ref'][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict['reverse_SNV_forward_ref'][dir_dict[True]['qbase']] += 1
                # Forward ref reverse QF
                else:
                    if dir_dict[True]['qbase'] not in filtered_read_dict['forward_ref_reverse_UM_QF']:
                        filtered_read_dict['forward_ref_reverse_UM_QF'][dir_dict[True]['qbase']] = 1
                    else:
                        filtered_read_dict['forward_ref_reverse_UM_QF'][dir_dict[True]['qbase']] += 1
            # Both reads match the reference sequence - don't quality filter, and add two matches (for forward and
            # reverse reads)
            else:
                if dir_dict[True]['qbase'] not in filtered_read_dict['congruent_ref']:
                    filtered_read_dict['congruent_ref'][dir_dict[True]['qbase']] = 2
                else:
                    filtered_read_dict['congruent_ref'][dir_dict[True]['qbase']] += 2
        # Either the reads are unpaired, or only a single read aligns to this position on the gene (no overlap)
        else:
            for direction in dir_dict:
                # SNV supported by a single read
                if not dir_dict[direction]['match']:
                    if dir_dict[direction]['qual'] >= quality_cutoff:
                        # Forward
                        if direction:
                            if dir_dict[direction]['qbase'] not in filtered_read_dict['forward_SNV_reverse_UM_QF']:
                                filtered_read_dict['forward_SNV_reverse_UM_QF'][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict['forward_SNV_reverse_UM_QF'][dir_dict[direction]['qbase']] += 1
                        # Reverse
                        else:
                            if dir_dict[direction]['qbase'] not in filtered_read_dict['reverse_SNV_forward_UM_QF']:
                                filtered_read_dict['reverse_SNV_forward_UM_QF'][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict['reverse_SNV_forward_UM_QF'][dir_dict[direction]['qbase']] += 1
                    else:
                        # Forward QF
                        if dir_dict[direction]['qbase'] not in filtered_read_dict['forward_quality_filtered']:
                            filtered_read_dict['forward_quality_filtered'][dir_dict[direction]['qbase']] = 1
                        else:
                            filtered_read_dict['forward_quality_filtered'][dir_dict[direction]['qbase']] += 1

                # Match to the reference sequence supported by a single read
                else:
                    if dir_dict[direction]['qual'] >= quality_cutoff:
                        # Forward
                        if direction:
                            if dir_dict[direction]['qbase'] not in filtered_read_dict['forward_ref_reverse_UM_QF']:
                                filtered_read_dict['forward_ref_reverse_UM_QF'][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict['forward_ref_reverse_UM_QF'][dir_dict[direction]['qbase']] += 1
                        # Reverse
                        else:
                            if dir_dict[direction]['qbase'] not in filtered_read_dict['reverse_ref_forward_UM_QF']:
                                filtered_read_dict['reverse_ref_forward_UM_QF'][dir_dict[direction]['qbase']] = 1
                            else:
                                filtered_read_dict['reverse_ref_forward_UM_QF'][dir_dict[direction]['qbase']] += 1
                    else:
                        # Reverse QF
                        if dir_dict[direction]['qbase'] not in filtered_read_dict['reverse_quality_filtered']:
                            filtered_read_dict['reverse_quality_filtered'][dir_dict[direction]['qbase']] = 1
                        else:
                            filtered_read_dict['reverse_quality_filtered'][dir_dict[direction]['qbase']] += 1
    return filtered_read_dict, qualities


def determine_cutoff(qualities, reference_sequence, base_cutoff, error_cutoff=1.0):
    """
    Determine the base cutoff value to use for SNV determination
    :param qualities: List of phred scores for every base in the pileup
    :param reference_sequence: String of the FASTA reference gene sequence
    :param error_cutoff: Float of the error cutoff value to use. Default is 1.0%
    :return: base_cutoff: Int of calculated gene-specific cutoff value to use
    :return: The percentage error value corresponding to the base_cutoff value
    """
    depth_prob = 0
    # Only find the cutoff if there are qualities (the gene is present in the sample)
    if qualities:
        # Convert the gene-specific phred score to its probability e.g phred score of 20 is 1x10^(-20/10) = 1x10^(-2)
        # = 0.01
        phred_prob = 1 * 10 ** (-(mean(qualities) / 10))
        # Determine the probability that any base in the pileup will be incorrect. Multiply the phred probability by
        # the average depth of the sample (number of bases / length of the gene)
        depth_prob = phred_prob * len(qualities) / len(reference_sequence)
        # Find the lowest value of the base cutoff that gives a false positive error value under the error cutoff value
        # for the entire length of the gene. Take the depth_prob to the power of the current base_cutoff value.
        # Multiply this by the length of the gene to find find that value for the entire gene, and by 100 to convert
        # this to a percent value. e.g. 0.01^2 * 1674 * 100 = 16.74%. Note that this is above the default 1% threshold.
        # The base_cutoff will be incremented by one, and the calculation will be run again: 0.01^3 * 1674 * 100 =
        # 0.1674%, which is below the threshold, and would pass
        while depth_prob ** base_cutoff * len(reference_sequence) > error_cutoff * 100:
            base_cutoff += 1
    return base_cutoff, depth_prob ** base_cutoff * 100 * len(reference_sequence)


def find_multibase_positions(ref_base, filtered_read_dict, base_cutoff, base_fraction_cutoff):
    """

    :param ref_base: String of the sequence of the reference gene at the current position
    :param filtered_read_dict: Dictionary with of base types: read direction: count
    :param base_cutoff: Integer of the number of identical mismatches in a column required for a SNV call
    :param base_fraction_cutoff: Float fraction of bases necessary to support a SNV call
    :return: snv_dict: Dictionary of characterised bases of the current column e.g. total, number matching reference
    sequence, number of SNVs in both forward and reverse reads, etc.
    :return: passing_snv_dict: Dictionary summarising counts of categories of bases e.g. congruent, forward, reverse,
    paired
    :return: total_coverage: Integer of the total number of bases passing filter in the pileup
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
    # Initialise a dictionary to store the count for each individual nucleotide at this position e.g. A:24, G:2
    base_count = dict()
    # Iterate through the categories e.g. congruent_ref in the dictionary
    for category, base_dict in filtered_read_dict.items():
        # Iterate through the sequence of each query base, and the corresponding count
        for base, count in base_dict.items():
            # Do not process the base if it has been flagged as 'filtered'
            if 'filtered' not in category:
                # Populate the base counting dictionary with the base sequence and increment the count
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
        'congruent': dict(),
        'forward': dict(),
        'reverse': dict(),
        'paired': dict()
    }
    # Boolean of whether there are bases passing filter, and the passing_snv_dict should be used
    return_dict = False
    # Iterate through the categories in the dictionary
    for category, base_dict in filtered_read_dict.items():
        # Iterate through each query base and its corresponding count
        for base, count in base_dict.items():
            # Congruent SNVs
            if 'congruent' in category and base != ref_base:
                # If the base_cutoff_fraction has been provided, use it in making SNV calls
                if base_fraction_cutoff:
                    # Ensure that the number of SNVs in the pileup is greater than the base_cutoff and the
                    # fraction of SNVs in the pileup is greater than the base_cutoff_fraction
                    if float(base_count[base] / snv_dict['total']) >= base_fraction_cutoff and base_count[base] \
                            >= base_cutoff:
                        # Update the summary dictionary with the base sequence and count
                        if base not in passing_snv_dict['forward']:
                            passing_snv_dict['congruent'][base] = count
                        else:
                            passing_snv_dict['congruent'][base] += count
                        return_dict = True
                # If base_cutoff_fraction is not supplied, only the number of SNVs in the pileup must be greater than
                # the base_cutoff value in order for a SNV call
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['congruent']:
                            passing_snv_dict['congruent'][base] = count
                        else:
                            passing_snv_dict['congruent'][base] += count
                    return_dict = True
            # SNVs in the forward read
            if category.startswith('forward_SNV') and base != ref_base:
                if base_fraction_cutoff:
                    if float(base_count[base] / snv_dict['total']) >= base_fraction_cutoff and base_count[base] \
                            >= base_cutoff:
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
                    if float(base_count[base] / snv_dict['total']) >= base_fraction_cutoff and base_count[base] \
                            >= base_cutoff:
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
                    if float(base_count[base] / snv_dict['total']) >= base_fraction_cutoff and base_count[base] \
                            >= base_cutoff:
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
    if return_dict:
        return snv_dict, passing_snv_dict, total_coverage
    else:
        return snv_dict, dict(), total_coverage


def position_details(actual_position, passing_snv_dict, contig_name, ref_base, total_coverage, base_cutoff, error_perc):
    """
    Create a string summarising the ConFindr results for current gene
    :param actual_position: Integer of the adjusted position of the current base
    :param passing_snv_dict: Dictionary summarising counts of categories of bases
    :param contig_name: Name of current gene
    :param ref_base: Sequence of the current reference base
    :param total_coverage: Integer of the total depth at this position
    :param base_cutoff: Integer of the SNV depth cutoff value
    :param error_perc: Float of the calculated error value
    :return: to_write: String containing formatted outputs to be entered into the final report
    """
    # List of the base categories present in passing_snv_dict
    read_types = ['congruent', 'paired', 'forward', 'reverse']
    # Update the to_write string with basic information
    to_write = '{contig_name},{pos},{ref_base},'.format(contig_name=contig_name,
                                                        pos=actual_position,
                                                        ref_base=ref_base)
    # Initialise the coverage value to zero
    snv_coverage = 0
    # Iterate through all the read types
    for read_type in read_types:
        # Boolean of whether a semi-colon needs to be added to the string to separate all the base sequences
        semi_colon = False
        for base, coverage in passing_snv_dict[read_type].items():
            # Ensure that the base isn't empty
            if base:
                if semi_colon:
                    to_write += ';'
                to_write += '{base}:{depth}'.format(base=base,
                                                    depth=coverage)
                semi_colon = True
                # Increment the coverage only if the read type isn't 'paired'
                if read_type != 'paired':
                    snv_coverage += coverage
        to_write += ','
    # Format the error_perc to be more readable
    error_perc = '{:0.2f}'.format(error_perc) if error_perc else 'ND'
    to_write += '{snv_cov},{total_cov},{base_cutoff},{error_perc}\n'\
        .format(snv_cov=snv_coverage,
                total_cov=total_coverage,
                base_cutoff=base_cutoff,
                error_perc=error_perc)
    return to_write


def read_contig(contig_name, bamfile_name, reference_fasta, allele_records, 
                fastq_records, quality_cutoff=20, base_cutoff=None, 
                base_fraction_cutoff=None, fasta=False, error_cutoff=1.0,
                nanopore=False):
    """
    Examines a contig to find if there are positions where more than one base is present.
    :param contig_name: Name of contig as a string.
    :param bamfile_name: Full path to bamfile. Must be sorted/indexed
    :param reference_fasta: Full path to fasta file that was used to generate the bamfile.
    :param allele_records:
    :param fastq_records: Dictionary of SeqIO records parsed from filtered FASTQ reads
    :param quality_cutoff: Bases must have at least this phred score to be considered (INT)
    :param base_cutoff: At least this many bases must support a minor variant (INT)
    :param base_fraction_cutoff: At least this percentage of bases must support minor variant (FLOAT)
    :param fasta: Boolean on whether the samples are in FASTA format. Default is False
    :param error_cutoff: Float of the error cutoff value to use. Default is 1.0%
    :param nanopore: Boolean of whether Nanopore reads were provided. Default is False.
    :return: Dictionary of positions where more than one base is present. Keys are contig name, values are positions
    """
    pysam_fasta = pysam.FastaFile(reference_fasta)
    multibase_position_dict = dict()
    to_write = str()
    # If analysing FASTA files, a single base difference is all that is expected
    if fasta:
        base_cutoff = 1
    reference_sequence = str(allele_records[contig_name].seq)
    # Parse the BAM file with pysam to create AlignmentFile, and AlignmentFile.pileup objects
    bamfile, pileup = parse_bam(bamfile_name=bamfile_name,
                                contig_name=contig_name,
                                pysam_fasta=pysam_fasta)
    filtered_read_dict = dict()
    quality_list = list()
    for i, column in enumerate(pileup):
        filtered_reads, qualities = characterise_read(column=column,
                                                      reference_sequence=reference_sequence,
                                                      fastq_records=fastq_records,
                                                      quality_cutoff=quality_cutoff,
                                                      fasta=fasta,
                                                      nanopore=nanopore)
        filtered_read_dict[i] = filtered_reads
        quality_list += qualities
    # Initialise the calculated error percentage to zero
    error_perc = None
    # If the base_cutoff is not manually specified on the command-line, adjust it based
    # upon overall sequence quality
    if base_cutoff != 3:
        base_cutoff, error_perc = determine_cutoff(qualities=quality_list,
                                                   reference_sequence=reference_sequence,
                                                   base_cutoff = base_cutoff,
                                                   error_cutoff=error_cutoff)
    bamfile.close()
    # It seems that the pileup (generator?) is used up above, so it must be recreated
    bamfile, pileup = parse_bam(bamfile_name=bamfile_name,
                                contig_name=contig_name,
                                pysam_fasta=pysam_fasta)
    for i, column in enumerate(pileup):
        # Extract the sequence of the reference gene at the current position
        ref_base = reference_sequence[column.pos]
        # Summarise the pileup
        snv_dict, passing_snv_dict, total_coverage = \
            find_multibase_positions(ref_base=ref_base,
                                     filtered_read_dict=filtered_read_dict[i],
                                     base_cutoff=base_cutoff,
                                     base_fraction_cutoff=base_fraction_cutoff)
        # If there are any SNVs called for the gene, update the multibase_position_dict, and to_write string
        if passing_snv_dict:
            # Pysam starts counting at 0, whereas we actually want to start counting at 1.
            actual_position = column.pos + 1
            # Initialise the gene name in the dictionary as required
            if column.reference_name not in multibase_position_dict:
                multibase_position_dict[column.reference_name] = dict()
            # Update the dictionary with the actual position:
            multibase_position_dict[column.reference_name].update({actual_position: passing_snv_dict})
            to_write += position_details(actual_position=actual_position,
                                         passing_snv_dict=passing_snv_dict,
                                         contig_name=contig_name,
                                         ref_base=ref_base,
                                         total_coverage=total_coverage,
                                         base_cutoff=base_cutoff,
                                         error_perc=error_perc)
    bamfile.close()
    return multibase_position_dict, to_write


def find_rmlst_type(kma_report, rmlst_report):
    """
    Uses a report generated by KMA to determine what allele is present for each rMLST gene.
    :param kma_report: The .res report generated by KMA.
    :param rmlst_report: rMLST report file to write information to.
    :return: a sorted list of loci present, in format 'gene_allele'
    """
    genes_to_use = dict()
    score_dict = dict()
    gene_alleles = list()
    with open(kma_report) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            gene_allele = row['#Template']
            score = int(row['Score'])
            gene = gene_allele.split('_')[0]
            allele = gene_allele.split('_')[1]
            if gene not in score_dict:
                score_dict[gene] = score
                genes_to_use[gene] = allele
            else:
                if score > score_dict[gene]:
                    score_dict[gene] = score
                    genes_to_use[gene] = allele
    for gene in genes_to_use:
        gene_alleles.append(gene + '_' + genes_to_use[gene].replace(' ', ''))
    gene_alleles = sorted(gene_alleles)
    with open(rmlst_report, 'w') as f:
        f.write('Gene,Allele\n')
        for gene_allele in gene_alleles:
            gene = gene_allele.split('_')[0]
            allele = gene_allele.split('_')[1]
            f.write('{},{}\n'.format(gene, allele))
    return gene_alleles


def base_dict_to_string(base_dict):
    """
    Converts a dictionary to a string. {'C': 12, 'A':4} gets converted to C:12;A:4
    :param base_dict: Dictionary of bases and counts created by find_if_multibase
    :return: String representing that dictionary.
    """
    outstr = ''
    # First, sort base_dict so that major allele always comes first - makes output report nicer to look at.
    base_list = sorted(base_dict.items(), key=lambda kv: kv[1], reverse=True)
    for base in base_list:
        outstr += '{}:{};'.format(base[0], base[1])
    return outstr[:-1]


def find_total_sequence_length(fasta_file):
    """
    Totals up number of bases in a fasta file.
    :param fasta_file: Path to an uncompressed, fasta-formatted file.
    :return: Number of total bases in file, as an int.
    """
    total_length = 0
    for sequence in SeqIO.parse(fasta_file, 'fasta'):
        total_length += len(sequence.seq)
    return total_length


# Removing percentage contamination estimation until it can be made more 
# accurate
#-------------------------------------------------------------------------------
# def estimate_percent_contamination(contamination_report_file):
#     """
#     Estimates the percent contamination of a sample (and standard deviation).
#     :param contamination_report_file: File created by read_contig,
#     :return: Estimated percent contamination and standard deviation.
#     """
#     # Initialise a list to store the percentage of SNVs out of the total bases
#     contam_levels = list()
#     with open(contamination_report_file) as csvfile:
#         reader = csv.DictReader(csvfile)
#         for row in reader:
#             base_counts = int(row['SNVCoverage'])
#             total_coverage = int(row['TotalCoverage'])
#             contam_levels.append(base_counts * 100 / total_coverage)
#     return '%.2f' % (np.mean(contam_levels)), '%.2f' % np.std(contam_levels)


def load_fastq_records(gz, paired, forward):
    """
    Use SeqIO to load FASTQ records from file
    :param gz: Name and path of FASTQ file
    :param paired: Boolean of whether the reads are paired
    :param forward: Boolean of whether the current reads are in the forward direction
    :return: records: Dictionary of SeqIO parsed FASTQ records with consistent ID naming
    """
    # Initialise a dictionary to store the FASTQ records
    records = dict()
    # Iterate through the reads
    for record in SeqIO.parse(gz, 'fastq'):
        # Only update the naming scheme for paired reads
        if paired:
            if forward:
                # Change a :1: to /1 in the record.id
                if ':1:' in record.id:
                    record.id = record.id + '/1'
                # Don't worry if the record.id already has a /1
                elif '/1' in record.id:
                    pass
                # If the record.id doesn't have a read direction, add /1
                else:
                    record.id = record.id + '/1'
            # Process reverse reads in a similar fashion to forward reads
            else:
                if ':2:' in record.id:
                    record.id = record.id + '/2'
                elif '/2' in record.id:
                    pass
                else:
                    record.id = record.id + '/2'
        else:
            pass
        records.update(SeqIO.to_dict([record]))
    return records


def index_databases(sample_database):
    """
    Index the database file with pysam and kma
    :param sample_database:
    """
    if not os.path.isfile(sample_database + '.fai'):  # Don't bother re-indexing, this only needs to happen once.
        try:
            pysam.faidx(sample_database)
        except pysam.utils.SamtoolsError:
            pass
    kma_database = sample_database.replace('.fasta', '') + '_kma'

    if not os.path.isfile(kma_database + '.name'):  # The .name is one of the files KMA creates when making a database.
        logging.info('Since this is the first time you are using this database, it needs to be indexed by KMA. '
                     'This might take a while')
        cmd = 'kma index -i {} -o {}'.format(sample_database, kma_database)  # NOTE: Need KMA >=1.2.0 for this to work.
        out = str()
        try:
            out, err = run_cmd(cmd)
        except subprocess.CalledProcessError as e:
            err = e
        log = sample_database + '_log.txt'
        write_to_logfile(log, out, err, cmd)


# noinspection PyUnresolvedReferences
def find_contamination(pair, output_folder, databases_folder, base_cutoff, forward_id='_R1', threads=1,
                       keep_files=False,
                       quality_cutoff=20, base_fraction_cutoff=0.05, cgmlst_db=None, xmx=None,
                       tmpdir=None, data_type='Illumina', use_rmlst=False, min_matching_hashes=40,
                       fasta=False, error_cutoff=1.0, debug=False):
    """
    This needs some documentation fairly badly, so here we go.
    :param pair: This has become a misnomer. If the input reads are actually paired, needs to be a list
    with the full filepath to forward reads at index 0 and full path to reverse reads at index 1.
    If reads are unpaired, should be a list of length 1 with the only entry being the full filepath to read set.
    :param output_folder: Folder where outputs (confindr log and report, and other stuff) will be stored.
    This will be created if it does not exist. (I think - should write a test that double checks this).
    :param databases_folder: Full path to folder where ConFindr's databases live. These files can be
    downloaded from figshare in .tar.gz format (https://ndownloader.figshare.com/files/11864267), and
    will be automatically downloaded if the script is run from the command line.
    :param forward_id: Identifier that marks reads as being in the forward direction for paired reads.
    Defaults to _R1
    :param threads: Number of threads to run analyses with. All parts of this pipeline scale pretty well,
    so more is better.
    :param keep_files: Boolean that says whether or not to keep temporary files.
    :param quality_cutoff: Integer of the phred score required to have a base count towards a multiallelic site.
    :param base_cutoff: Integer of number of bases needed to have a base be part of a multiallelic site.
    :param base_fraction_cutoff: Float of fraction of bases needed to have a base be part of a multiallelic site.
    If specified will be used in parallel with base_cutoff
    :param cgmlst_db: if None, we're using rMLST, if True, using some sort of custom cgMLST database. This requires some
    custom parameters.
    :param xmx: if None, BBTools will use auto memory detection. If string, BBTools will use what's specified as their
    memory request.
    :param tmpdir: if None, any genus-specifc databases that need to be created will be written to ConFindr DB location.
    :param data_type: Either Illumina or Nanopore, depending on what type your reads are. (STR)
    :param use_rmlst: If False, use cgderived data instead of rMLST where possible. If True, always use rMLST. (BOOL)
    :param min_matching_hashes: Minimum number of matching hashes in a MASH screen in order for a genus to be
    considered present in a sample. Default is 40
    :param fasta: Boolean on whether the samples are in FASTA format. Default is False
    :param error_cutoff: Float of the error cutoff value to use. Default is 1.0%
    :param debug: Run the find_contamination with multi-processing to allow easier debugging
    """
    if os.path.isfile(os.path.join(databases_folder, 'download_date.txt')):
        with open(os.path.join(databases_folder, 'download_date.txt')) as f:
            database_download_date = f.readline().rstrip()
    else:
        database_download_date = 'ND'
    log = os.path.join(output_folder, 'confindr_log.txt')
    if len(pair) == 2:
        sample_name = os.path.split(pair[0])[-1].split(forward_id)[0]
        paired = True
        logging.debug('Sample is paired. Sample name is {}'.format(sample_name))
    else:
        sample_name = os.path.split(pair[0])[-1].split('.')[0]
        paired = False
        logging.debug('Sample is unpaired. Sample name is {}'.format(sample_name))
    sample_tmp_dir = os.path.join(output_folder, sample_name)
    if not os.path.isdir(sample_tmp_dir):
        os.makedirs(sample_tmp_dir)

    logging.info('Checking for cross-species contamination...')
    if paired:
        genus = find_cross_contamination(databases_folder,
                                         reads=pair,
                                         sample_name=sample_name,
                                         tmpdir=sample_tmp_dir,
                                         log=log,
                                         threads=threads,
                                         min_matching_hashes=min_matching_hashes)
    else:
        genus = find_cross_contamination(databases_folder,
                                         reads=pair[0],
                                         sample_name=sample_name,
                                         tmpdir=sample_tmp_dir,
                                         log=log,
                                         threads=threads,
                                         min_matching_hashes=min_matching_hashes)
    # if len(genus.split(':')) > 1:
    #     write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
    #                     sample_name=sample_name,
    #                     multi_positions=0,
    #                     genus=genus,
    #                     total_gene_length=0,
    #                     database_download_date=database_download_date)
    #     logging.info('Found cross-contamination! Skipping rest of analysis...\n')
    #     if keep_files is False:
    #         shutil.rmtree(sample_tmp_dir)
    #     return
    # Setup genus-specific databases, if necessary.
    if cgmlst_db is not None:
        # Sanity check that the DB specified is actually a file, otherwise, quit with appropriate error message.
        if not os.path.isfile(cgmlst_db):
            logging.error('ERROR: Specified cgMLST file ({}) does not exist. Please check the path and try again.'
                          .format(cgmlst_db))
            quit(code=1)
        sample_database = cgmlst_db
    else:
        db_folder = databases_folder if tmpdir is None else tmpdir
        if not os.path.isdir(db_folder):
            os.makedirs(db_folder)
        if genus != 'ND':
            # Logic here is as follows: users can either have both rMLST databases, which cover all of bacteria,
            # cgmlst-derived databases, which cover only Escherichia, Salmonella, and Listeria (may add more at some
            # point), or they can have both. They can also set priority to either always use rMLST, or to use my
            # core-genome derived stuff and fall back on rMLST if they're trying to look at a genus I haven't created
            # a scheme for.
            #
            if len(genus.split(':')) > 1:
                predominant_genus = genus.split(':')[0]
            else:
                predominant_genus = genus
            # In the event rmlst databases have priority, always use them.
            if use_rmlst is True:
                sample_database = os.path.join(db_folder, '{}_db.fasta'.format(predominant_genus))
                if not os.path.isfile(sample_database):

                    if os.path.isfile(os.path.join(db_folder, 'gene_allele.txt')) and \
                            os.path.isfile(os.path.join(db_folder, 'rMLST_combined.fasta')):
                        logging.info('Setting up rMLST genus-specific database for genus {}...'
                                     .format(predominant_genus))
                        allele_list = find_genus_specific_allele_list(os.path.join(db_folder, 'gene_allele.txt'),
                                                                      predominant_genus)
                        # Create the allele-specific database
                        setup_allelespecific_database(fasta_file=sample_database,
                                                      database_folder=db_folder,
                                                      allele_list=allele_list)
            else:
                # Check if a cgderived database is available. If not, try to use rMLST database.
                sample_database = os.path.join(db_folder, '{}_db_cgderived.fasta'.format(predominant_genus))
                if not os.path.isfile(sample_database):
                    sample_database = os.path.join(db_folder, '{}_db.fasta'.format(predominant_genus))
                    # Create genus specific database if it doesn't already exist and we have the necessary rMLST files.
                    if os.path.isfile(os.path.join(db_folder, 'rMLST_combined.fasta')) and \
                            os.path.isfile(os.path.join(db_folder, 'gene_allele.txt')) and not \
                            os.path.isfile(sample_database):
                        logging.info('Setting up core genome genus-specific database for genus {}...'
                                     .format(predominant_genus))
                        allele_list = find_genus_specific_allele_list(os.path.join(db_folder, 'gene_allele.txt'),
                                                                      predominant_genus)
                        setup_allelespecific_database(fasta_file=sample_database,
                                                      database_folder=db_folder,
                                                      allele_list=allele_list)

        else:
            sample_database = os.path.join(db_folder, 'rMLST_combined.fasta')

    # If a user has gotten to this point and they don't have any database available to do analysis because
    # they don't have rMLST downloaded and we don't have a cg-derived database available, boot them with a helpful
    # message.
    if not os.path.isfile(sample_database):
        write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
                     sample_name=sample_name,
                     multi_positions=0,
                     genus=genus,
                     total_gene_length=0,
                     database_download_date=database_download_date)
        logging.info('Did not find databases for genus {genus}. You can download the rMLST database to get access to '
                     'all genera (see https://olc-bioinformatics.github.io/ConFindr/install/). Alternatively, if you '
                     'have a high-quality core-genome derived database for your genome of interest, we would be happy '
                     'to add it - open an issue at https://github.com/OLC-Bioinformatics/ConFindr/issues with the '
                     'title "Add genus-specific database: {genus}"\n'.format(genus=genus))
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
    if paired:
        forward_bait = os.path.join(sample_tmp_dir, '{sn}_baited_R1.fastq.gz'.format(sn=sample_name))
        reverse_bait = forward_bait.replace('_R1', '_R2')
        if not os.path.isfile(forward_bait):
            if xmx is None:
                out, err, cmd = bbtools.bbduk_bait(reference=sample_database,
                                                   forward_in=pair[0],
                                                   reverse_in=pair[1],
                                                   forward_out=forward_bait,
                                                   reverse_out=reverse_bait,
                                                   threads=threads,
                                                   returncmd=True)
            else:
                out, err, cmd = bbtools.bbduk_bait(reference=sample_database,
                                                   forward_in=pair[0],
                                                   reverse_in=pair[1],
                                                   forward_out=forward_bait,
                                                   reverse_out=reverse_bait,
                                                   threads=threads,
                                                   Xmx=xmx,
                                                   returncmd=True)
    else:
        # Still name the file '_baited_trimmed' even if the file won't be trimmed
        if data_type == 'Nanopore' or fasta:
            unpaired_bait = os.path.join(sample_tmp_dir, '{sn}_baited_trimmed.fastq.gz'.format(sn=sample_name))
        else:
            unpaired_bait = os.path.join(sample_tmp_dir, '{sn}_baited.fastq.gz'.format(sn=sample_name))
        if not os.path.isfile(unpaired_bait):
            if xmx is None:
                out, err, cmd = bbtools.bbduk_bait(reference=sample_database,
                                                   forward_in=pair[0],
                                                   forward_out=unpaired_bait,
                                                   returncmd=True, threads=threads)
            else:
                out, err, cmd = bbtools.bbduk_bait(reference=sample_database,
                                                   forward_in=pair[0],
                                                   forward_out=unpaired_bait,
                                                   Xmx=xmx,
                                                   returncmd=True, threads=threads)
    if out:
        write_to_logfile(log, out, err, cmd)
    logging.info('Quality trimming...')
    out, err, cmd = '', '', ''
    if data_type == 'Illumina':
        if paired:
            forward_trimmed = os.path.join(sample_tmp_dir, '{sn}_baited_trimmed_R1.fastq.gz'.format(sn=sample_name))
            reverse_trimmed = forward_trimmed.replace('_R1', '_R2')
            if not os.path.isfile(forward_trimmed):
                if xmx is None:
                    out, err, cmd = bbtools.bbduk_trim(forward_in=forward_bait,
                                                       reverse_in=reverse_bait,
                                                       forward_out=forward_trimmed,
                                                       reverse_out=reverse_trimmed,
                                                       threads=str(threads), returncmd=True)
                else:
                    out, err, cmd = bbtools.bbduk_trim(forward_in=forward_bait,
                                                       reverse_in=reverse_bait,
                                                       forward_out=forward_trimmed,
                                                       reverse_out=reverse_trimmed,
                                                       Xmx=xmx,
                                                       threads=str(threads),
                                                       returncmd=True)

            with gzip.open(forward_trimmed, 'rt') as gz:
                fastq_records = load_fastq_records(gz=gz,
                                                   paired=True,
                                                   forward=True)
            with gzip.open(reverse_trimmed, 'rt') as gz:
                # fastq_records.update(SeqIO.to_dict(SeqIO.parse(gz, 'fastq')))
                fastq_records.update(load_fastq_records(gz=gz,
                                                        paired=True,
                                                        forward=False))

        else:
            unpaired_trimmed = os.path.join(sample_tmp_dir, '{sn}_baited_trimmed.fastq.gz'.format(sn=sample_name))

            if not fasta:
                if not os.path.isfile(unpaired_trimmed):
                    if xmx is None:
                        out, err, cmd = bbtools.bbduk_trim(forward_in=unpaired_bait,
                                                           forward_out=unpaired_trimmed,
                                                           returncmd=True,
                                                           threads=threads)
                    else:
                        out, err, cmd = bbtools.bbduk_trim(forward_in=unpaired_bait,
                                                           forward_out=unpaired_trimmed,
                                                           returncmd=True,
                                                           threads=threads,
                                                           Xmx=xmx)
                with gzip.open(unpaired_trimmed, 'rt') as gz:
                    # fastq_records = SeqIO.to_dict(SeqIO.parse(gz, 'fastq'))
                    fastq_records = load_fastq_records(gz=gz,
                                                       paired=False,
                                                       forward=True)
            else:
                # '''unpaired_bait'''
                with gzip.open(unpaired_bait, 'rt') as gz:
                    # fastq_records = SeqIO.to_dict(SeqIO.parse(gz, 'fastq'))
                    fastq_records = load_fastq_records(gz=gz,
                                                       paired=False,
                                                       forward=True)
        write_to_logfile(log, out, err, cmd)
    else:
        if paired:
            with gzip.open(forward_bait, 'rt') as gz:
                # fastq_records = SeqIO.to_dict(SeqIO.parse(gz, 'fastq'))
                fastq_records = load_fastq_records(gz=gz,
                                                   paired=True,
                                                   forward=True)
            with gzip.open(reverse_bait, 'rt') as gz:
                # fastq_records.update(SeqIO.to_dict(SeqIO.parse(gz, 'fastq')))
                fastq_records.update(load_fastq_records(gz=gz,
                                                        paired=True,
                                                        forward=False))
        else:
            with gzip.open(unpaired_bait, 'rt') as gz:
                # fastq_records = SeqIO.to_dict(SeqIO.parse(gz, 'fastq'))
                fastq_records = load_fastq_records(gz=gz,
                                                   paired=False,
                                                   forward=True)
    logging.info('Detecting contamination...')
    # Now do mapping in two steps - first, map reads back to database with ambiguous reads matching all - this
    # will be used to get a count of number of reads aligned to each gene/allele so we can create a custom rmlst file
    # with only the most likely allele for each gene.
    kma_report = os.path.join(sample_tmp_dir, '{sn}_kma'.format(sn=sample_name))
    # if not os.path.isfile(sample_database + '.fai'):  # Don't bother re-indexing, this only needs to happen once.
    #     pysam.faidx(sample_database)
    kma_database = sample_database.replace('.fasta', '') + '_kma'
    #
    # if not os.path.isfile(kma_database + '.name'):  # The .name is one of the files KMA creates when making a database
    #     logging.info('Since this is the first time you are using this database, it needs to be indexed by KMA. '
    #                  'This might take a while')
    #     cmd = 'kma index -i {} -o {}'.format(sample_database, kma_database)  # NOTE: Need KMA >=1.2.0 for this to work
    #     out, err = run_cmd(cmd)
    #     write_to_logfile(log, out, err, cmd)
    index_databases(sample_database=sample_database)
    # Run KMA.
    if paired:
        if not os.path.isfile(kma_report + '.res'):
            cmd = 'kma -ipe {forward_in} {reverse_in} -t_db {kma_database} -o {kma_report} ' \
                  '-t {threads}'.format(forward_in=forward_trimmed,
                                        reverse_in=reverse_trimmed,
                                        kma_database=kma_database,
                                        kma_report=kma_report,
                                        threads=threads)
            out, err = run_cmd(cmd)
            write_to_logfile(log, out, err, cmd)
    else:
        if not os.path.isfile(kma_report + '.res'):
            if data_type == 'Illumina':
                # Use the FASTA file (rather than the reads) as the input
                if fasta:
                    cmd = 'kma -i {input_reads} -t_db {kma_database} -mem_mode -ID 100 -ConClave 2 -ex_mode ' \
                          '-o {kma_report} -t {threads}' \
                        .format(input_reads=pair[0],
                                kma_database=kma_database,
                                kma_report=kma_report,
                                threads=threads)
                else:
                    cmd = 'kma -i {input_reads} -t_db {kma_database} -o {kma_report} ' \
                          '-t {threads}'.format(input_reads=unpaired_trimmed,
                                                kma_database=kma_database,
                                                kma_report=kma_report,
                                                threads=threads)
            else:
                # Recommended Nanopore settings from KMA repo: https://bitbucket.org/genomicepidemiology/kma
                cmd = 'kma -i {input_reads} -t_db {kma_database} -o {kma_report} -mem_mode -mp 20 -mrs 0.0 -bcNano ' \
                      '-t {threads}'.format(input_reads=unpaired_bait,
                                            kma_database=kma_database,
                                            kma_report=kma_report,
                                            threads=threads)
            out, err = run_cmd(cmd)
            write_to_logfile(log, out, err, cmd)

    rmlst_report = os.path.join(output_folder, sample_name + '_alleles.csv')
    gene_alleles = find_rmlst_type(kma_report=kma_report + '.res',
                                   rmlst_report=rmlst_report)

    rmlst_fasta = os.path.join(sample_tmp_dir, '{sn}_alleles.fasta'.format(sn=sample_name))
    if not os.path.isfile(rmlst_fasta):
        with open(rmlst_fasta, 'w') as f:
            for contig in SeqIO.parse(sample_database, 'fasta'):
                if contig.id in gene_alleles:
                    SeqIO.write(contig, f, 'fasta')
                    # f.write('>{}\n'.format(contig.id))
                    # f.write(str(contig.seq) + '\n')
    rmlst_gene_length = find_total_sequence_length(rmlst_fasta)
    logging.debug('Total gene length is {}'.format(rmlst_gene_length))
    pysam_pass = True
    # Second step of mapping - Do a mapping of our baited reads against a fasta file that has only one allele per
    # rMLST gene.
    try:
        if not os.path.isfile(rmlst_fasta + '.fai'):
            pysam.faidx(rmlst_fasta)
        outbam = os.path.join(sample_tmp_dir, '{sn}_contamination.bam'.format(sn=sample_name))
        sorted_bam = os.path.join(sample_tmp_dir, '{sn}_contamination_sorted.bam'.format(sn=sample_name))
        outsam = outbam.replace('.bam', '.sam')
        if not os.path.isfile(sorted_bam):
            if paired:
                cmd = 'bbmap.sh ref={ref} in={forward_in} in2={reverse_in} out={outbam} threads={threads} mdtag ' \
                      'nodisk'.format(ref=rmlst_fasta,
                                      forward_in=forward_trimmed,
                                      reverse_in=reverse_trimmed,
                                      outbam=outbam,
                                      threads=threads)
                if cgmlst_db is not None:
                    # Lots of core genes seem to have relatives within a genome that are at ~70% identity. This means
                    # that reads that shouldn't map do, and cause false positives. Adding in this sub-filter means that
                    # reads can only have one mismatch, so they have to be from the right gene for this to work.
                    cmd += ' subfilter=1'
                if xmx:
                    cmd += ' -Xmx{}'.format(xmx)
                out, err = run_cmd(cmd)
                write_to_logfile(log, out, err, cmd)
            else:
                #
                if data_type == 'Illumina' and not fasta:
                    cmd = 'bbmap.sh ref={ref} in={forward_in} out={outbam} threads={threads} mdtag ' \
                          'nodisk'.format(ref=rmlst_fasta,
                                          forward_in=unpaired_trimmed,
                                          outbam=outbam,
                                          threads=threads)
                    if cgmlst_db is not None:
                        # Core genes can have relatives within a genome that are at ~70 percent identity. This means
                        # that reads that shouldn't map do, and cause false positives. Adding in this sub-filter means
                        # reads can only have one mismatch, so they have to be from the right gene for this to work.
                        cmd += ' subfilter=1'
                    if xmx:
                        cmd += ' -Xmx{}'.format(xmx)
                    out, err = run_cmd(cmd)
                    write_to_logfile(log, out, err, cmd)
                else:
                    if fasta:
                        ax = 'asm5'
                        # ax = 'sr'
                    # If Nanopore FASTQ:
                    else:
                        ax = 'map-ont'
                    # '> {outsam}'
                    cmd = 'minimap2 --MD -t {threads} -ax {ax} {ref} {reads}' \
                        .format(ax=ax,
                                ref=rmlst_fasta,
                                reads=unpaired_bait,
                                outsam=outsam,
                                threads=threads)
                    # ' | samtools rmdup - -S -' \
                    cmd += ' | samtools view -@ {threads} -h -bT {abs_ref_link} -' \
                           ' | samtools sort - -@ {threads} -o {sorted_bam}' \
                        .format(threads=threads,
                                abs_ref_link=rmlst_fasta,
                                sorted_bam=sorted_bam)
                    out, err = run_cmd(cmd)
                    write_to_logfile(log, out, err, cmd)

        if not os.path.isfile(sorted_bam):
            # noinspection PyUnresolvedReferences
            pysam.sort('-o', sorted_bam, outbam)
        if not os.path.isfile(sorted_bam + '.bai'):
            pysam.index(sorted_bam)
        # Now find number of multi-positions for each rMLST gene/allele combination
        multi_positions = 0

        # Run the BAM parsing in parallel! Some refactoring of the code would likely be a good idea so this
        # isn't quite so ugly, but it works.
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
        multibase_dict_list = list()
        report_write_list = list()
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
                    nanopore = nanopore_list[i],
                    error_cutoff=error_cutoff_list[i])
                multibase_dict_list.append(multibase_dict)
                report_write_list.append(report_write)
        else:
            for multibase_dict, report_write in p.starmap(read_contig,
                                                          zip(gene_alleles,
                                                              bamfile_list,
                                                              reference_fasta_list,
                                                              records_list,
                                                              fastq_records_list,
                                                              quality_cutoff_list,
                                                              base_cutoff_list,
                                                              base_fraction_list,
                                                              fasta_list,
                                                              error_cutoff_list,
                                                              nanopore_list),
                                                          chunksize=1):
                multibase_dict_list.append(multibase_dict)
                report_write_list.append(report_write)
            p.close()
            p.join()
    except SamtoolsError:
        pysam_pass = False
        multi_positions = 0
        multibase_dict_list = list()
        report_write_list = list()
    # Write out report info.
    report_file = os.path.join(output_folder, sample_name + '_contamination.csv')
    with open(report_file, 'w') as r:
        r.write('Gene,Position,RefBase,CongruentSNVs,TotalSNVs,ForwardSNVs,ReverseSNVs,SNVCoverage,TotalCoverage,'
                'BaseCutoff,ErrorPercent\n')
        for item in report_write_list:
            for contamination_info in item:
                r.write(contamination_info)
    # Total up the number of multibase positions.
    for multibase_position_dict in multibase_dict_list:
        multi_positions += sum([len(snp_positions) for snp_positions in multibase_position_dict.values()])
    if cgmlst_db is None:
        snp_cutoff = math.ceil(rmlst_gene_length / 10000) + 1
    elif fasta:
        snp_cutoff = 1
    else:
        snp_cutoff = 10

    logging.info('Done! Number of contaminating SNVs found: {}\n'.format(multi_positions))
    write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
                 sample_name=sample_name,
                 multi_positions=multi_positions,
                 genus=genus,
                 total_gene_length=rmlst_gene_length,
                 snp_cutoff=snp_cutoff,
                 database_download_date=database_download_date,
                 pysam_pass=pysam_pass)
    if keep_files is False:
        shutil.rmtree(sample_tmp_dir)


def write_output(output_report, sample_name, multi_positions, genus, total_gene_length,
                 database_download_date, snp_cutoff=3, pysam_pass=True):
    """
    Function that writes the output generated by ConFindr to a report file. Appends to a file that already exists,
    or creates the file if it doesn't already exist.
    :param output_report: Path to CSV output report file. Should have headers SampleName,Genus,NumContamSNVs, and
    ContamStatus, in that order.
    :param sample_name: string - name of sample
    :param multi_positions: integer - number of positions that were found to have more than one base present.
    :param genus: string - The genus of your sample
    :param total_gene_length: integer - number of bases examined to make a contamination call.
    :param database_download_date:
    :param snp_cutoff: Number of cSNVs to use to call a sample contaminated. Default 3. (INT)
    :param pysam_pass: Boolean of whether pysam encountered an error
    """
    # If the report file hasn't been created, make it, with appropriate header.
    if not os.path.isfile(output_report):
        with open(os.path.join(output_report), 'w') as f:
            f.write('Sample,Genus,NumContamSNVs,ContamStatus,'
                    'BasesExamined,DatabaseDownloadDate\n')
    if pysam_pass:
        if multi_positions >= snp_cutoff or len(genus.split(':')) > 1:
            contaminated = True
        else:
            contaminated = False
    else:
        contaminated = 'Pysam SamtoolsError'
        multi_positions = 'ND'
    with open(output_report, 'a+') as f:
        f.write('{samplename},{genus},{numcontamsnvs},{contamstatus},'
                '{gene_length},{database_download_date}\n'.format(samplename=sample_name,
                                                                  genus=genus,
                                                                  numcontamsnvs=multi_positions,
                                                                  contamstatus=contaminated,
                                                                  gene_length=total_gene_length,
                                                                  database_download_date=database_download_date))


def check_for_databases_and_download(database_location):
    # Check for the files necessary - should have rMLST_combined.fasta, gene_allele.txt, profiles.txt, and refseq.msh
    necessary_files = ['Escherichia_db_cgderived.fasta', 'Listeria_db_cgderived.fasta',
                       'Salmonella_db_cgderived.fasta', 'refseq.msh']
    optional_files = ['rMLST_combined.fasta', 'gene_allele.txt', 'profiles.txt']
    all_files_present = True
    for necessary_file in necessary_files:
        if not os.path.isfile(os.path.join(database_location, necessary_file)):
            logging.warning('Could not find {}'.format(necessary_file))
            all_files_present = False

    if not all_files_present:
        logging.warning('Databases not present - downloading basic databases now...')
        if not os.path.isdir(database_location):
            os.makedirs(database_location)
        download_mash_sketch(database_location)
        download_cgmlst_derived_data(database_location)

    optional_files_present = True
    for optional_file in optional_files:
        if not os.path.isfile(os.path.join(database_location, optional_file)):
            optional_files_present = False
    if not optional_files_present:
        logging.warning('Did not find rMLST databases, if you want to use ConFindr on genera other than Listeria, '
                        'Salmonella, and Escherichia, you\'ll need to download them. Instructions are available at '
                        'https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases\n')


def check_valid_base_fraction(base_fraction):
    """
    Checks that base fraction specified is either None, and therefore won't be used, or is between 0 and 1.
    :param base_fraction: Base fraction, which should be either None or a float.
    :return: True if base fraction is valid, False if not.
    """
    if base_fraction is None:
        return True
    if 0 <= base_fraction <= 1:
        return True
    else:
        return False


def check_acceptable_xmx(xmx_string):
    """
    BBTools can have their memory set manually. This will check that the memory setting is actually valid
    :param xmx_string: The users requested XMX, as a string.
    :return: True if the Xmx string will be accepted by BBTools, otherwise false.
    """
    acceptable_xmx = True
    acceptable_suffixes = ['K', 'M', 'G']
    if xmx_string[-1].upper() not in acceptable_suffixes:
        acceptable_xmx = False
        logging.error('ERROR: Memory must be specified as K (kilobytes), M (megabytes), or G (gigabytes). '
                      'Your specified suffix was {}.'.format(xmx_string[-1]))
    if '.' in xmx_string:
        acceptable_xmx = False
        logging.error('ERROR: Xmx strings must be integers, floating point numbers are not accepted.')
    if not str.isdigit(xmx_string[:-1]):
        acceptable_xmx = False
        logging.error('ERROR: The amount of memory requested was not an integer.')
    return acceptable_xmx


def get_version():
    try:
        version = 'ConFindr {}'.format(pkg_resources.get_distribution('confindr').version)
    except pkg_resources.DistributionNotFound:
        version = 'ConFindr (Unknown version)'
    return version
