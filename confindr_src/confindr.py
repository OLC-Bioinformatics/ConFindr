#!/usr/bin/env python
from pysam.utils import SamtoolsError
import multiprocessing
import pkg_resources
import numpy as np
import subprocess
import traceback
import argparse
import logging
import shutil
import glob
import gzip
import csv
import os
import pysam
from Bio import SeqIO
from confindr_src.database_setup import download_cgmlst_derived_data, download_mash_sketch
from confindr_src.wrappers import mash
from confindr_src.wrappers import bbtools


def run_cmd(cmd):
    """
    Runs a command using subprocess, and returns both the stdout and stderr from that command
    If exit code from command is non-zero, raises subproess.CalledProcessError
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


def find_genusspecific_allele_list(profiles_file, target_genus):
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
    :param allele_list: allele list generated by find_genusspecific_allele_list
    """
    index = SeqIO.index(os.path.join(database_folder, 'rMLST_combined.fasta'), 'fasta')
    seqs = list()
    for s in allele_list:
        try:
            seqs.append(index[s])
        except KeyError:
            logging.warning('Tried to add {} to allele-specific database, but could not find it.'.format(s))
    SeqIO.write(seqs, fasta_file, 'fasta')


def extract_rmlst_genes(pair, database, forward_out, reverse_out, threads=12, logfile=None):
    """
    Given a pair of reads and an rMLST database, will extract reads that contain sequence from the database.
    :param pair: List containing path to forward reads at index 0 and path to reverse reads at index 1.
    :param database: Path to rMLST database, in FASTA format.
    :param forward_out: Filepath to write forward reads to (STR)
    :param reverse_out: Filepath to write reverse reads to (STR)
    :param threads: Number of threads to use (INT)
    :param logfile: Logfile to write command/stdout/stderr to. If None, nothing will be written.
    """
    out, err, cmd = bbtools.bbduk_bait(database, pair[0], forward_out, reverse_in=pair[1],
                                       reverse_out=reverse_out, threads=str(threads), returncmd=True)
    if logfile:
        write_to_logfile(logfile, out, err, cmd)


def find_cross_contamination(databases, reads, sample_name, tmpdir='tmp', log='log.txt', threads=1,
                             min_matching_hashes=40):
    """
    Uses mash to find out whether or not a sample has more than one genus present, indicating cross-contamination.
    :param databases: A databases folder, which must contain refseq.msh, a mash sketch that has one representative
    per genus from refseq.
    :param reads: Relative path(s) to either unpaired (type STR) or paired (type LIST) FASTQ reads
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
        if mash_genus == 'Shigella':
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
    Finds if a site has at least two bases of  high quality, enough that it can be considered
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
        bases_above_threshold = {base: float(count) / total_hq_base_count >= base_fraction_cutoff and
                                       count >= base_count_cutoff for (base, count) in high_quality_base_count.items()}
    else:
        bases_above_threshold = {base: count >= base_count_cutoff for (base, count) in high_quality_base_count.items()}

    # True is equal to 1 so sum of the number of Trues in the bases_above_threshold dict is the number of bases
    # passing threshold
    total = sum(bases_above_threshold.values())
    return sum(bases_above_threshold.values())


def characterise_read(column, reference_sequence, fastq_records, quality_cutoff, fasta=False):
    """
    Finds if a position in a pileup has more than one base present.
    :param column: A pileupColumn generated by pysam
    :param quality_cutoff: Desired min phred quality for a base in order to be counted towards a multi-allelic column
    :param base_cutoff: Minimum number of bases needed to support presence of a base.
    :param base_fraction_cutoff: Minimum fraction of bases needed to support presence of a base.
    If specified, noth the base_cutoff and base_fraction_cutoff will have to be met
    :return: If position has more than one base, a dictionary with counts for the bases. Otherwise, returns
    empty dictionary
    """
    # Sometimes the qualities come out to ridiculously high (>70) values. Looks to be because sometimes reads
    # are overlapping and the qualities get summed for overlapping bases. Issue opened on pysam.
    # unfiltered_base_qualities = dict()
    read_dict = {
        'forward': dict(),
        'reverse': dict(),
        'unmapped': dict()
    }
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
    # print('!')
    # quit()
    unfiltered_read_details = dict()
    read_list = list()
    snv_dict = dict()
    for read in column.pileups:
        if read.query_position is not None:  # Not entirely sure why this is sometimes None, but it causes bad stuff

            if read.alignment.qname not in unfiltered_read_details:
                unfiltered_read_details[read.alignment.qname] = dict()
            # reference_sequence = read.alignment.get_reference_sequence()
            previous_position = column.pos - 1 if column.pos >= 1 else 0
            # This causes index errors. Fix at some point soon. -- FIXED?
            next_position = column.pos + 1 if column.pos < len(reference_sequence) - 1 else len(reference_sequence) - 1
            previous_positions = column.pos - 5 if column.pos >= 5 else 0
            # This causes index errors. Fix at some point soon. -- FIXED?
            next_positions = column.pos + 5 if column.pos < len(reference_sequence) - 5 else column.pos + len(
                reference_sequence) - column.pos
            prev_base_range = range(previous_positions, column.pos)
            next_base_range = range(next_position, next_positions)
            # print(column.reference_name, read.alignment.qname, column.pos, previous_position, next_position)
            # quit()
            from itertools import chain
            # Another stringency check - to make sure that we're actually looking at a point mutation, check that the
            # base before and after the one we're looking at match the reference. With Nanopore data, lots of indels and
            # the like cause false positives, so this filters those out.
            # try:  # Need to actually handle this at some point. For now, be lazy
            # print(read.alignment.qname)
            # previous_reference_base = reference_sequence[previous_position]
            # next_reference_base = reference_sequence[next_position]
            # previous_base = read.alignment.query_sequence[read.query_position - 1 if read.query_position > 1 else 0]
            # next_base = \
            #     read.alignment.query_sequence[read.query_position + 1
            #                                   if read.query_position < read.alignment.query_alignment_end - 1
            #                                   else read.alignment.query_alignment_end - 1]
            query_base = read.alignment.query_sequence[read.query_position]
            previous_read_position = read.query_position - 1 if read.query_position >= 1 else 0
            previous_read_positions = read.query_position - 5 if read.query_position >= 5 else 0
            next_read_position = read.query_position + 1 \
                if read.query_position < read.alignment.query_alignment_end - 1 \
                else read.alignment.query_alignment_end - 1
            next_read_positions = read.query_position + 5 \
                if read.query_position < read.alignment.query_alignment_end - 5 \
                else read.query_position + read.alignment.query_alignment_end - read.query_position

            # base = read.alignment.query_alignment_sequence[read.query_position]
            # quality = read.alignment.query_qualities[column.pos]
            # try:
            if not fasta:
                read_name = read.alignment.qname.split(' ')[0] + '/1' if read.alignment.is_read1 else \
                    read.alignment.qname.split(' ')[0] + '/2'
            else:
                read_name = read.alignment.qname
            quality = fastq_records[read_name].letter_annotations["phred_quality"][read.query_position]
            # raise
            # print(column.reference_name, read.alignment.qname, column.pos, previous_position, next_position)
            # quit()
            # range_dict = dict()
            # for i, contig_pos in enumerate(chain(prev_base_range, next_base_range)):
            #     try:
            #         read_pos = list(chain(range(previous_read_positions, read.query_position), range(next_read_position, next_read_positions)))[i]
            #     except IndexError:
            #         read_pos = int()
            #     if read_pos:
            #         range_dict[contig_pos] = read_pos
            # print(i, contig_pos, j)
            range_dict = dict()
            for iterator, contig_pos in enumerate(
                    chain(range(column.pos - 5, column.pos), range(column.pos + 1, column.pos + 6))):

                if 0 <= contig_pos < len(reference_sequence) - 1:
                    read_pos = list(chain(range(read.query_position - 5, read.query_position),
                                          range(read.query_position + 1, read.query_position + 6)))[iterator]
                    if 0 <= read_pos < read.alignment.query_alignment_end - 1:
                        range_dict[contig_pos] = read_pos

            # NC_013353.1-1216812
            # NC_013353.1-1124022
            # NC_002695.1-290602
            # if read.alignment.is_reverse:
            # ref_seq = reference_sequence[read.query_position]
            ref_base = reference_sequence[column.pos]
            match = query_base == ref_base
            #
            add_base = True
            for contig_pos, read_pos in range_dict.items():
                try:
                    reference_base = reference_sequence[contig_pos]
                    base = read.alignment.query_sequence[read_pos]

                    if not match and reference_base != base:
                        # print(column.reference_name, read_name)
                        # print(ref_base, query_base)
                        # print(column.pos, contig_pos, read_pos, reference_base, base)
                        # print(range_dict)
                        # print(prev_base_range)
                        # print(previous_reference_base, previous_base)
                        # print(next_reference_base, next_base)
                        # print(next_base_range)
                        add_base = False
                except KeyError:
                    pass
            # if not add_base:
            # print(column.reference_name)
            # print(column.pos, previous_position, previous_positions, next_position, next_positions)
            # print(prev_base_range)
            # print(next_base_range)
            # print(column.reference_name, read.alignment.qname, column.pos, previous_position, next_position)
            # quit()
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
            # else:
            # print(column.reference_name, read.alignment.qname, column.pos, previous_position, next_position)
            # # # else:
            # # #     ref_seq = reference_sequence[column.pos]
            # #
            # match = base == ref_seq
            # if not match:
            #     ref_seq = reference_sequence[column.pos]
            # match = base == ref_seq
            # if read.alignment.mate_is_unmapped:
            #     if query_base not in read_dict['unmapped']:
            #         read_dict['unmapped'][query_base] = [quality]
            #     else:
            #         read_dict['unmapped'][query_base].append(quality)
            # else:
            #     if read.alignment.is_read1:
            #         if query_base not in read_dict['forward']:
            #             read_dict['forward'][query_base] = [quality]
            #         else:
            #             read_dict['forward'][query_base].append(quality)
            #         if not match:
            #             if read.alignment.qname not in snv_dict:
            #                 snv_dict[read.alignment.qname].append('forward')
            #     else:
            #         if query_base not in read_dict['reverse']:
            #             read_dict['reverse'][query_base] = [quality]
            #         else:
            #             read_dict['reverse'][query_base].append(quality)
            #         if not match:
            #             if read.alignment.qname not in snv_dict:
            #                 snv_dict[read.alignment.qname].append('reverse')
            # if not match:
            # # #     # pass
            # # # # if read.alignment.qname == 'NC_013353.1-1057328':
            # #  and read.alignment.qname == 'NC_002695.2-935545'
            #     if column.reference_name == 'b0014_11' and read.alignment.qname == 'NC_002695.2-965583':
            #         read_list.append(
            #             ('col pos {:03d}'.format(column.pos),
            #              'qp {:03d}'.format(read.query_position),
            #              'ref ' + column.reference_name,
            #              'qname ' + read.alignment.qname,
            #              'paired ' + str(read.alignment.is_paired),
            #              # 'proper pair? ' + str(read.alignment.is_proper_pair),
            #              'rev ' + str(read.alignment.is_reverse),
            #              'read1 ' + str(read.alignment.is_read1),
            #              'read2 ' + str(read.alignment.is_read2),
            #              'qseq ' + query_base,
            #              'rseq ' + ref_base,
            #              'match ' + str(match),
            #              'qual {:02d}'.format(quality),
            #              # str(reference_sequence)[:5],
            #              'mate unmap ' + str(read.alignment.mate_is_unmapped)
            #              )
            #         )
            # # if previous_reference_base == previous_base and next_reference_base == next_base:
            # if query_base not in unfiltered_base_qualities:
            #     unfiltered_base_qualities[query_base] = [quality]
            # else:
            #     unfiltered_base_qualities[query_base].append(quality)
            # except IndexError:
            #     print('ewje')
    # if read_list:
    #     for reado in read_list:
    #         print(reado)
    # print(len(read_list))

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
                    if dir_dict[True]['qual'] > quality_cutoff and dir_dict[False]['qual'] > quality_cutoff:
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
                    elif dir_dict[True]['qual'] > quality_cutoff:
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
                    elif dir_dict[False]['qual'] > quality_cutoff:
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
                if dir_dict[True]['qual'] > quality_cutoff:
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
                if dir_dict[False]['qual'] > quality_cutoff:
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
                    if dir_dict[direction]['qual'] > quality_cutoff:
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
                    if dir_dict[direction]['qual'] > quality_cutoff:
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
    return filtered_read_dict


# def get_contig_names(fasta_file):
#     """
#     Gets contig names from a fasta file using SeqIO.
#     :param fasta_file: Full path to uncompressed, fasta-formatted file
#     :return: List of contig names.
#     """
#     contig_names = list()
#     for contig in SeqIO.parse(fasta_file, 'fasta'):
#         contig_names.append(contig.id)
#     return contig_names


def find_multibase_positions(column, ref_base, filtered_read_dict, base_cutoff, base_fraction_cutoff):
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
    total_coverage = 0
    base_count = dict()
    for category, base_dict in filtered_read_dict.items():
        # print(category)
        # if len(base_dict) > 1:
        # print(category)
        for base, count in base_dict.items():
            if 'filtered' not in category:
                if base not in base_count:
                    base_count[base] = count
                else:
                    base_count[base] += count
                total_coverage += count
                snv_dict['total'] += count
            # Forward and reverse reads agree on base
            if 'congruent' in category:
                snv_dict['total_congruent'] += count
                snv_dict['total_forward'] += int(count / 2)
                snv_dict['total_reverse'] += int(count / 2)
                if 'SNV' in category and base != ref_base:
                    snv_dict['total_congruent_SNV'] += count
                    snv_dict['total_forward_SNV'] += int(count / 2)
                    snv_dict['total_reverse_SNV'] += int(count / 2)
                    snv_dict['total_SNV'] += count
            elif category.startswith('forward_SNV'):
                snv_dict['total_forward'] += count
                if base != ref_base:
                    snv_dict['total_forward_SNV'] += count
                    snv_dict['total_SNV'] += count
            elif category.startswith('forward_ref'):
                snv_dict['total_forward'] += count
            elif category.startswith('reverse_SNV'):
                snv_dict['total_reverse'] += count
                if base != ref_base:
                    snv_dict['total_reverse_SNV'] += count
                    snv_dict['total_SNV'] += count
            elif category.startswith('reverse_ref'):
                snv_dict['total_reverse'] += count
            # elif 'SNV2' in category:
            #     snv_dict['total_forward'] += int(count / 2)
            #     snv_dict['total_forward_SNV'] += int(count / 2)
            #     snv_dict['total_reverse'] += int(count / 2)
            #     snv_dict['total_reverse_SNV'] += int(count / 2)
            #     snv_dict['total_SNV'] += count
            # print(column.reference_name, column.pos, base_dict)
            # print(base, count)
    # print(snv_dict)
    # multibase_list = list()
    # if snv_dict['total_congruent_SNV'] > 0:
    #     multibase_list.append('congruent')
    # try:
    #     if base_fraction_cutoff:
    #         if float(snv_dict['total_forward_SNV'] / snv_dict['total_forward']) >= base_fraction_cutoff and snv_dict['total_forward_SNV'] >= base_cutoff:
    #             multibase_list.append('forward')
    #     else:
    #         if snv_dict['total_forward_SNV'] >= base_cutoff:
    #             multibase_list.append('forward')
    # except ZeroDivisionError:
    #     pass
    # try:
    #     if base_fraction_cutoff:
    #         if float(snv_dict['total_reverse_SNV'] / snv_dict['total_reverse']) >= base_fraction_cutoff and snv_dict['total_reverse_SNV'] >= base_cutoff:
    #             multibase_list.append('reverse')
    #     else:
    #         if snv_dict['total_reverse_SNV'] >= base_cutoff:
    #             multibase_list.append('reverse')
    # except ZeroDivisionError:
    #     pass
    # try:
    #     if base_fraction_cutoff:
    #         if float(snv_dict['total_SNV'] / snv_dict['total']) >= base_fraction_cutoff and snv_dict['total_SNV'] >= base_cutoff:
    #             multibase_list.append('paired')
    #     else:
    #         if snv_dict['total_SNV'] >= base_cutoff:
    #             multibase_list.append('paired')
    # except ZeroDivisionError:
    #     pass
    passing_snv_dict = {
        'congruent': dict(),
        'forward': dict(),
        'reverse': dict(),
        'paired': dict()
    }

    '''
    empty_filtered_read_dict = {
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
    '''
    # if multibase_list:
    # print(column.reference_name, column.pos, multibase_list)
    return_dict = False
    for category, base_dict in filtered_read_dict.items():
        for base, count in base_dict.items():
            if 'congruent' in category and base != ref_base:
                # passing_snv_dict['congruent'][base] = count
                if base_fraction_cutoff:
                    if float(base_count[base] / snv_dict['total']) >= base_fraction_cutoff and base_count[
                        base] >= base_cutoff:
                        if base not in passing_snv_dict['forward']:
                            passing_snv_dict['congruent'][base] = count
                        else:
                            passing_snv_dict['congruent'][base] += count
                        return_dict = True
                else:
                    if base_count[base] >= base_cutoff:
                        if base not in passing_snv_dict['congruent']:
                            passing_snv_dict['congruent'][base] = count
                        else:
                            passing_snv_dict['congruent'][base] += count
                    return_dict = True
            if category.startswith('forward_SNV') and base != ref_base:
                # if column.reference_name == 'b0014_11' and column.pos == 425:
                #     print(category, base, count)
                if base_fraction_cutoff:
                    if float(base_count[base] / snv_dict['total']) >= base_fraction_cutoff and base_count[
                        base] >= base_cutoff:
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
            if category.startswith('reverse_SNV') and base != ref_base:
                if base_fraction_cutoff:
                    if float(base_count[base] / snv_dict['total']) >= base_fraction_cutoff and base_count[
                        base] >= base_cutoff:
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
            # and category.startswith('forward_SNV') or category.startswith('reverse_SNV')
            if base != ref_base and 'filtered' not in category:
                if base_fraction_cutoff:
                    if float(base_count[base] / snv_dict['total']) >= base_fraction_cutoff and base_count[
                        base] >= base_cutoff:
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
    # print(passing_snv_dict)
    # print(base_cutoff, base_fraction_cutoff, snv_dict['total_forward_SNV'])
    if return_dict:
        return snv_dict, passing_snv_dict, total_coverage
    else:
        return snv_dict, dict(), total_coverage
    # else:
    #     return snv_dict, dict(), total_coverage
    # for category, base_dict in filtered_read_dict.items():
    #     # if len(base_dict) > 1:
    #     # print(category)
    #     for base, count in base_dict.items():
    #         if 'congruent' in category:
    #             passing_snv_dict['congruent'].append(base)
    #         elif category.startswith('forward_SNV_'):
    # return snv_dict, passing_snv_dict


def position_details(actual_position, cutoff_dict, contig_name, ref_base, total_coverage):
    read_types = ['congruent', 'paired', 'forward', 'reverse']
    to_write = '{contig_name},{pos},{ref_base},'.format(contig_name=contig_name,
                                                        pos=actual_position,
                                                        ref_base=ref_base)
    snv_coverage = 0
    for read_type in read_types:
        # print(read_type, actual_position, cutoff_dict[read_type])
        semi_colon = False
        for base, coverage in cutoff_dict[read_type].items():
            if base:
                if semi_colon:
                    to_write += ';'
                to_write += '{base}:{depth}'.format(base=base,
                                                    depth=coverage)
                semi_colon = True
                if read_type != 'paired':
                    snv_coverage += coverage
        to_write += ','
    to_write += '{snv_cov},{total_cov}\n'.format(snv_cov=snv_coverage,
                                                 total_cov=total_coverage)
    return to_write


def read_contig(contig_name, bamfile_name, reference_fasta, allele_records, fastq_records, quality_cutoff=20,
                base_cutoff=2, base_fraction_cutoff=None, fasta=False):
    """
    Examines a contig to find if there are positions where more than one base is present.
    :param contig_name: Name of contig as a string.
    :param bamfile_name: Full path to bamfile. Must be sorted/indexed
    :param reference_fasta: Full path to fasta file that was used to generate the bamfile.
    :param quality_cutoff: Bases must have at least this phred score to be considered (INT)
    :param base_cutoff: At least this many bases must support a minor variant (INT)
    :param base_fraction_cutoff: At least this percentage of bases must support minor variant (FLOAT)
    :param fasta: Boolean on whether the samples are in FASTA format. Default is False
    :return: Dictionary of positions where more than one base is present. Keys are contig name, values are positions
    """
    bamfile = pysam.AlignmentFile(bamfile_name, 'rb')
    multibase_position_dict = dict()
    to_write = str()
    # If analysing FASTA files, a single base difference is all that is expected
    if fasta:
        base_cutoff = 1
    reference_sequence = str(allele_records[contig_name].seq)
    pysam_fasta = pysam.FastaFile(reference_fasta)
    # These parameters seem to be fairly undocumented with pysam, but I think that they should make the output
    # that I'm getting to match up with what I'm seeing in Tablet.
    pileup = bamfile.pileup(contig_name,
                            stepper='samtools',
                            ignore_orphans=False,
                            fastafile=pysam_fasta,
                            min_base_quality=0)
    # pileup = bamfile.pileup(contig_name,
    #                              stepper='samtools',
    #                              ignore_orphans=False,
    #                              fastafile=pysam_fasta,
    #                              min_base_quality=0)

    # for read in bamfile.pileup(contig_name,
    #                            stepper='samtools',
    #                            fastafile=pysam_fasta, ignore_orphans=False, min_base_quality=0):
    #     print(dir(read))
    #     quit()
    # # quit()
    for column in pileup:
        # print(column.reference_name, column.pos)
        # quit()

        filtered_read_dict = characterise_read(column=column,
                                               reference_sequence=reference_sequence,
                                               fastq_records=fastq_records,
                                               quality_cutoff=quality_cutoff,
                                               fasta=fasta)
        ref_base = reference_sequence[column.pos]
        snv_dict, cutoff_dict, total_coverage = find_multibase_positions(column=column,
                                                                         ref_base=ref_base,
                                                                         filtered_read_dict=filtered_read_dict,
                                                                         base_cutoff=base_cutoff,
                                                                         base_fraction_cutoff=base_fraction_cutoff)
        if cutoff_dict:
            # Pysam starts counting at 0, whereas we actually want to start counting at 1.
            actual_position = column.pos + 1
            if column.reference_name not in multibase_position_dict:
                multibase_position_dict[column.reference_name] = dict()

            multibase_position_dict[column.reference_name].update({actual_position: cutoff_dict})
            to_write += position_details(actual_position=actual_position,
                                         cutoff_dict=cutoff_dict,
                                         contig_name=contig_name,
                                         ref_base=ref_base,
                                         total_coverage=total_coverage)
            # to_write.append('{reference},{position},{bases},{coverage}\n'
            #                 .format(reference=column.reference_name,
            #                         position=actual_position,
            #                         bases=base_dict_to_string(base_dict),
            #                         coverage=sum(base_dict.values())))
    bamfile.close()
    return multibase_position_dict, to_write


def find_rmlst_type(kma_report, rmlst_report):
    """
    Uses a report generated by KMA to determine what allele is present for each rMLST gene.
    :param kma_report: The .res report generated by KMA.
    :param rmlst_report: rMLST report file to write information to.
    :return: a sorted list of loci present, in format gene_allele
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


def estimate_percent_contamination(contamination_report_file):
    """
    Estimates the percent contamination of a sample (and standard deviation).
    :param contamination_report_file: File created by read_contig,
    :return: Estimated percent contamination and standard deviation.
    """
    contam_levels = list()
    with open(contamination_report_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # lowest_count = 99999
            # base_counts = row['Bases'].split(';')
            base_counts = int(row['SNVCoverage'])
            total_coverage = int(row['TotalCoverage'])
            # for count in base_counts:
            #     num_bases = int(count.split(':')[1])
            #     if num_bases < lowest_count:
            #         lowest_count = num_bases
            # total_coverage = int(row['Coverage'])
            contam_levels.append(base_counts * 100 / total_coverage)
    return '%.2f' % (np.mean(contam_levels)), '%.2f' % np.std(contam_levels)


def load_fastq_records(gz, paired, forward):
    records = dict()
    for record in SeqIO.parse(gz, 'fastq'):
        if paired:
            if forward:
                if ':1:' in record.id:
                    # print(record.id, type(record.id))
                    record.id = record.id + '/1'
                    # print(record.id, type(record.id))
                    # records.update(SeqIO.to_dict([record]))
                elif '/1' in record.id:
                    pass
                    # records.update(SeqIO.to_dict(record))
                else:
                    record.id = record.id + '/1'
                    # records.update(SeqIO.to_dict(record))
            else:
                if ':2:' in record.id:
                    record.id = record.id + '/2'
                    # records.update(SeqIO.to_dict(record))
                elif '/2' in record.id:
                    pass
                    # records.update(SeqIO.to_dict(record))
                else:
                    record.id = record.id + '/2'
                    # records.update(SeqIO.to_dict(record))
        else:
            pass
            # records.update(SeqIO.to_dict(record))
        records.update(SeqIO.to_dict([record]))
    return records


def find_contamination(pair, output_folder, databases_folder, sample_name, forward_id='_R1', threads=1,
                       keep_files=False,
                       quality_cutoff=20, base_cutoff=2, base_fraction_cutoff=0.05, cgmlst_db=None, xmx=None,
                       tmpdir=None, data_type='Illumina', use_rmlst=False, cross_details=False, min_matching_hashes=40,
                       fasta=False):
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
    :param cross_details: If False, stop workflow when cross contamination is detected. If True, continue so estimates
    of percent contamination can be found (BOOL)
    :param min_matching_hashes: Minimum number of matching hashes in a MASH screen in order for a genus to be
    considered present in a sample. Default is 40
    :param fasta: Boolean on whether the samples are in FASTA format. Default is False
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
    if len(genus.split(':')) > 1:
        if not cross_details:
            write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
                         sample_name=sample_name,
                         multi_positions=0,
                         genus=genus,
                         percent_contam='ND',
                         contam_stddev='ND',
                         total_gene_length=0,
                         database_download_date=database_download_date)
            logging.info('Found cross-contamination! Skipping rest of analysis...\n')
            if keep_files is False:
                shutil.rmtree(sample_tmp_dir)
            return
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
                        allele_list = find_genusspecific_allele_list(os.path.join(db_folder, 'gene_allele.txt'),
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
                        allele_list = find_genusspecific_allele_list(os.path.join(db_folder, 'gene_allele.txt'),
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
                     percent_contam='ND',
                     contam_stddev='ND',
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
                # fastq_records = SeqIO.to_dict(SeqIO.parse(gz, 'fastq'))
                # for record in SeqIO.parse(gz, 'fastq'):
                #     print(record.id)
                fastq_records = load_fastq_records(gz=gz,
                                                   paired=True,
                                                   forward=True)
            with gzip.open(reverse_trimmed, 'rt') as gz:
                # fastq_records.update(SeqIO.to_dict(SeqIO.parse(gz, 'fastq')))
                fastq_records.update(load_fastq_records(gz=gz,
                                                        paired=True,
                                                        forward=False))
            # for record in fastq_records:
            #     print(record)
            #     break
            # for record in SeqIO.parse(gz, 'fasta'):
            #     fastq

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
    if not os.path.isfile(sample_database + '.fai'):  # Don't bother re-indexing, this only needs to happen once.
        pysam.faidx(sample_database)
    kma_database = sample_database.replace('.fasta', '') + '_kma'
    kma_report = os.path.join(sample_tmp_dir, '{sn}_kma'.format(sn=sample_name))
    if not os.path.isfile(kma_database + '.name'):  # The .name is one of the files KMA creates when making a database.
        logging.info('Since this is the first time you are using this database, it needs to be indexed by KMA. '
                     'This might take a while')
        cmd = 'kma index -i {} -o {}'.format(sample_database, kma_database)  # NOTE: Need KMA >=1.2.0 for this to work.
        out, err = run_cmd(cmd)
        write_to_logfile(log, out, err, cmd)

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
                # Use the FASTA file (rather than the readsd) as the input
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
                      '-t {threads}'.format(input_reads=unpaired_trimmed,
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
            # cmd = 'samtools faidx {ref}'.format(ref=rmlst_fasta)
            # if not os.path.isfile(rmlst_fasta + '.fai'):
            #     out, err = run_cmd(cmd)
            #     write_to_logfile(log, out, err, cmd)
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
                    # Lots of core genes seem to have relatives within a genome that are at ~70 percent identity. This means
                    # that reads that shouldn't map do, and cause false positives. Adding in this sub-filter means that
                    # reads can only have one mismatch, so they actually have to be from the right gene for this to work.
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
                    else:
                        ax = 'map-ont'
                    # '> {outsam}'
                    cmd = 'minimap2 --MD -t {threads} -ax {ax} {ref} {reads}' \
                        .format(ax=ax,
                                ref=rmlst_fasta,
                                reads=unpaired_trimmed,
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
                    # Apparently have to perform equivalent of a touch on this file for this to work.
                    # fh = open(outbam, 'w')
                    # fh.close()
                    # pysam.view('-b', '-o', outbam, outsam, save_stdout=outbam)

        if not os.path.isfile(sorted_bam):
            pysam.sort('-o', sorted_bam, outbam)
        if not os.path.isfile(sorted_bam + '.bai'):
            pysam.index(sorted_bam)
        # Now find number of multi-positions for each rMLST gene/allele combination
        multi_positions = 0

        # Run the BAM parsing in parallel! Some refactoring of the code would likely be a good idea so this
        # isn't quite so ugly, but it works.
        p = multiprocessing.Pool(processes=threads)
        allele_records = SeqIO.to_dict(SeqIO.parse(rmlst_fasta, 'fasta'))
        bamfile_list = [sorted_bam] * len(gene_alleles)
        # bamfile_list = [os.path.join(sample_tmp_dir, 'rmlst.bam')] * len(gene_alleles)
        reference_fasta_list = [rmlst_fasta] * len(gene_alleles)
        fasta_list = [fasta] * len(gene_alleles)
        quality_cutoff_list = [quality_cutoff] * len(gene_alleles)
        base_cutoff_list = [base_cutoff] * len(gene_alleles)
        base_fraction_list = [base_fraction_cutoff] * len(gene_alleles)
        records_list = [allele_records] * len(gene_alleles)
        fastq_records_list = [fastq_records] * len(gene_alleles)
        multibase_dict_list = list()
        report_write_list = list()
        # gene_alleles = ['b0169_2']
        '''
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
                fasta=fasta_list[i])
            multibase_dict_list.append(multibase_dict)
            report_write_list.append(report_write)
        '''
        for multibase_dict, report_write in p.starmap(read_contig,
                                                      zip(gene_alleles,
                                                          bamfile_list,
                                                          reference_fasta_list,
                                                          records_list,
                                                          fastq_records_list,
                                                          quality_cutoff_list,
                                                          base_cutoff_list,
                                                          base_fraction_list,
                                                          fasta_list),
                                                      chunksize=1):
            multibase_dict_list.append(multibase_dict)
            report_write_list.append(report_write)
        p.close()
        p.join()
        # '''
    except SamtoolsError:
        pysam_pass = False
        multi_positions = 0
        multibase_dict_list = list()
        report_write_list = list()
    # Write out report info.
    report_file = os.path.join(output_folder, sample_name + '_contamination.csv')
    with open(report_file, 'w') as r:
        r.write('Gene,Position,RefBase,CongruentSNVs,TotalSNVs,ForwardSNVs,ReverseSNVs,SNVCoverage,TotalCoverage\n')
        # r.write('{reference},{position},{congruent},{},{coverage}\n'.format(reference='Gene',
        #                                                              position='Position',
        #                                                              bases='CongruentReads',
        #                                                              coverage='Coverage'))
        for item in report_write_list:
            for contamination_info in item:
                r.write(contamination_info)

    # Total up the number of multibase positions.
    for multibase_position_dict in multibase_dict_list:
        multi_positions += sum([len(snp_positions) for snp_positions in multibase_position_dict.values()])
    if cgmlst_db is None:
        snp_cutoff = int(rmlst_gene_length / 10000) + 1
    elif fasta:
        snp_cutoff = 1
    else:
        snp_cutoff = 10
    if multi_positions >= snp_cutoff:
        percent_contam, contam_stddev = estimate_percent_contamination(contamination_report_file=report_file)
    else:
        percent_contam = 0
        contam_stddev = 0
    logging.info('Done! Number of contaminating SNVs found: {}\n'.format(multi_positions))
    write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
                 sample_name=sample_name,
                 multi_positions=multi_positions,
                 genus=genus,
                 percent_contam=percent_contam,
                 contam_stddev=contam_stddev,
                 total_gene_length=rmlst_gene_length,
                 snp_cutoff=snp_cutoff,
                 database_download_date=database_download_date,
                 pysam_pass=pysam_pass)
    if keep_files is False:
        shutil.rmtree(sample_tmp_dir)


def write_output(output_report, sample_name, multi_positions, genus, percent_contam, contam_stddev, total_gene_length,
                 database_download_date, snp_cutoff=3, pysam_pass=True):
    """
    Function that writes the output generated by ConFindr to a report file. Appends to a file that already exists,
    or creates the file if it doesn't already exist.
    :param output_report: Path to CSV output report file. Should have headers SampleName,Genus,NumContamSNVs,
    ContamStatus,PercentContam, and PercentContamStandardDeviation, in that order.
    :param sample_name: string - name of sample
    :param multi_positions: integer - number of positions that were found to have more than one base present.
    :param genus: string - The genus of your sample
    :param percent_contam: float - Estimated percentage contamination
    :param contam_stddev: float - Standard deviation of percentage contamination
    :param total_gene_length: integer - number of bases examined to make a contamination call.
    :param database_download_date:
    :param snp_cutoff: Number of cSNVs to use to call a sample contaminated. Default 3. (INT)
    :param pysam_pass: Boolean of whether pysam encountered an error
    """
    # If the report file hasn't been created, make it, with appropriate header.
    if not os.path.isfile(output_report):
        with open(os.path.join(output_report), 'w') as f:
            f.write('Sample,Genus,NumContamSNVs,ContamStatus,PercentContam,PercentContamStandardDeviation,'
                    'BasesExamined,DatabaseDownloadDate\n')
    if pysam_pass:
        if multi_positions >= snp_cutoff or len(genus.split(':')) > 1:
            contaminated = True
        else:
            contaminated = False
    else:
        contaminated = 'Pysam SamtoolsError'
        multi_positions = 'ND'
        percent_contam = 'ND'
        contam_stddev = 'ND'
    with open(output_report, 'a+') as f:
        f.write('{samplename},{genus},{numcontamsnvs},'
                '{contamstatus},{percent_contam},{contam_stddev},'
                '{gene_length},{database_download_date}\n'.format(samplename=sample_name,
                                                                  genus=genus,
                                                                  numcontamsnvs=multi_positions,
                                                                  contamstatus=contaminated,
                                                                  percent_contam=percent_contam,
                                                                  contam_stddev=contam_stddev,
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
                               sample_name=sample_name,
                               keep_files=args.keep_files,
                               quality_cutoff=args.quality_cutoff,
                               base_cutoff=args.base_cutoff,
                               base_fraction_cutoff=args.base_fraction_cutoff,
                               cgmlst_db=args.cgmlst,
                               xmx=args.Xmx,
                               tmpdir=args.tmp,
                               data_type=args.data_type,
                               use_rmlst=args.rmlst,
                               cross_details=args.cross_details,
                               min_matching_hashes=min_matching_hashes,
                               fasta=args.fasta)
        except subprocess.CalledProcessError:
            # If something unforeseen goes wrong, traceback will be printed to screen.
            # We then add the sample to the report with a note that it failed.
            multi_positions = 0
            genus = 'Error processing sample'
            write_output(output_report=os.path.join(args.output_name, 'confindr_report.csv'),
                         sample_name=sample_name,
                         multi_positions=multi_positions,
                         genus=genus,
                         percent_contam='ND',
                         contam_stddev='ND',
                         total_gene_length=0,
                         database_download_date='ND')
            logging.warning('Encountered error when attempting to run ConFindr on sample '
                            '{sample}. Skipping...'.format(sample=sample_name))
            logging.warning('Error encounted was:\n{}'.format(traceback.format_exc()))
            if args.keep_files is False:
                shutil.rmtree(os.path.join(args.output_name, sample_name))
    if args.keep_files is False and args.tmp is not None:
        shutil.rmtree(args.tmp)
    logging.info('Contamination detection complete!')


def get_version():
    try:
        version = 'ConFindr {}'.format(pkg_resources.get_distribution('confindr').version)
    except pkg_resources.DistributionNotFound:
        version = 'ConFindr (Unknown version)'
    return version


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
                        default=2,
                        help='Number of bases necessary to support a multiple allele call. Defaults to 2.')
    parser.add_argument('-bf', '--base_fraction_cutoff',
                        type=float,
                        default=0.05,
                        help='Fraction of bases necessary to support a multiple allele call. Particularly useful when '
                             'dealing with very high coverage samples. Default is 0.05.')
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
    parser.add_argument('-cross_details', '--cross_details',
                        action='store_true',
                        help='Continue ConFindr analyses on samples with two or more genera identified. Default is '
                             'False')
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
