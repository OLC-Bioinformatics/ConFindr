#!/usr/bin/env python

import multiprocessing
import urllib.request
import numpy as np
import subprocess
import argparse
import tarfile
import logging
import shutil
import glob
import csv
import os
import pysam
from Bio import SeqIO
from confindr_wrappers import mash
from confindr_wrappers import bbtools
from confindr_wrappers import nanopore_methods


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


def find_unpaired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
    """
    Looks at a directory to find unpaired fastq files.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default _R1.
    :param reverse_id: Identifier for forward reads. Default _R2.
    :return: List of files that appear to be unpaired reads.
    """
    read_list = list()
    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
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
    :param target_genus:
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


def setup_allelespecific_database(database_folder, genus, allele_list):
    """
    Since some genera have some rMLST genes missing, or two copies of some genes, genus-specific databases are needed.
    This will take only the alleles known to be part of each genus and write them to a genus-specific file.
    :param database_folder: Path to folder where confindr databases are stored.
    :param genus: Genus of organism, as a string. First letter should be capitalized, everything else lowercase
    :param allele_list: allele list generated by find_genusspecific_allele_list
    """
    with open(os.path.join(database_folder, '{}_db.fasta'.format(genus)), 'w') as f:
        sequences = SeqIO.parse(os.path.join(database_folder, 'rMLST_combined.fasta'), 'fasta')
        for item in sequences:
            if item.id in allele_list:
                f.write('>' + item.id + '\n')
                f.write(str(item.seq) + '\n')


def extract_rmlst_genes(pair, database, forward_out, reverse_out, threads=12, logfile=None):
    """
    Given a pair of reads and an rMLST database, will extract reads that contain sequence from the database.
    :param pair: List containing path to forward reads at index 0 and path to reverse reads at index 1.
    :param database: Path to rMLST database, in FASTA format.
    :param forward_out:
    :param reverse_out:
    :param threads:
    """
    out, err, cmd = bbtools.bbduk_bait(database, pair[0], forward_out, reverse_in=pair[1],
                                       reverse_out=reverse_out, threads=str(threads), returncmd=True)
    if logfile:
        write_to_logfile(logfile, out, err, cmd)


def find_cross_contamination(databases, pair, tmpdir='tmp', log='log.txt', threads=1):
    """
    Usese mash to find out whether or not a sample has more than one genus present, indicating cross-contamination.
    :param databases: A databases folder, which must contain refseq.msh, a mash sketch that has one representative
    per genus from refseq.
    :param tmpdir: Temporary directory to store mash result files in.
    :param pair: Array with path to forward reads at index 0 and path to reverse reads at index o
    :param log: Logfile to write to.
    :param threads: Number of threads to run mash wit.
    :return: cross_contam: a bool that is True if more than one genus is found, and False otherwise.
    :return: genera_present: A string. If only one genus is found, string is just genus. If more than one genus is found,
    the string is a list of genera present, separated by colons (i.e. for Escherichia and Salmonella found, string would
    be 'Escherichia:Salmonella'. If no genus found, return 'NA'
    """
    genera_present = list()
    out, err, cmd = mash.screen('{}/refseq.msh'.format(databases), pair[0],
                                pair[1], threads=threads, w='', i='0.95',
                                output_file=os.path.join(tmpdir, 'screen.tab'), returncmd=True)
    write_to_logfile(log, out, err, cmd)
    screen_output = mash.read_mash_screen(os.path.join(tmpdir, 'screen.tab'))
    for item in screen_output:
        mash_genus = item.query_id.split('/')[-3]
        if mash_genus == 'Shigella':
            mash_genus = 'Escherichia'
        if mash_genus not in genera_present:
            genera_present.append(mash_genus)
    if len(genera_present) == 1:
        genera_present = genera_present[0]
    elif len(genera_present) == 0:
        genera_present = 'NA'
    else:
        tmpstr = ''
        for mash_genus in genera_present:
            tmpstr += mash_genus + ':'
        genera_present = tmpstr[:-1]
    return genera_present


def find_cross_contamination_unpaired(databases, reads, tmpdir='tmp', log='log.txt', threads=1):
    """
    Usese mash to find out whether or not a sample has more than one genus present, indicating cross-contamination.
    :param databases: A databases folder, which must contain refseq.msh, a mash sketch that has one representative
    per genus from refseq.
    :param tmpdir: Temporary directory to store mash result files in.
    :param log: Logfile to write to.
    :param threads: Number of threads to run mash wit.
    :return: cross_contam: a bool that is True if more than one genus is found, and False otherwise.
    :return: genera_present: A string. If only one genus is found, string is NA. If more than one genus is found,
    the string is a list of genera present, separated by colons (i.e. for Escherichia and Salmonella found, string would
    be 'Escherichia:Salmonella'
    """
    genera_present = list()
    out, err, cmd = mash.screen('{}/refseq.msh'.format(databases), reads,
                                threads=threads, w='', i='0.95',
                                output_file=os.path.join(tmpdir, 'screen.tab'), returncmd=True)
    write_to_logfile(log, out, err, cmd)
    screen_output = mash.read_mash_screen(os.path.join(tmpdir, 'screen.tab'))
    for item in screen_output:
        mash_genus = item.query_id.split('/')[-3]
        if mash_genus == 'Shigella':
            mash_genus = 'Escherichia'
        if mash_genus not in genera_present:
            genera_present.append(mash_genus)
    if len(genera_present) == 1:
        genera_present = genera_present[0]
    elif len(genera_present) == 0:
        genera_present = 'NA'
    else:
        tmpstr = ''
        for mash_genus in genera_present:
            tmpstr += mash_genus + ':'
        genera_present = tmpstr[:-1]
    return genera_present


# TODO: Maybe rename this, seeing how the number of high quality bases can now be user-defined.
def has_two_high_quality_bases(list_of_scores, quality_cutoff=20, base_count_cutoff=2):
    """
    Finds if a site has at least two bases of  high quality, enough that it can be considered
    fairly safe to say that base is actually there.
    :param list_of_scores: List of quality scores as integers.
    :param quality_cutoff: Phred quality bases need to have in order to support multiple allele presence.
    :param base_count_cutoff: Number of bases needed to support multiple allele presence.
    :return: True if site has at least two bases >= 20 phred score, false otherwise (changeable by user)
    """
    quality_bases_count = 0
    for score in list_of_scores:
        if score >= quality_cutoff:
            quality_bases_count += 1
        if quality_bases_count >= base_count_cutoff:
            return True
    return False


def find_if_multibase(column, quality_cutoff, base_cutoff):
    """
    Finds if a position in a pileup has more than one base present.
    :param column: A pileupColumn generated by pysam
    :param quality_cutoff: Desired minimum phred quality for a base in order to be counted towards a multi-allelic column
    :param base_cutoff: Minimum number of bases needed to support presence of a base.
    :return: If position has more than one base, a dictionary with counts for the bases. Otherwise, returns
    empty dictionary
    """
    base_dict = dict()
    # Sometimes the qualities come out to ridiculously high (>70) values. Looks to be because sometimes reads
    # are overlapping and the qualities get summed for overlapping bases. Issue opened on pysam.
    base_qualities = dict()
    for read in column.pileups:
        if read.query_position is not None:  # Not entirely sure why this is sometimes None, but it causes bad stuff
            base = read.alignment.query_sequence[read.query_position]
            quality = read.alignment.query_qualities[read.query_position]
            if base not in base_qualities:
                base_qualities[base] = [quality]
            else:
                base_qualities[base].append(quality)
            if base in base_dict:
                base_dict[base] += 1
            else:
                base_dict[base] = 1
    # If we only have one base, ignore things.
    if len(base_dict) == 1:
        base_dict = dict()
    # Now check that at least two bases for each of the bases present high quality.
    if base_dict:
        high_quality_bases = True
        for base in base_qualities:
            if has_two_high_quality_bases(base_qualities[base], base_count_cutoff=base_cutoff, quality_cutoff=quality_cutoff) is False:
                high_quality_bases = False
        if high_quality_bases is False:
            base_dict = dict()
    return base_dict


def get_contig_names(fasta_file):
    """
    Gets contig names from a fasta file using SeqIO.
    :param fasta_file: Full path to uncompressed, fasta-formatted file
    :return: List of contig names.
    """
    contig_names = list()
    for contig in SeqIO.parse(fasta_file, 'fasta'):
        contig_names.append(contig.id)
    return contig_names


def read_contig(contig_name, bamfile_name, reference_fasta, report_file, quality_cutoff=20, base_cutoff=2):
    """
    Examines a contig to find if there are positions where more than one base is present.
    :param contig_name: Name of contig as a string.
    :param bamfile_name: Full path to bamfile. Must be sorted/indexed
    :param reference_fasta: Full path to fasta file that was used to generate the bamfile.
    :param report_file: File where information about each position found to have multiple bases
    will be written.
    :return: Dictionary of positions where more than one base is present. Keys are contig name, values are positions
    """
    bamfile = pysam.AlignmentFile(bamfile_name, 'rb')
    multibase_position_dict = dict()
    # These parameters seem to be fairly undocumented with pysam, but I think that they should make the output
    # that I'm getting to match up with what I'm seeing in Tablet.
    for column in bamfile.pileup(contig_name,
                                 stepper='samtools', ignore_orphans=False, fastafile=pysam.FastaFile(reference_fasta),
                                 min_base_quality=0):
        base_dict = find_if_multibase(column, quality_cutoff=quality_cutoff, base_cutoff=base_cutoff)

        if base_dict:
            if column.reference_name in multibase_position_dict:
                multibase_position_dict[column.reference_name].append(column.pos)
            else:
                multibase_position_dict[column.reference_name] = [column.pos]
            with open(report_file, 'a+') as r:
                r.write('{reference},{position},{bases},{coverage}\n'.format(reference=column.reference_name,
                                                                             position=column.pos,
                                                                             bases=base_dict_to_string(base_dict),
                                                                             coverage=sum(base_dict.values())))
    bamfile.close()
    return multibase_position_dict


def find_rmlst_type(bamfile, sample_database, rmlst_report):
    """
    Approximates finding the rMLST type by counting the number of reads aligned to each allele and reporting
    which allele has the most reads aligned.
    :param bamfile: Path to bamfile where reads where aligned to rmlst file. It's assumed these reads were aligned
    using bbmap, with ambig=all set as an option.
    :param sample_database: Path to fasta-formatted sample database (usually a genus-specific rMLST)
    :param rmlst_report: Path to file where rmlst report (what allele each gene is) is written.
    :return: gene_alleles: Dictionary where keys are genes, and values are alleles
    """
    # Find which rMLST allele has the most reads mapped back to it for each gene.
    gene_alleles_to_use = dict()
    bamfile = pysam.AlignmentFile(os.path.join(bamfile), 'rb')
    contig_names = get_contig_names(fasta_file=sample_database)
    for contig in contig_names:
        gene = contig.split('_')[0]
        allele = contig.split('_')[1]
        if gene not in gene_alleles_to_use:
            gene_alleles_to_use[gene] = [allele, bamfile.count(contig=contig)]
        else:
            read_count = bamfile.count(contig=contig)
            if gene_alleles_to_use[gene][1] < read_count:
                gene_alleles_to_use[gene] = [allele, read_count]
    # Extract only the genes with most reads aligned, then align reads to those.
    gene_alleles = list()
    for gene in gene_alleles_to_use:
        contig_name = gene + '_' + gene_alleles_to_use[gene][0]
        gene_alleles.append(contig_name)
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
    for base in base_dict:
        outstr += base + ':' + str(base_dict[base]) + ';'
    return outstr[:-1]


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
            lowest_count = 99999
            base_counts = row['Bases'].split(';')
            for count in base_counts:
                num_bases = int(count.split(':')[1])
                if num_bases < lowest_count:
                    lowest_count = num_bases
            total_coverage = int(row['Coverage'])
            contam_levels.append(lowest_count*100/total_coverage)
    return '%.2f' % (np.mean(contam_levels)), '%.2f' % np.std(contam_levels)


def find_contamination(pair, output_folder, databases_folder, forward_id='_R1', threads=1, keep_files=False,
                       quality_cutoff=20, base_cutoff=2):
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
    :return:
    """
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

    # Get the output report file set up
    if not os.path.isfile(os.path.join(output_folder, 'confindr_report.csv')):
        with open(os.path.join(output_folder, 'confindr_report.csv'), 'w') as f:
            f.write('Sample,Genus,NumContamSNVs,ContamStatus,PercentContam,PercentContamStandardDeviation\n')
    logging.info('Checking for cross-species contamination...')
    if paired:
        genus = find_cross_contamination(databases_folder, pair, tmpdir=sample_tmp_dir, log=log, threads=threads)
    else:
        genus = find_cross_contamination_unpaired(databases_folder, reads=pair[0], tmpdir=sample_tmp_dir, log=log, threads=threads)
    if len(genus.split(':')) > 1:
        write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
                     sample_name=sample_name,
                     multi_positions=0,
                     genus=genus,
                     percent_contam='NA',
                     contam_stddev='NA')
        logging.info('Found cross-contamination! Skipping rest of analysis...\n')
        shutil.rmtree(sample_tmp_dir)
        return
    # Main method for finding contamination - works on one pair at a time.
    # Need to:
    # Setup genus-specific databases, if necessary.
    if genus != 'NA':
        sample_database = os.path.join(databases_folder, '{}_db.fasta'.format(genus))
        if not os.path.isfile(os.path.join(databases_folder, '{}_db.fasta'.format(genus))):
            logging.info('Setting up genus-specific database for genus {}...'.format(genus))
            allele_list = find_genusspecific_allele_list(os.path.join(databases_folder, 'gene_allele.txt'), genus)
            setup_allelespecific_database(databases_folder, genus, allele_list)
    else:
        sample_database = os.path.join(databases_folder, 'rMLST_combined.fasta')
    # Extract rMLST reads and quality trim.
    logging.info('Extracting rMLST genes...')
    if paired:
        out, err, cmd = bbtools.bbduk_bait(reference=sample_database,
                                           forward_in=pair[0],
                                           reverse_in=pair[1],
                                           forward_out=os.path.join(sample_tmp_dir, 'rmlst_R1.fastq.gz'),
                                           reverse_out=os.path.join(sample_tmp_dir, 'rmlst_R2.fastq.gz'),
                                           threads=threads,
                                           returncmd=True)
    else:
        out, err, cmd = bbtools.bbduk_bait(reference=sample_database, forward_in=pair[0],
                                           forward_out=os.path.join(sample_tmp_dir, 'rmlst.fastq.gz'),
                                           returncmd=True, threads=threads)
    write_to_logfile(log, out, err, cmd)
    logging.info('Quality trimming...')
    if paired:
        out, err, cmd = bbtools.bbduk_trim(forward_in=os.path.join(sample_tmp_dir, 'rmlst_R1.fastq.gz'),
                                           reverse_in=os.path.join(sample_tmp_dir, 'rmlst_R2.fastq.gz'),
                                           forward_out=os.path.join(sample_tmp_dir, 'trimmed_R1.fastq.gz'),
                                           reverse_out=os.path.join(sample_tmp_dir, 'trimmed_R2.fastq.gz'),
                                           threads=str(threads), returncmd=True)
    else:
        out, err, cmd = bbtools.bbduk_trim(forward_in=os.path.join(sample_tmp_dir, 'rmlst.fastq.gz'),
                                           forward_out=os.path.join(sample_tmp_dir, 'trimmed.fastq.gz'),
                                           returncmd=True, threads=threads)
    write_to_logfile(log, out, err, cmd)

    logging.info('Detecting contamination...')
    # Now do mapping in two steps - first, map reads back to database with ambiguous reads matching all - this
    # will be used to get a count of number of reads aligned to each gene/allele so we can create a custom rmlst file
    # with only the most likely allele for each gene.
    cmd = 'samtools faidx {}'.format(sample_database)
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    if paired:
        cmd = 'bbmap.sh ref={ref} in={forward_in} in2={reverse_in} out={outbam} ' \
              'nodisk ambig=all threads={threads}'.format(ref=sample_database,
                                                          forward_in=os.path.join(sample_tmp_dir, 'trimmed_R1.fastq.gz'),
                                                          reverse_in=os.path.join(sample_tmp_dir, 'trimmed_R2.fastq.gz'),
                                                          outbam=os.path.join(sample_tmp_dir, 'out.bam'),
                                                          threads=threads)
    else:
        cmd = 'bbmap.sh ref={ref} in={forward_in} out={outbam} threads={threads} ' \
              'nodisk ambig=all'.format(ref=sample_database,
                                        forward_in=os.path.join(sample_tmp_dir, 'trimmed.fastq.gz'),
                                        outbam=os.path.join(sample_tmp_dir, 'out.bam'),
                                        threads=threads)
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    cmd = 'samtools sort {inbam} -o {sorted_bam}'.format(inbam=os.path.join(sample_tmp_dir, 'out.bam'),
                                                         sorted_bam=os.path.join(sample_tmp_dir, 'rmlst.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    cmd = 'samtools index {sorted_bam}'.format(sorted_bam=os.path.join(sample_tmp_dir, 'rmlst.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)

    rmlst_report = os.path.join(output_folder, sample_name + '_rmlst.csv')
    gene_alleles = find_rmlst_type(bamfile=os.path.join(sample_tmp_dir, 'rmlst.bam'),
                                   sample_database=sample_database,
                                   rmlst_report=rmlst_report)

    with open(os.path.join(sample_tmp_dir, 'rmlst.fasta'), 'w') as f:
        for contig in SeqIO.parse(sample_database, 'fasta'):
            if contig.id in gene_alleles:
                f.write('>{}\n'.format(contig.id))
                f.write(str(contig.seq) + '\n')

    # Second step of mapping - Do a mapping of our baited reads against a fasta file that has only one allele per
    # rMLST gene.
    cmd = 'samtools faidx {}'.format(os.path.join(sample_tmp_dir, 'rmlst.fasta'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    if paired:
        cmd = 'bbmap.sh ref={ref} in={forward_in} in2={reverse_in} out={outbam} threads={threads} ' \
              'nodisk'.format(ref=os.path.join(sample_tmp_dir, 'rmlst.fasta'),
                              forward_in=os.path.join(sample_tmp_dir, 'trimmed_R1.fastq.gz'),
                              reverse_in=os.path.join(sample_tmp_dir, 'trimmed_R2.fastq.gz'),
                              outbam=os.path.join(sample_tmp_dir, 'out_2.bam'),
                              threads=threads)
    else:
        cmd = 'bbmap.sh ref={ref} in={forward_in} out={outbam} threads={threads} ' \
              'nodisk'.format(ref=os.path.join(sample_tmp_dir, 'rmlst.fasta'),
                              forward_in=os.path.join(sample_tmp_dir, 'trimmed.fastq.gz'),
                              outbam=os.path.join(sample_tmp_dir, 'out_2.bam'),
                              threads=threads)
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    cmd = 'samtools sort {inbam} -o {sorted_bam}'.format(inbam=os.path.join(sample_tmp_dir, 'out_2.bam'),
                                                         sorted_bam=os.path.join(sample_tmp_dir, 'contamination.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    cmd = 'samtools index {sorted_bam}'.format(sorted_bam=os.path.join(sample_tmp_dir, 'contamination.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    # Now find number of multi-positions for each rMLST gene/allele combination
    report_file = os.path.join(output_folder, sample_name + '_contamination.csv')
    with open(report_file, 'w') as r:
        r.write('{reference},{position},{bases},{coverage}\n'.format(reference='Gene',
                                                                     position='Position',
                                                                     bases='Bases',
                                                                     coverage='Coverage'))
    multi_positions = 0
    for contig_name in gene_alleles:
        multibase_position_dict = read_contig(contig_name=contig_name,
                                              bamfile_name=os.path.join(sample_tmp_dir, 'contamination.bam'),
                                              reference_fasta=os.path.join(sample_tmp_dir, 'rmlst.fasta'),
                                              report_file=report_file,
                                              quality_cutoff=quality_cutoff,
                                              base_cutoff=base_cutoff)
        multi_positions += len(multibase_position_dict)
    if multi_positions >= 3:
        percent_contam, contam_stddev = estimate_percent_contamination(contamination_report_file=report_file)
    else:
        percent_contam = 'NA'
        contam_stddev = 'NA'
    logging.info('Done! Number of contaminating SNVs found: {}\n'.format(multi_positions))
    write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
                 sample_name=sample_name,
                 multi_positions=multi_positions,
                 genus=genus,
                 percent_contam=percent_contam,
                 contam_stddev=contam_stddev)
    if keep_files is False:
        shutil.rmtree(sample_tmp_dir)


def write_output(output_report, sample_name, multi_positions, genus, percent_contam, contam_stddev):
    if multi_positions > 2 or len(genus.split(':')) > 1:
        contaminated = True
    else:
        contaminated = False
    with open(output_report, 'a+') as f:
        f.write('{samplename},{genus},{numcontamsnvs},'
                '{contamstatus},{percent_contam},{contam_stddev}\n'.format(samplename=sample_name,
                                                                           genus=genus,
                                                                           numcontamsnvs=multi_positions,
                                                                           contamstatus=contaminated,
                                                                           percent_contam=percent_contam,
                                                                           contam_stddev=contam_stddev))


# TODO: This can be deleted soon, as it's no longer used. Keep it around until nanopore stuff gets tested more.
def find_contamination_unpaired(reads, output_folder, databases_folder, threads=1, keep_files=False, read_type='Illumina'):
    # Setup log file.
    log = os.path.join(output_folder, 'confindr_log.txt')
    # Setup a sample name - may want to improve this at some point, currently takes everything before the .fastq.gz
    sample_name = os.path.split(reads)[-1].split('.')[0]
    sample_tmp_dir = os.path.join(output_folder, sample_name)
    if not os.path.isdir(sample_tmp_dir):
        os.makedirs(sample_tmp_dir)
    logging.info('Checking for cross-species contamination...')
    genus = find_cross_contamination_unpaired(databases_folder, reads, tmpdir=sample_tmp_dir, log=log, threads=threads)
    if len(genus.split(':')) > 1:
        multi_positions = 0
        write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
                     sample_name=sample_name,
                     multi_positions=multi_positions,
                     genus=genus)
        logging.info('Found cross-contamination! Skipping rest of analysis...\n')
        shutil.rmtree(sample_tmp_dir)
        return
    # Setup a genusspecfic database, if necessary.
    if genus != 'NA':
        sample_database = os.path.join(databases_folder, '{}_db.fasta'.format(genus))
        if not os.path.isfile(os.path.join(databases_folder, '{}_db.fasta'.format(genus))):
            logging.info('Setting up genus-specific database for genus {}...'.format(genus))
            allele_list = find_genusspecific_allele_list(os.path.join(databases_folder, 'gene_allele.txt'), genus)
            setup_allelespecific_database(databases_folder, genus, allele_list)
    else:
        sample_database = os.path.join(databases_folder, 'rMLST_combined.fasta')
    # Get tmpdir for this sample created.
    sample_tmp_dir = os.path.join(output_folder, sample_name)
    if not os.path.isdir(sample_tmp_dir):
        os.makedirs(sample_tmp_dir)
    # With everything set up, time to start the workflow.
    # First thing to do: Extract rMLST genes.
    logging.info('Extracting rMLST genes...')
    out, err, cmd = bbtools.bbduk_bait(reference=sample_database, forward_in=reads,
                                       forward_out=os.path.join(sample_tmp_dir, 'rmlst.fastq.gz'),
                                       returncmd=True, threads=threads)
    logging.debug('rMLST extraction command used: {}'.format(cmd))
    write_to_logfile(log, out, err, cmd)
    if read_type == 'Illumina':
        logging.info('Quality trimming...')
        # With rMLST genes extracted, get our quality trimming done.
        out, err, cmd = bbtools.bbduk_trim(forward_in=os.path.join(sample_tmp_dir, 'rmlst.fastq.gz'),
                                           forward_out=os.path.join(sample_tmp_dir, 'trimmed.fastq.gz'),
                                           returncmd=True, threads=threads)
        logging.debug('Quality trim command used: {}'.format(cmd))
        write_to_logfile(log, out, err, cmd)
    logging.info('Detecting contamination...')
    # Now do mapping in two steps - first, map reads back to database with ambiguous reads matching all - this
    # will be used to get a count of number of reads aligned to each gene/allele so we can create a custom rmlst file
    # with only the most likely allele for each gene.
    cmd = 'samtools faidx {}'.format(sample_database)
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    if read_type == 'Nanopore':
        cmd = 'minimap2 -ax map-ont -t {threads} {reference} {reads} | samtools view -bS > ' \
              '{output_bam}'.format(threads=threads,
                                    reference=sample_database,
                                    reads=os.path.join(sample_tmp_dir, 'rmlst.fastq.gz'),
                                    output_bam=os.path.join(sample_tmp_dir, 'out.bam'))
    elif read_type == 'Illumina':
        cmd = 'bbmap.sh ref={ref} in={forward_in} out={outbam} threads={threads} ' \
              'nodisk ambig=all'.format(ref=sample_database,
                                        forward_in=os.path.join(sample_tmp_dir, 'trimmed.fastq.gz'),
                                        outbam=os.path.join(sample_tmp_dir, 'out.bam'),
                                        threads=threads)
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    cmd = 'samtools sort {inbam} -o {sorted_bam}'.format(inbam=os.path.join(sample_tmp_dir, 'out.bam'),
                                                         sorted_bam=os.path.join(sample_tmp_dir, 'rmlst.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    cmd = 'samtools index {sorted_bam}'.format(sorted_bam=os.path.join(sample_tmp_dir, 'rmlst.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)

    rmlst_report = os.path.join(output_folder, sample_name + '_rmlst.csv')
    gene_alleles = find_rmlst_type(bamfile=os.path.join(sample_tmp_dir, 'rmlst.bam'),
                                   sample_database=sample_database,
                                   rmlst_report=rmlst_report)

    with open(os.path.join(sample_tmp_dir, 'rmlst.fasta'), 'w') as f:
        for contig in SeqIO.parse(sample_database, 'fasta'):
            if contig.id in gene_alleles:
                f.write('>{}\n'.format(contig.id))
                f.write(str(contig.seq) + '\n')

    # Second step of mapping - Do a mapping of our baited reads against a fasta file that has only one allele per
    # rMLST gene.
    cmd = 'samtools faidx {}'.format(os.path.join(sample_tmp_dir, 'rmlst.fasta'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    if read_type == 'Nanopore':
        cmd = 'minimap2 -L -ax map-ont -t {threads} {reference} {reads} | samtools view -bS > ' \
              '{output_bam}'.format(threads=threads,
                                    reference=os.path.join(sample_tmp_dir, 'rmlst.fasta'),
                                    reads=os.path.join(sample_tmp_dir, 'rmlst.fastq.gz'),
                                    output_bam=os.path.join(sample_tmp_dir, 'out_2.bam'))
    elif read_type == 'Illumina':
        cmd = 'bbmap.sh ref={ref} in={forward_in} out={outbam} threads={threads} ' \
              'nodisk'.format(ref=os.path.join(sample_tmp_dir, 'rmlst.fasta'),
                              forward_in=os.path.join(sample_tmp_dir, 'trimmed.fastq.gz'),
                              outbam=os.path.join(sample_tmp_dir, 'out_2.bam'),
                              threads=threads)
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    cmd = 'samtools sort {inbam} -o {sorted_bam}'.format(inbam=os.path.join(sample_tmp_dir, 'out_2.bam'),
                                                         sorted_bam=os.path.join(sample_tmp_dir, 'contamination.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    cmd = 'samtools index {sorted_bam}'.format(sorted_bam=os.path.join(sample_tmp_dir, 'contamination.bam'))
    out, err = run_cmd(cmd)
    write_to_logfile(log, out, err, cmd)
    # Now find number of multi-positions for each rMLST gene/allele combination
    report_file = os.path.join(output_folder, sample_name + '_contamination.csv')
    with open(report_file, 'w') as r:
        r.write('{reference},{position},{bases},{coverage}\n'.format(reference='Gene',
                                                                     position='Position',
                                                                     bases='Bases',
                                                                     coverage='Coverage'))
    multi_positions = 0
    for contig_name in gene_alleles:
        # TODO: This should be parallelizable (and ConFindr on MinION data is slooow). Get this done at some
        # point
        logging.debug(contig_name)
        multibase_position_dict = read_contig(contig_name=contig_name,
                                              bamfile_name=os.path.join(sample_tmp_dir, 'contamination.bam'),
                                              reference_fasta=os.path.join(sample_tmp_dir, 'rmlst.fasta'),
                                              report_file=report_file)
        multi_positions += len(multibase_position_dict)
    logging.info('Done! Number of contaminating SNVs found: {}\n'.format(multi_positions))
    write_output(output_report=os.path.join(output_folder, 'confindr_report.csv'),
                 sample_name=sample_name,
                 multi_positions=multi_positions,
                 genus=genus)
    if keep_files is False:
        shutil.rmtree(sample_tmp_dir)


def check_for_databases_and_download(database_location, tmpdir):
    # Check for the files necessary - should have rMLST_combined.fasta, gene_allele.txt, profiles.txt, and refseq.msh
    necessary_files = ['rMLST_combined.fasta', 'gene_allele.txt', 'profiles.txt', 'refseq.msh']
    all_files_present = True
    for necessary_file in necessary_files:
        if not os.path.isfile(os.path.join(database_location, necessary_file)):
            all_files_present = False

    if not all_files_present:
        logging.info('Databases not present. Downloading to {}. This may take a few minutes...'.format(database_location))
        if not os.path.isdir(tmpdir):
            os.makedirs(tmpdir)
        # Download
        # TODO: Progress bar?
        urllib.request.urlretrieve('https://ndownloader.figshare.com/files/11864267', os.path.join(tmpdir, 'confindr_db.tar.gz'))
        # Extract.
        tar = tarfile.open(os.path.join(tmpdir, 'confindr_db.tar.gz'))
        tar.extractall(path=database_location)
        tar.close()
        # We now have the files extracted to database_location/databases - move them to database_location
        for item in glob.glob(os.path.join(database_location, 'databases', '*')):
            shutil.move(item, database_location)
        # Cleanup!
        shutil.rmtree(os.path.join(database_location, 'databases'))
        os.remove(os.path.join(tmpdir, 'confindr_db.tar.gz'))
        logging.info('Databases successfully downloaded...')


if __name__ == '__main__':
    version = 'ConFindr 0.4.2'
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
                        default=os.environ.get('CONFINDR_DB', os.path.expanduser('~/.confindr_db')),
                        help='Databases folder. If you don\'t already have databases, they will be downloaded '
                             'automatically. You may also specify the full path to the databases.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=cpu_count,
                        help='Number of threads to run analysis with.')
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
                        choices=['Illumina', 'Nanopore', 'auto'],
                        default='auto',
                        help='Type of input data. Default is to guess which type of data based on read length. '
                             'Currently has no effect, but future versions of ConFindr will support nanopore data '
                             'as well.')
    parser.add_argument('-verbosity', '--verbosity',
                        choices=['debug', 'info', 'warning'],
                        default='info',
                        help='Amount of output you want printed to the screen. Defaults to info, which should be good '
                             'for most users.')

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
    # Check for dependencies.
    logging.info('Welcome to {version}! Beginning analysis of your samples...'.format(version=version))
    all_dependencies_present = True
    # Re-enable minimap2 as dependency once nanopore stuff actually works.
    # dependencies = ['bbmap.sh', 'bbduk.sh', 'mash', 'minimap2']
    dependencies = ['bbmap.sh', 'bbduk.sh', 'mash']
    for dependency in dependencies:
        if dependency_check(dependency) is False:
            logging.error('Dependency {} not found. Please make sure it is installed and present'
                          ' on your $PATH.'.format(dependency))
            all_dependencies_present = False
    if not all_dependencies_present:
        quit(code=1)

    # Make the output directory.
    if not os.path.isdir(args.output_name):
        os.makedirs(args.output_name)

    # Check if databases necessary to run are present, and download them if they aren't
    check_for_databases_and_download(database_location=args.databases,
                                     tmpdir=args.output_name)

    # Figure out what pairs of reads, as well as unpaired reads, are present.
    paired_reads = find_paired_reads(args.input_directory, forward_id=args.forward_id, reverse_id=args.reverse_id)
    unpaired_reads = find_unpaired_reads(args.input_directory, forward_id=args.forward_id, reverse_id=args.reverse_id)
    # Process paired reads, one sample at a time.
    for pair in paired_reads:
        sample_name = os.path.split(pair[0])[-1].split(args.forward_id)[0]
        logging.info('Beginning analysis of sample {}...'.format(sample_name))
        try:
            find_contamination(pair=pair,
                               forward_id=args.forward_id,
                               threads=args.threads,
                               output_folder=args.output_name,
                               databases_folder=args.databases,
                               keep_files=args.keep_files,
                               quality_cutoff=args.quality_cutoff,
                               base_cutoff=args.base_cutoff)
        except subprocess.CalledProcessError:
            # If something unforeseen goes wrong, traceback will be printed to screen.
            # We then add the sample to the report with a note that it failed.
            multi_positions = 0
            genus = 'Error processing sample'
            write_output(output_report=os.path.join(args.output_name, 'confindr_report.csv'),
                         sample_name=sample_name,
                         multi_positions=multi_positions,
                         genus=genus,
                         percent_contam='NA',
                         contam_stddev='NA')
            logging.warning('Encountered error when attempting to run ConFindr on sample '
                            '{sample}. Skipping...'.format(sample=sample_name))
            if args.keep_files is False:
                shutil.rmtree(os.path.join(args.output_name, sample_name))
    # Process unpaired reads, also one sample at a time.
    for pair in unpaired_reads:
        sample_name = os.path.split(pair[0])[-1].split('.')[0]
        logging.info('Beginning analysis of sample {}...'.format(sample_name))
        try:
            find_contamination(pair=pair,
                               forward_id=args.forward_id,
                               threads=args.threads,
                               output_folder=args.output_name,
                               databases_folder=args.databases,
                               keep_files=args.keep_files,
                               quality_cutoff=args.quality_cutoff,
                               base_cutoff=args.base_cutoff)
        except subprocess.CalledProcessError:
            # If something unforeseen goes wrong, traceback will be printed to screen.
            # We then add the sample to the report with a note that it failed.
            multi_positions = 0
            genus = 'Error processing sample'
            write_output(output_report=os.path.join(args.output_name, 'confindr_report.csv'),
                         sample_name=sample_name,
                         multi_positions=multi_positions,
                         genus=genus,
                         percent_contam='NA',
                         contam_stddev='NA')
            logging.warning('Encountered error when attempting to run ConFindr on sample '
                            '{sample}. Skipping...'.format(sample=sample_name))
            if args.keep_files is False:
                shutil.rmtree(os.path.join(args.output_name, sample_name))
    # Process unpaired reads, also one sample at a time.
    logging.info('Contamination detection complete!')
