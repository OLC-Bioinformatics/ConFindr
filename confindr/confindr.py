#!/usr/bin/env python

from io import StringIO
import multiprocessing
import statistics
import subprocess
import argparse
import shutil
import glob
import time
import csv
import os
import pysam
from Bio import SeqIO
from biotools import mash
from biotools import bbtools
from Bio.Blast import NCBIXML
from biotools import jellyfish
from Bio.Blast.Applications import NcbiblastnCommandline
from accessoryFunctions.accessoryFunctions import printtime


def write_to_logfile(logfile, out, err, cmd):
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

    fastq_files = glob.glob(fastq_directory + '/*.f*q*')
    for name in fastq_files:
        if forward_id in name and os.path.isfile(name.replace(forward_id, reverse_id)):
            pair_list.append([name, name.replace(forward_id, reverse_id)])
    return pair_list


def run_mashsippr(sequence_dir, output_dir, database_dir, logfile=None):
    """
    Runs the mashsippr on a directory containing fastq files. Creates a report in output_dir/reports/mash.csv that
    contains information on what genus each sample in the directory is.
    :param sequence_dir: Path to a directory containing fastq files.
    :param output_dir: Path to an output directory. Created if it doesn't exist.
    :param database_dir: Path to the directory where RefSeqSketchesDefault.msh is found.
    :return: True if mash.csv was created, False if not created.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    cmd = 'python -m confindr.mashsippr -s {sequence_dir} -t {database_dir} {output_dir}'.format(sequence_dir=sequence_dir,
                                                                                                 database_dir=database_dir,
                                                                                                 output_dir=output_dir)
    if logfile:
        write_to_logfile(logfile, '', '', cmd)
        with open(logfile, 'a+') as f:
            subprocess.call(cmd, shell=True, stdout=f, stderr=f)
    else:
        with open(os.devnull, 'w') as f:
            subprocess.call(cmd, shell=True, stdout=f, stderr=f)

    if os.path.isfile(os.path.join(output_dir, 'reports/mash.csv')):
        return True
    else:
        return False


def read_mashsippr_output(mashsippr_result_file, sample):
    """
    Method to read the output of a mashsippr run. Given a sample name and a result file, finds out what genus the sample
    is from.
    :param mashsippr_result_file: Path to a mashsippr result csv.
    :param sample: Sample name, as seen in first column of the mash.csv file
    :return: A string containing the genus the sample belongs to. If the sample can't be found in the csv file, returns
    NA
    """
    genus = 'NA'
    with open(mashsippr_result_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Strain'] == sample:
                genus = row['ReferenceGenus']
    return genus


def cleanup_mashsippr_files(pairs, args):
    """
    Cleans up mashsippr files that get left behind when mashsippr runs.
    :param pairs: The paired files found using find_paired_reads.
    :param args: List of arguments.
    """
    samples = list()
    for pair in pairs:  # Get the sample name for each pair.
        samples.append(os.path.split(pair[0])[-1].split(args.forward_id)[0])
    # For each sample, find the folder left over by mashsippr and get rid of it.
    for sample in samples:
        shutil.rmtree(os.path.join(args.input_directory, sample))


def find_genusspecific_alleles(profiles_file, target_genus):
    """
    Given a genus name, will parse a custom-made profiles file to find which alleles should be excluded for a target
    genus.
    :param profiles_file: Path to profiles file.
    :param target_genus: Genus that needs to be found.
    :return: List of genes that should be excluded for the genus in question. If the specified genus can't be found
    in the profiles file, list will be empty.
    """
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
    """
    Sets up genus-specific databases (aka databases that exclude genes that are known to have multiple copies, such as
    BACT000060 and BACT000065 for Escherichia)
    :param database_folder: Folder containing rMLST database (will also end up including the genus-specfic database)
    :param genus: Name of genus that database is being created for.
    :param genes_to_exclude: List of genes to exclude for that genus database, generated by find_genusspecific_alleles()
    """
    with open(os.path.join(database_folder, '{}_db.fasta'.format(genus)), 'w') as f:
        sequences = SeqIO.parse(os.path.join(database_folder, 'rMLST_combined.fasta'), 'fasta')
        for item in sequences:
            if item.id.split('_')[0] not in genes_to_exclude:
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


def subsample_reads(forward_in, reverse_in, coverage_level, genome_size, forward_out, reverse_out,
                    threads=12, logfile=None):
    """
    Will subsample reads to a desired coverage level, given the coverage level and genome size.
    :param forward_in: Forward input reads.
    :param reverse_in: Reverse input reads.
    :param coverage_level: Desired coverage depth, as an int.
    :param genome_size: Estimated genome size, as an int.
    :param forward_out: Forward output reads.
    :param reverse_out: Reverse output reads.
    :param threads: Number of threads to use.
    """
    bases_target = coverage_level * genome_size
    out, err, cmd = bbtools.subsample_reads(forward_in, forward_out, bases_target, reverse_in=reverse_in,
                                            reverse_out=reverse_out, returncmd=True, threads=str(threads))
    if logfile:
        write_to_logfile(logfile, out, err, cmd)


def generate_kmers(forward_reads, reverse_reads, counts_file, kmer_size, tmpdir, logfile=None):
    """
    Generates a set of kmers given a set of forward and reverse reads using jellyfish.
    Output will be a fasta-formatted file, with the title of each kmer being its count.
    :param forward_reads: Path to forward input reads.
    :param reverse_reads: Path to reverse input reads.
    :param counts_file: Output fasta file.
    :param kmer_size: Kmer size, should be an int.
    :param tmpdir: Temporary directory where intermediary files get stored. Deleted when method finishes.
    """
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    out, err, cmd = jellyfish.count(forward_reads, reverse_in=reverse_reads, count_file=os.path.join(tmpdir, 'mer_counts.jf'),
                                    kmer_size=kmer_size, options='--bf-size 100M', returncmd=True)
    if logfile:
        write_to_logfile(logfile, out, err, cmd)
    jellyfish.dump(os.path.join(tmpdir, 'mer_counts.jf'), counts_file)  # TODO: Add logging for this too.
    shutil.rmtree(tmpdir)


def rename_kmers(input_kmers, output_kmers, cutoff):
    """
    Given a fasta-formatted kmer count file generated by jellyfish, renames the kmers so that
    the names for each are unique in format >count_uniqueid
    :param input_kmers: Path to Kmers created by jellyfish (the generate_kmers() method works to make these)
    :param output_kmers: Path to output location. Can be the same as input location if you want to overwrite things.
    :param cutoff: Kmers must appear at least this many times to make it into the renamed file.
    :return: Total number of unique kmers present, as an int.
    """
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


def check_db_presence(database):
    extensions = ['.nhr', '.nin', '.nsq']
    is_present = True
    for extension in extensions:
        if not os.path.isfile(database + extension):
            is_present = False
    return is_present


def make_blast_database(database, logfile='log.txt'):
    cmd = 'makeblastdb -in {} -dbtype nucl'.format(database)
    with open(logfile, 'a+') as f:
        subprocess.call(cmd, shell=True, stdout=f, stderr=f)


def present_in_db(query_sequence, database, kmer_size):
    # Blast the sequence against our database.
    blastn = NcbiblastnCommandline(db=database, outfmt=5)
    stdout, stderr = blastn(stdin=query_sequence)
    # If there's any full-length result, the sequence is present. No result means not present.
    if stdout:
        for record in NCBIXML.parse(StringIO(stdout)):
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.align_length == kmer_size:
                        return True
                    else:
                        return False
        # Sometimes despite something being in stdout there aren't any records to iterate through.
        # Not how I thought it worked, but apparently the case.
        return False
    # Given that apparently stdout always gets created, I don't think this is actually reachable,
    # but it's left here just in case I've totally misunderstood how things work.
    else:
        return False


def find_contamination(pair, args, genus):
    log = os.path.join(args.output_name, 'confindr_log.txt')
    sample_start = time.time()
    snv_list = list()
    max_kmers = 0
    # Main method for finding contamination - works on one pair at a time.
    sample_name = os.path.split(pair[0])[-1].split(args.forward_id)[0]
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
                        reverse_out=os.path.join(sample_tmp_dir, 'rmlst_R2.fastq.gz'),
                        threads=args.threads, logfile=log)
    printtime('Quality trimming...', sample_start)
    out, err, cmd = bbtools.bbduk_trim(forward_in=os.path.join(sample_tmp_dir, 'rmlst_R1.fastq.gz'),
                                       reverse_in=os.path.join(sample_tmp_dir, 'rmlst_R2.fastq.gz'),
                                       forward_out=os.path.join(sample_tmp_dir, 'trimmed_R1.fastq.gz'),
                                       reverse_out=os.path.join(sample_tmp_dir, 'trimmed_R2.fastq.gz'),
                                       threads=str(args.threads), returncmd=True)
    write_to_logfile(log, out, err, cmd)
    # Now do the actual contamination detection cycle the number of times specified by arguments.
    printtime('Beginning {} cycles of contamination detection...'.format(str(args.number_subsamples)), sample_start)
    for i in range(args.number_subsamples):
        printtime('Working on cycle {} of {}...'.format(str(i + 1), str(args.number_subsamples)), sample_start, '\033[0;35m')
        # Subsample
        subsample_reads(forward_in=os.path.join(sample_tmp_dir, 'trimmed_R1.fastq.gz'),
                        reverse_in=os.path.join(sample_tmp_dir, 'trimmed_R2.fastq.gz'),
                        coverage_level=args.subsample_depth,
                        genome_size=35000,  # This is the sum of the longest allele for each rMLST gene.
                        forward_out=os.path.join(sample_tmp_dir, 'subsample_{}_R1.fastq.gz'.format(str(i))),
                        reverse_out=os.path.join(sample_tmp_dir, 'subsample_{}_R2.fastq.gz'.format(str(i))),
                        threads=args.threads, logfile=log)
        # Kmerize subsampled reads.
        generate_kmers(forward_reads=os.path.join(sample_tmp_dir, 'subsample_{}_R1.fastq.gz'.format(str(i))),
                       reverse_reads=os.path.join(sample_tmp_dir, 'subsample_{}_R2.fastq.gz'.format(str(i))),
                       counts_file=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                       kmer_size=args.kmer_size,
                       tmpdir=os.path.join(sample_tmp_dir, 'tmp'), logfile=log)
        # Rename kmers so each has a unique ID, and count the number of kmers.
        num_kmers = rename_kmers(input_kmers=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                                 output_kmers=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                                 cutoff=args.kmer_cutoff)
        if num_kmers > max_kmers:
            max_kmers = num_kmers
        # Find mismatches.

        # Step 1 of mismatch finding: Run bbmap with the kmer file.
        out, err, cmd = bbtools.bbmap(reference=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                                      forward_in=os.path.join(sample_tmp_dir, 'kmer_counts_{}.fasta'.format(str(i))),
                                      ambig='all',
                                      overwrite='true',
                                      out_bam=os.path.join(sample_tmp_dir, 'subsample_{}.bam'.format(str(i))),
                                      threads=str(args.threads), returncmd=True)
        write_to_logfile(log, out, err, cmd)

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
        # First part of this step: Check that the blast database is actually present. If it isn't, make one.
        if not check_db_presence(sample_database):
            make_blast_database(sample_database, logfile=log)
        # Now set up the blast.
        # Create list of sequences to blast.
        to_blast = list()
        for fasta_id in fasta_ids:
            to_blast.append(mer_dict[fasta_id])
        # Setup the multiprocessing pool.
        pool = multiprocessing.Pool(processes=args.threads)
        db_list = [sample_database] * len(fasta_ids)
        kmer_size_list = [args.kmer_size] * len(fasta_ids)
        results = pool.starmap(present_in_db, zip(to_blast, db_list, kmer_size_list))
        pool.close()
        pool.join()
        snv_count = 0
        for result in results:
            if result:
                snv_count += 1
        snv_list.append(snv_count)

    # Find cross contamination.
    printtime('Finding cross contamination...', sample_start)
    genera_present = list()
    out, err, cmd = mash.screen('{}/refseq.msh'.format(args.databases), pair[0],
                                pair[1], threads=args.threads, w='', i='0.95', p=str(args.threads),
                                output_file=os.path.join(sample_tmp_dir, 'screen.tab'), returncmd=True)
    write_to_logfile(log, out, err, cmd)
    screen_output = mash.read_mash_screen(os.path.join(sample_tmp_dir, 'screen.tab'))
    for item in screen_output:
        mash_genus = item.query_id.split('/')[-3]
        if mash_genus not in genera_present:
            genera_present.append(mash_genus)
    if len(genera_present) <= 1:
        genera_present = 'NA'
        cross_contam = False
    else:
        tmpstr = ''
        for mash_genus in genera_present:
            tmpstr += mash_genus + ':'
        genera_present = tmpstr[:-1]
        cross_contam = True
    # Create contamination report.
    if statistics.median(snv_list) > 2 or cross_contam or max_kmers > 45000:
        contamination = True
    else:
        contamination = False
    with open(os.path.join(args.output_name, 'confindr_report.csv'), 'a+') as f:
        f.write('{samplename},{genus},{numcontamsnvs},{numuniquekmers},{crosscontamination},'
                '{contamstatus}\n'.format(samplename=sample_name,
                                          genus=genus,
                                          numcontamsnvs=statistics.median(snv_list),
                                          numuniquekmers=max_kmers,
                                          crosscontamination=genera_present,
                                          contamstatus=contamination))
    shutil.rmtree(sample_tmp_dir)
    printtime('Finished analysis of sample {}!'.format(sample_name), sample_start)


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
    dependencies = ['jellyfish', 'bbmap.sh', 'bbduk.sh', 'blastn', 'mash']
    for dependency in dependencies:
        if dependency_check(dependency) is False:
            print('WARNING: Dependency {} not found. ConFindr will likely crash!'.format(dependency))
    args = parser.parse_args()
    # Make the output directory.
    if not os.path.isdir(args.output_name):
        os.makedirs(args.output_name)
    # Open the output report file.
    with open(os.path.join(args.output_name, 'confindr_report.csv'), 'w') as f:
        f.write('Sample,Genus,NumContamSNVs,NumUniqueKmers,CrossContamination,ContamStatus\n')
    # Run mashsippr on all sample, and then read the file individually for each sample.
    printtime('Determining genus of each sample...', start)
    run_mashsippr(args.input_directory, args.output_name, args.databases)
    paired_reads = find_paired_reads(args.input_directory, forward_id=args.forward_id, reverse_id=args.reverse_id)
    cleanup_mashsippr_files(paired_reads, args)
    for pair in paired_reads:
        sample_name = os.path.split(pair[0])[-1].split(args.forward_id)[0]
        print('\n\n')
        printtime('Beginning analysis of sample {}...\n'.format(sample_name), start, '\033[1;34m')
        genus = read_mashsippr_output(os.path.join(args.output_name, 'reports/mash.csv'), sample_name)
        find_contamination(pair, args, genus)
    printtime('Contamination detection complete!', start, '\033[0;32m')
