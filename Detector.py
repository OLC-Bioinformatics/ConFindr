from detection import ContamDetect
import argparse
import time
import multiprocessing
from accessoryFunctions.accessoryFunctions import printtime
import os
import shutil

if __name__ == '__main__':

    # Check the number of CPUs available on the system to be used by bbmap/other multithreaded things.
    cpu_count = multiprocessing.cpu_count()
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_folder", help="Folder that contains fastq files you want to check for contamination. "
                                             "Will find any fastq file that contains .fq or .fastq in the filename.")
    parser.add_argument("output_file", help="Base name of the output csv you want to create. (.csv extension is added"
                                            "by the program).")
    parser.add_argument("-k", "--kmer_size", type=int, default=31, help="Size of kmer to use. Experimental feature. "
                                                                        "Probably don't mess with it.")
    parser.add_argument("-t", "--threads", type=int, default=cpu_count, help="Number of CPUs to run analysis on."
                                                                             " Defaults to number of CPUs on the system.")
    parser.add_argument("-c", "--classify", default=False, action='store_true', help='If cross-contamination is suspected, '
                                                                                     'try to classify using CLARK-l. '
                                                                                     'Off by default.')
    parser.add_argument('-tr', '--trim_reads', default=False, action='store_true', help='If enabled, trims reads to'
                                                                                        'remove low quality bases '
                                                                                        'before kmer-izing. Off by '
                                                                                        'default, but highly recommended.')
    parser.add_argument('-x', '--remove_bad_reads', default=False, action='store_true', help='Create cleaned up fastq'
                                                                                             'files with contaminants'
                                                                                             'removed. Still a work'
                                                                                             'in progress.')
    parser.add_argument('-d', '--database', default='database.fasta', type=str, help='rMLST database, minus BACT000014,'
                                                                                       'BACT000060, and BACT000065, which'
                                                                                       'tend to have multiple copies in'
                                                                                       'some species.')
    arguments = parser.parse_args()
    # Get our contamination detector object going.
    detector = ContamDetect(arguments, start)
    # Get lists of single and paired files.
    paired_files, single_files = ContamDetect.parse_fastq_directory(arguments.fastq_folder)
    # If we're trimming reads, do that for each file and put the trimmed reads into a tmp folder.
    # Then, parse the tmp folder to get new lists of filenames to do work on.
    printtime('Extracting rMLST reads...', start)
    ContamDetect.extract_rmlst_reads(detector, paired_files, single_files)
    paired_files, single_files = ContamDetect.parse_fastq_directory(arguments.output_file + 'rmlsttmp/')
    if arguments.trim_reads:
        printtime('Trimming input fastq files...', start)
        for pair in paired_files:
            for i in range(len(pair)):
                pair[i] = os.path.abspath(pair[i])
        for single in single_files:
            single = os.path.abspath(single)
        ContamDetect.trim_fastqs(detector, paired_files, single_files)
        paired_files, single_files = ContamDetect.parse_fastq_directory(arguments.output_file + 'tmp/')
    else:
        shutil.rmtree(arguments.output_file + 'tmp/')
        shutil.copytree(arguments.output_file + 'rmlsttmp/', arguments.output_file + 'tmp')
        paired_files, single_files = ContamDetect.parse_fastq_directory(arguments.output_file + 'tmp/')

    # Subsample reads to desired coverage level (20X or so, hopefully).
    printtime('Subsampling reads...', start)
    ContamDetect.subsample_reads(detector, paired_files, single_files)
    paired_files, single_files = ContamDetect.parse_fastq_directory(arguments.output_file + 'tmp/')
    # Get a counter started so that we can tell the user how far along we are.
    sample_num = 1
    # Do contamination detection on paired files first.
    for pair in paired_files:
        # Make sure paths are absolute, otherwise bad stuff tends to happen.
        for i in range(len(pair)):
            pair[i] = os.path.abspath(pair[i])
        printtime('Working on sample ' + str(sample_num) + ' of ' + str(len(paired_files) + len(single_files)), start)
        # Run jellyfish to split into mers.
        printtime('Running jellyfish...', start)
        ContamDetect.run_jellyfish(detector, pair, arguments.threads)
        # Write the mers to a file that can be used by bbmap.
        printtime('Writing mers to file...', start)
        num_mers = ContamDetect.write_mer_file(detector, arguments.output_file + 'tmp/mer_counts.jf', pair)
        # Run bbmap on the output mer file
        printtime('Finding mismatching mers...', start)
        ContamDetect.run_bbmap(detector, pair, arguments.threads)
        # Read through bbmap's samfile output to generate our statistics.
        printtime('Generating contamination statistics...', start)
        ContamDetect.read_samfile(detector, num_mers, pair)
        sample_num += 1
    # Essentially the exact same as our paired file parsing.
    for single in single_files:
        # Make sure paths are absolute, otherwise bad stuff tends to happen.
        single = os.path.abspath(single)
        printtime('Working on sample ' + str(sample_num) + ' of ' + str(len(paired_files) + len(single_files)), start)
        # Split reads into mers.
        printtime('Running jellyfish...', start)
        ContamDetect.run_jellyfish(detector, [single], arguments.threads)
        # Write mers to fasta file
        printtime('Writing mers to file...', start)
        num_mers = ContamDetect.write_mer_file(detector, arguments.output_file + 'tmp/mer_counts.jf', [single])
        # Run bbmap on fasta file.
        printtime('Finding mismatching mers...', start)
        ContamDetect.run_bbmap(detector, [single], arguments.threads)
        # Read samfile to generate statistics.
        printtime('Generating contamination statistics...', start)
        ContamDetect.read_samfile(detector, num_mers, [single])
        sample_num += 1

    end = time.time()
    shutil.rmtree(arguments.output_file + 'tmp')
    shutil.rmtree(arguments.output_file + 'rmlsttmp')
    printtime("Finished contamination detection!", start)
