#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import printtime, GenObject, make_path, run_subprocess, write_to_logfile
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from glob import glob
import multiprocessing
import pysam
import statistics
from subprocess import TimeoutExpired
__author__ = 'adamkoziol'


class PipelineContaminationDetection(object):

    def main(self):
        """

        """
        printtime('Calculating contamination in reads', self.start)
        self.extract_rmlst_reads()
        for i in range(self.number_subsamples):
            printtime('Working on subsampling {} of {}'.format(str(i + 1), self.number_subsamples), self.start)
            self.subsample_reads()
            self.run_jellyfish()
            self.write_mer_file()
            self.run_bbmap()
            self.read_bamfile()
        printtime('Writing output and cleaning up temporary files!', self.start)
        self.write_output()

    def extract_rmlst_reads(self):
        """
        rMLST read extraction. Should be the first thing called after parsing the fastq directory.
        """
        for sample in self.metadata:
            # Create the object to store the variables
            setattr(sample, self.analysistype, GenObject())
            # Initialise variables
            sample[self.analysistype].snv_count = list()
            # Initialise a starting value for the number of unique kmers found in each sample
            sample[self.analysistype].unique_kmers = -1
            # Set and create the output directory
            try:
                sample[self.analysistype].outputdir = os.path.join(sample.run.outputdirectory, self.analysistype)
            except KeyError:
                sample[self.analysistype].outputdir = os.path.join(sample.general.outputdirectory, self.analysistype)
            make_path(sample[self.analysistype].outputdir)
            sample[self.analysistype].logout = os.path.join(sample[self.analysistype].outputdir, 'logout.txt')
            sample[self.analysistype].logerr = os.path.join(sample[self.analysistype].outputdir, 'logerr.txt')
            sample[self.analysistype].baitedfastq = os.path.join(
                sample[self.analysistype].outputdir, '{}_targetMatches.fastq.gz'.format(self.analysistype))
            # Create the command to run the baiting - paired inputs and a single, zipped output
            sample[self.analysistype].bbdukcmd = 'bbduk.sh ref={} in1={} in2={} threads={} outm={}'\
                .format(self.database,
                        sample.general.trimmedcorrectedfastqfiles[0],
                        sample.general.trimmedcorrectedfastqfiles[1],
                        str(self.threads),
                        sample[self.analysistype].baitedfastq)
            # Sometimes bbduk hangs forever, so that needs to be handled. Give it a very generous timeout.
            try:
                # Run the call, and write any errors to the logfile
                command = sample[self.analysistype].bbdukcmd
                if self.analyse:
                    out, err = run_subprocess(command)
                else:
                    out = str()
                    err = str()
                write_to_logfile(command, command, self.logfile, sample.general.logout, sample.general.logerr,
                                 sample[self.analysistype].logout, sample[self.analysistype].logerr)
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                                 sample[self.analysistype].logout, sample[self.analysistype].logerr)
            except TimeoutExpired:
                print('ERROR: Could not extract rMLST reads from sample {}'.format(sample.name))

    def subsample_reads(self):
        """
        Subsampling of reads to 20X coverage of rMLST genes (roughly).
        To be called after rMLST extraction and read trimming, in that order.
        """
        for sample in self.metadata:
            # Create the name of the subsampled read file
            sample[self.analysistype].subsampledreads = os.path.join(
                sample[self.analysistype].outputdir, '{}_targetMatches_subsampled.fastq'.format(self.analysistype))
            # Set the reformat.sh command - as this command will be run multiple times, overwrite previous iterations
            # each time. Use samplebasestarget to provide an approximation of the number of bases to include in the
            # subsampled reads e.g. for rMLST: 700000 (approx. 35000 bp total length of genes x 20X coverage)
            sample[self.analysistype].subsamplecmd = 'reformat.sh in={} out={} overwrite samplebasestarget={}' \
                .format(sample[self.analysistype].baitedfastq,
                        sample[self.analysistype].subsampledreads,
                        self.samplebasestarget)
            # Run the call, and write any errors to the logfile
            command = sample[self.analysistype].subsamplecmd
            if self.analyse:
                out, err = run_subprocess(command)
            else:
                out = str()
                err = str()
            write_to_logfile(command, command, self.logfile, sample.general.logout, sample.general.logerr,
                             sample[self.analysistype].logout, sample[self.analysistype].logerr)
            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                             sample[self.analysistype].logout, sample[self.analysistype].logerr)

    def run_jellyfish(self):
        """
        Runs jellyfish to split subsample reads into kmers. Runs kmers through a bloom filter to get rid of singletons
        that are likely just sequencing errors. Should be run after subsampling reads.
        """
        for sample in self.metadata:
            # Set the name of the jellyfish count file
            sample[self.analysistype].jellyfish_file = os.path.join(sample[self.analysistype].outputdir,
                                                                    sample.name + '_jellyfish')
            # Set the system call
            sample[self.analysistype].jellyfishcountcmd \
                = 'jellyfish count -m 31 -s 100M --bf-size 100M -C -F 2 {} -o {} -t {}'\
                .format(sample[self.analysistype].subsampledreads,
                        sample[self.analysistype].jellyfish_file,
                        str(self.threads))
            # Run the call, and write any errors to the logfile
            command = sample[self.analysistype].jellyfishcountcmd
            if self.analyse:
                out, err = run_subprocess(command)
            else:
                out = str()
                err = str()
            write_to_logfile(command, command, self.logfile, sample.general.logout, sample.general.logerr,
                             sample[self.analysistype].logout, sample[self.analysistype].logerr)
            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                             sample[self.analysistype].logout, sample[self.analysistype].logerr)

    def write_mer_file(self):
        """
        Writes the mer file created by jellyfish to a fasta format to be used by other things downstream.
        Only writes kmers that have been seen at least twice to attempt to get rid of sequencing erros.
        """
        for sample in self.metadata:
            # Set the name of the kmer file dumped from jellyfish
            sample[self.analysistype].mer_fasta = sample[self.analysistype].jellyfish_file + '.fasta'
            sample[self.analysistype].solid_mers = sample[self.analysistype].jellyfish_file + '_solid.fasta'
            # Set the system call
            sample[self.analysistype].jellyfishdumpcmd =\
                'jellyfish dump {} > {}'\
                .format(sample[self.analysistype].jellyfish_file,
                        sample[self.analysistype].mer_fasta)
            # Run the system call
            command = sample[self.analysistype].jellyfishdumpcmd
            if self.analyse:
                out, err = run_subprocess(command)
            else:
                out = str()
                err = str()
            write_to_logfile(command, command, self.logfile, sample.general.logout, sample.general.logerr,
                             sample[self.analysistype].logout, sample[self.analysistype].logerr)
            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                             sample[self.analysistype].logout, sample[self.analysistype].logerr)
            # Read in the dumped file to a list
            with open(sample[self.analysistype].mer_fasta, 'r') as mers:
                fastas = mers.readlines()
            # Initialise variables for use in parsing outputs
            num_mers = 0
            sequences = list()
            # Iterate through the list of the fasta outputs. Output is a multifasta e.g.:
            # >8
            # GCCTGGAAAACTGGCCACCGGCAAGCATCGC
            # where the header, >8, indicates that the sequence is present 8 times in the sample
            for i in range(len(fastas)):
                # Find the headers
                if '>' in fastas[i]:
                    # If the number of times the sequence is present in the sample is greater than one, increment
                    # the total number of kmers observed
                    if int(fastas[i].replace('>', '')) > 1:
                        num_mers += 1
                        # Append a string of the header plus the total number of mers, and the sequence information
                        # to the list of sequences e.g. ['>8_1\nGCCTGGAAAACTGGCCACCGGCAAGCATCGC\n']
                        sequences.append(fastas[i].rstrip() + '_' + str(num_mers) + '\n' + fastas[i + 1])
            # Write out our solid kmers to file to be used later.
            with open(sample[self.analysistype].solid_mers, 'w') as solidmers:
                solidmers.write(''.join(sequences))
            # Update the number of unique kmers
            if num_mers > sample[self.analysistype].unique_kmers:
                sample[self.analysistype].unique_kmers = num_mers

    def run_bbmap(self):
        """
        Runs bbmap on kmer fasta file, against kmer fasta file to generate a samfile which can then be parsed to find
        low frequency kmers that have one mismatch to high frequency kmers, indicating that they're from contaminating
        alleles.
        """
        for sample in self.metadata:
            # Create the name for the output bam file
            sample[self.analysistype].bamfile = sample[self.analysistype].mer_fasta.replace('.fasta', '.bam')
            # Set the bbmap call - use the overwrite option to overwrite previous files that were created on previous
            # iterations, ambig=all to use all highest scoring mappings, nodisk to build index in memory, and only write
            # ouput to disk, local to allow soft-clipping
            sample[self.analysistype].bbmapcmd = \
                'bbmap.sh ref={} in={} outm={} overwrite ambig=all nodisk local threads={}'\
                .format(sample[self.analysistype].solid_mers,
                        sample[self.analysistype].solid_mers,
                        sample[self.analysistype].bamfile,
                        str(self.threads))
            # Run the call, and write any errors to the logfile
            command = sample[self.analysistype].bbmapcmd
            if self.analyse:
                out, err = run_subprocess(command)
            else:
                out = str()
                err = str()
            write_to_logfile(command, command, self.logfile, sample.general.logout, sample.general.logerr,
                             sample[self.analysistype].logout, sample[self.analysistype].logerr)
            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                             sample[self.analysistype].logout, sample[self.analysistype].logerr)

    def read_bamfile(self):
        """
        Reads the bamfile generated by run_bbmap to find low frequency kmers that are mismatches of high frequency kmers
        BLASTS the kmers found against the rMLST database in order to ensure that the kmers come from rMLST genes,
        and weren't overhangs from rMLST genes into other stuff from the original reads.
        """
        # Make sure the blast database is setup. If it's already there, this won't do anything.
        self.make_db()
        for sample in self.metadata:
            # Create a list of sequences to be processed with BLAST
            to_blast = list()
            # Make a dict for all or our mer sequences so we know what sequence goes with which header
            mer_dict = dict()
            # Open the jellyfish dumped fasta file
            with open(sample[self.analysistype].solid_mers, 'r') as mers:
                # Read the lines into a list
                mers = mers.readlines()
            #
            for i in range(0, len(mers), 2):
                key = mers[i].replace('>', '')
                key = key.replace('\n', '')
                mer_dict[key] = mers[i + 1]
            i = 0
            try:
                # Read in the sam file.
                bamfile = pysam.AlignmentFile(sample[self.analysistype].bamfile, 'rb')
                for match in bamfile:
                    # Every time we get a match, put the sequence into a list of things we're going to BLAST to
                    # make sure they actually come from rMLST genes.
                    if "1X" in match.cigarstring and match.query_alignment_length == 31:
                        reference = bamfile.getrname(match.reference_id)
                        to_blast.append(mer_dict[reference])
                # BLAST our potentially contaminating sequences in parallel to speed things up because BLASTing
                # is painfully slow.
                pool = multiprocessing.Pool(processes=self.threads)
                results = pool.starmap(present_in_db, zip(to_blast, [self.database]))
                pool.close()
                pool.join()
                # Every time BLAST verifies that the potential contaminant is indeed part of our rMLST dataset,
                # increment our contamination counter.
                for item in results:
                    if item:
                        i += 1
                # Append the result to our snv_count list. We'll take the median of the list later.
                sample[self.analysistype].snv_count.append(i)
            except OSError:
                # If our samfile is empty because there were no rMLST genes, just say that nothing was found.
                sample[self.analysistype].snv_count.append(0)

    def make_db(self):
        """
        Makes the blast database if it isn't present. Doesn't do anything if we already have database files.
        """
        db_files = ['.nhr', '.nin', '.nsq']
        db_present = True
        for db_file in db_files:
            if not os.path.isfile(self.database + db_file):
                db_present = False
        if not db_present:
            printtime('Making database!', self.start)
            command = 'makeblastdb -dbtype nucl -in ' + self.database
            if self.analyse:
                out, err = run_subprocess(command)
            else:
                out = str()
                err = str()
            write_to_logfile(command, command, self.logfile, None, None, None, None)
            write_to_logfile(out, err, self.logfile, None, None, None, None)

    def write_output(self):
        """
        Writes our output. Takes the median number of contaminating SNVs across our subsample reps, and the highest
        number of unique kmers.
        """
        with open(self.reportfile, 'w') as report:
            header = 'Sample,NumContamSNVs,NumUniqueKmers,ContamStatus\n'
            report.write(header)
            for sample in self.metadata:
                with open(os.path.join(sample[self.analysistype].outputdir, 'confinder_results.csv'), 'w') as rep:
                    rep.write(header)
                    try:
                        snv_median = statistics.median(sample[self.analysistype].snv_count)
                    except statistics.StatisticsError:
                        snv_median = 0
                    if snv_median > 0 or sample[self.analysistype].unique_kmers > 50000:
                        sample[self.analysistype].contam_status = 'Contaminated'
                    else:
                        sample[self.analysistype].contam_status = 'Clean'
                    data = '{},{},{},{}\n'.format(sample.name,
                                                  str(statistics.median(sample[self.analysistype].snv_count)),
                                                  str(sample[self.analysistype].unique_kmers),
                                                  sample[self.analysistype].contam_status)
                    report.write(data)
                    rep.write(data)

    def __init__(self, inputobject, samplebasestarget=700000):
        self.metadata = inputobject.runmetadata.samples
        self.database = glob(os.path.join(inputobject.reffilepath, 'rMLST', '*.fasta'))[0]
        self.logfile = inputobject.logfile
        self.threads = inputobject.cpus
        self.analysistype = 'confinder'
        self.number_subsamples = 5
        self.start = inputobject.starttime
        self.reportpath = inputobject.reportpath
        make_path(self.reportpath)
        self.samplebasestarget = samplebasestarget
        self.reportfile = os.path.join(self.reportpath, self.analysistype + '.csv')
        if not os.path.isfile(self.reportfile):
            self.analyse = True
        else:
            self.analyse = False
        self.main()


def present_in_db(query_sequence, database):
    """
    Checks if a sequence is present in our rMLST database, as some overhangs on reads can be in repetitive regions
    that could cause false positives and should therefore be screened out.
    :param query_sequence: nucleotide sequence, as a string.
    :param database: rMLST database to be searched
    :return: True if sequence is found in rMLST database, False if it isn't.
    """
    # Blast the sequence against our database.
    blastn = NcbiblastnCommandline(db=database, outfmt=6)
    stdout, stderr = blastn(stdin=query_sequence)
    # If there's any result, the sequence is present. No result means not present.
    if stdout:
        return True
    else:
        return False
