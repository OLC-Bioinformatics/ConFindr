import statistics
from accessoryFunctions.accessoryFunctions import printtime
import time
import multiprocessing
import shutil
import pysam
import subprocess
import argparse
import glob
import os
from Bio.Blast.Applications import NcbiblastnCommandline


class ContamObject:
    """
    Object that keeps track of everything related to contamination detection for a sample.
    Created by telling it the name of a path which contains forward reads (or, if single ended sequence, the only
    set of reads)"
    """
    def __init__(self, forward_reads):
        self.forward_reads = forward_reads
        self.reverse_reads = 'NA'
        self.unique_kmers = -1
        self.snv_count = list()
        self.trimmed_forward_reads = 'NA'
        self.trimmed_reverse_reads = 'NA'
        self.forward_rmlst_reads = 'NA'
        self.reverse_rmlst_reads = 'NA'
        self.subsampled_forward = 'NA'
        self.subsampled_reverse = 'NA'
        self.jellyfish_file = 'NA'
        self.mer_fasta = 'NA'
        self.samfile = 'NA'

#  Would be nice to do error checking of some sort on input data to make sure they're properly formatted fastqs.


class Detector(object):
    def __init__(self, args):
        self.fastq_directory = args.fastq_directory
        self.samples = dict()
        self.tmpdir = args.output_name + 'tmp/'
        self.outfile = args.output_name + '.csv'
        self.logfile = args.output_name + '.log'
        self.database = args.rmlst_database
        if not os.path.exists(self.tmpdir):
            os.makedirs(self.tmpdir)
        self.threads = args.threads

    def parse_fastq_directory(self):
        """
        Parses the fastq directory specified in self.fastq_directory (which is passed in to the object by the
        argparser). Creates ContamObjects for each sample it finds, which are put into the self.samples dict.
        Keys for the dictionary are the paths to forward reads.
        """
        fastq_files = glob.glob(self.fastq_directory + '/*.f*q*')
        for name in fastq_files:
            name = os.path.abspath(name)
            if "_R1" in name and os.path.isfile(name.replace("_R1", "_R2")):
                self.samples[name] = ContamObject(name)
                self.samples[name].reverse_reads = name.replace('_R1', '_R2')
            # Other naming convention support.
            elif "_1" in name and os.path.isfile(name.replace("_1", "_2")):
                self.samples[name] = ContamObject(name)
                self.samples[name].reverse_reads = name.replace('_R1', '_R2')
            # Assume that if we can't find a mate reads are single ended, and add them to the appropriate list.
            elif '_R2' not in name and '_2' not in name:
                self.samples[name] = ContamObject(name)

    def quality_trim_reads(self):
        """
        To be called after extraction of rMLST sequences. Will quality trim the reads, which is necessary before
        subsampling occurs. Uses bbduk (tested with version 37.23, although others should work) to trim reads, which
        helps with excluding false positives quire a bit.
        """
        # Figure out where bbduk is so that we can use the adapter file.
        cmd = 'which bbduk.sh'
        bbduk_dir = subprocess.check_output(cmd.split()).decode('utf-8')
        bbduk_dir = bbduk_dir.split('/')[:-1]
        bbduk_dir = '/'.join(bbduk_dir)
        for sample in self.samples:
            # This captures paired reads.
            if self.samples[sample].reverse_reads != 'NA':
                self.samples[sample].trimmed_forward_reads = self.samples[sample].forward_rmlst_reads + '_trimmed'
                self.samples[sample].trimmed_reverse_reads = self.samples[sample].reverse_rmlst_reads + '_trimmed'
                cmd = 'bbduk.sh in1={} in2={} out1={} out2={} qtrim=w trimq=20 k=25 minlength=50 forcetrimleft=15' \
                    ' ref={}/resources/adapters.fa overwrite hdist=1 tpe tbo threads={}'.format(self.samples[sample].forward_rmlst_reads,
                                                                                   self.samples[sample].reverse_rmlst_reads,
                                                                                   self.samples[sample].trimmed_forward_reads,
                                                                                   self.samples[sample].trimmed_reverse_reads,
                                                                                   bbduk_dir, str(self.threads))
            else:
                # Unpaired reads.
                self.samples[sample].trimmed_forward_reads = self.samples[sample].forward_rmlst_reads + '_trimmed'
                cmd = 'bbduk.sh in={} out={} qtrim=w trimq=20 k=25 minlength=50 forcetrimleft=15' \
                      ' ref={}/resources/adapters.fa overwrite hdist=1 tpe tbo threads={}'.format(self.samples[sample].forward_rmlst_reads,
                                                                                       self.samples[sample].trimmed_forward_reads,
                                                                                       bbduk_dir, str(self.threads))
            # For some reason bbduk will occasionally hang and run forever, so we need to handle that.
            try:
                with open(self.logfile, 'a+') as logfile:
                    subprocess.call(cmd, shell=True, stderr=logfile, timeout=600)
            except subprocess.TimeoutExpired:
                # If BBDUK does hang forever, remove sample from further analysis.
                print('ERROR: Sample {} was not able to be quality trimmed. Excluding from further analyses...'.format(sample))
                del self.samples[sample]

    def extract_rmlst_reads(self):
        """
        rMLST read extraction. Should be the first thing called after parsing the fastq directory.
        """
        for sample in self.samples:
            # For paired reads.
            if self.samples[sample].reverse_reads != 'NA':
                self.samples[sample].forward_rmlst_reads = self.tmpdir + self.samples[sample].forward_reads.split('/')[-1] + '_rmlst'
                self.samples[sample].reverse_rmlst_reads = self.tmpdir + self.samples[sample].reverse_reads.split('/')[-1] + '_rmlst'
                cmd = 'bbduk.sh ref={} in1={} in2={} outm={} outm2={} threads={}'.format(self.database,
                                                                              self.samples[sample].forward_reads,
                                                                              self.samples[sample].reverse_reads,
                                                                              self.samples[sample].forward_rmlst_reads,
                                                                              self.samples[sample].reverse_rmlst_reads,
                                                                                         str(self.threads))
            else:
                # Unpaired reads.
                self.samples[sample].forward_rmlst_reads = self.tmpdir + self.samples[sample].forward_reads.split('/')[-1] + '_rmlst'
                cmd = 'bbduk.sh ref={} in={} outm={} threads={}'.format(self.database,
                                                                        self.samples[sample].forward_reads,
                                                                        self.samples[sample].forward_rmlst_reads,
                                                                        str(self.threads))
            # Again, sometimes bbduk hangs forever, so that needs to be handled. Give it a very generous timeout.
            try:
                with open(self.logfile, 'a+') as logfile:
                    subprocess.call(cmd, shell=True, stderr=logfile, timeout=1000)
            except subprocess.TimeoutExpired:
                # If BBDUK does hang forever, remove sample from further analysis.
                print('ERROR: Could not extract rMLST reads from sample {}. Excluding from further analyses...'.format(sample))
                del self.samples[sample]

    def subsample_reads(self):
        """
        Subsampling of reads to 20X coverage of rMLST genes (roughly).
        To be called after rMLST extraction and read trimming, in that order.
        """
        for sample in self.samples:
            if self.samples[sample].trimmed_reverse_reads != 'NA':
                self.samples[sample].subsampled_forward = self.samples[sample].trimmed_forward_reads + '_subsampled'
                self.samples[sample].subsampled_reverse = self.samples[sample].trimmed_reverse_reads + '_subsampled'
                cmd = 'reformat.sh in1={} in2={} out1={} out2={} overwrite samplebasestarget=700000'.format(self.samples[sample].trimmed_forward_reads,
                                                                                                            self.samples[sample].trimmed_reverse_reads,
                                                                                                            self.samples[sample].subsampled_forward,
                                                                                                            self.samples[sample].subsampled_reverse)
            else:
                self.samples[sample].subsampled_forward = self.samples[sample].trimmed_forward_reads + '_subsampled'
                cmd = 'reformat.sh in={} out={} overwrite samplebasestarget=700000'.format(self.samples[sample].trimmed_forward_reads,
                                                                                           self.samples[sample].subsampled_forward)
            with open(self.logfile, 'a+') as logfile:
                subprocess.call(cmd, shell=True, stderr=logfile)

    def run_jellyfish(self):
        """
        Runs jellyfish to split subsample reads into kmers. Runs kmers through a bloom filter to get rid of singletons
        that are likely just sequencing errors. Should be run after subsampling reads.
        """
        for sample in self.samples:
            self.samples[sample].jellyfish_file = self.tmpdir + self.samples[sample].forward_reads.split('/')[-1] + '_jellyfish'
            if self.samples[sample].subsampled_reverse != 'NA':
                cmd = 'jellyfish count -m 31 -s 100M --bf-size 100M -C -F 2 {} {} -o {} -t {}'.format(self.samples[sample].subsampled_forward,
                                                                                                      self.samples[sample].subsampled_reverse,
                                                                                                      self.samples[sample].jellyfish_file,
                                                                                                      str(self.threads))
            else:
                cmd = 'jellyfish count -m 31 -s 100M --bf-size 100M -C -F 1 {} -o {} -t {}'.format(self.samples[sample].subsampled_forward,
                                                                                                   self.samples[sample].jellyfish_file,
                                                                                                   str(self.threads))
            with open(self.logfile, 'a+') as logfile:
                subprocess.call(cmd, shell=True, stderr=logfile)

    def write_mer_file(self):
        """
        Writes the mer file created by jellyfish to a fasta format to be used by other things downstream.
        Only writes kmers that have been seen at least twice to attempt to get rid of sequencing erros.
        """
        for sample in self.samples:
            self.samples[sample].mer_fasta = self.samples[sample].jellyfish_file + '.fasta'
            cmd = 'jellyfish dump {} > {}'.format(self.samples[sample].jellyfish_file,
                                                  self.samples[sample].mer_fasta)
            subprocess.call(cmd, shell=True)

            f = open(self.samples[sample].mer_fasta)
            fastas = f.readlines()
            f.close()

            num_mers = 0
            sequences = list()
            for i in range(len(fastas)):
                if '>' in fastas[i]:
                    if int(fastas[i].replace('>', '')) > 1:
                        num_mers += 1
                        sequences.append(fastas[i].rstrip() + '_' + str(num_mers) + '\n' + fastas[i + 1])
            # Write out our solid kmers to file to be used later.
            f = open(self.samples[sample].mer_fasta, 'w')
            f.write(''.join(sequences))
            f.close()
            if num_mers > self.samples[sample].unique_kmers:
                self.samples[sample].unique_kmers = num_mers

    def run_bbmap(self):
        """
        Runs bbmap on kmer fasta file, against kmer fasta file to generate a samfile which can then be parsed to find
        low frequency kmers that have one mismatch to high frequency kmers, indicating that they're from contaminating
        alleles.
        """
        for sample in self.samples:
            self.samples[sample].samfile = self.samples[sample].mer_fasta.replace('.fasta', '.sam')
            cmd = 'bbmap.sh ref={} in={} outm={} overwrite ambig=all ' \
                  'nodisk threads={}'.format(self.samples[sample].mer_fasta, self.samples[sample].mer_fasta,
                                             self.samples[sample].samfile, str(self.threads))
            with open(self.logfile, 'a+') as logfile:
                subprocess.call(cmd, shell=True, stderr=logfile)

    def read_samfile(self):
        """
        Reads the samfile generated by run_bbmap to find low frequency kmers that are mismatches of high frequency kmers
        BLASTS the kmers found against the rMLST database in order to ensure that the kmers come from rMLST genes,
        and weren't overhangs from rMLST genes into other stuff from the original reads.
        """
        # Make sure the blast database is setup. If it's already there, this won't do anything.
        self.make_db()
        # Make a dict for all or our mer sequences so we know what sequence goes with which header.
        for sample in self.samples:
            to_blast = list()
            f = open(self.samples[sample].mer_fasta)
            mers = f.readlines()
            f.close()
            mer_dict = dict()
            for i in range(0, len(mers), 2):
                key = mers[i].replace('>', '')
                key = key.replace('\n', '')
                mer_dict[key] = mers[i + 1]
            i = 0
            try:
                # Read in the sam file.
                samfile = pysam.AlignmentFile(self.samples[sample].samfile, 'r')
                for match in samfile:
                    # Every time we get a match, put the sequence into a list of things we're going to BLAST to
                    # make sure they actually come from rMLST genes.
                    if "1X" in match.cigarstring and match.query_alignment_length == 31:
                        reference = samfile.getrname(match.reference_id)
                        to_blast.append(mer_dict[reference])
                # BLAST our potentially contaminating sequences in parallel to speed things up because BLASTing
                # is painfully slow.
                pool = multiprocessing.Pool(processes=self.threads)
                results = pool.map(self.present_in_db, to_blast)
                pool.close()
                pool.join()
                # Every time BLAST verifies that the potential contaminant is indeed part of our rMLST dataset,
                # increment our contamination counter.
                for item in results:
                    if item:
                        i += 1
                # Append the result to our snv_count list. We'll take the median of the list later.
                self.samples[sample].snv_count.append(i)
            except OSError:
                # If our samfile is empty because there were no rMLST genes, just say that nothing was found.
                self.samples[sample].snv_count.append(0)

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
            print('Making database!')
            cmd = 'makeblastdb -dbtype nucl -in ' + self.database
            with open(self.logfile, 'a+') as logfile:
                subprocess.call(cmd, shell=True, stderr=logfile)

    def present_in_db(self, query_sequence):
        """
        Checks if a sequence is present in our rMLST database, as some overhangs on reads can be in repetitive regions
        that could cause false positives and should therfore be screened out.
        :param query_sequence: nucleotide sequence, as a string.
        :return: True if sequence is found in rMLST database, False if it isn't.
        """
        # Blast the sequence against our database.
        blastn = NcbiblastnCommandline(db=self.database, outfmt=6)
        stdout, stderr = blastn(stdin=query_sequence)
        # If there's any result, the sequence is present. No result means not present.
        if stdout:
            return True
        else:
            return False

    def write_output(self):
        """
        Writes our output. Takes the median number of contaminating SNVs across our subsample reps, and the highest
        number of unique kmers.
        """
        f = open(self.outfile, 'w')
        f.write('Sample,NumContamSNVs,NumUniqueKmers,ContamStatus\n')
        for sample in self.samples:
            try:
                snv_median = statistics.median(self.samples[sample].snv_count)
            except statistics.StatisticsError:
                snv_median = 0
            if snv_median > 0 or self.samples[sample].unique_kmers > 50000:
                contam_status = 'Contaminated'
            else:
                contam_status = 'Clean'
            f.write('{},{},{},{}\n'.format(sample.split('/')[-1], str(statistics.median(self.samples[sample].snv_count)), str(self.samples[sample].unique_kmers),
                                                                      contam_status))
        f.close()

    def cleanup(self):
        """
        Gets rid of our tmp directory and all its files.
        """
        shutil.rmtree(self.tmpdir)


if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_directory", help="Folder that contains fastq files you want to check for contamination. "
                                                "Will find any fastq file that contains .fq or .fastq in the filename.")
    parser.add_argument('output_name', help='Base name for output/temporary directories.')
    parser.add_argument('rmlst_database', help='rMLST database, in fasta format.')
    cpu_count = multiprocessing.cpu_count()
    parser.add_argument('-t', '--threads', type=int, default=cpu_count, help='Number of threads to run analysis with.')
    parser.add_argument('-n', '--number_subsamples', type=int, default=5, help='Number of time to subsample.')
    arguments = parser.parse_args()
    detector = Detector(arguments)
    Detector.parse_fastq_directory(detector)
    printtime('Extracting rMLST reads...', start)
    Detector.extract_rmlst_reads(detector)
    printtime('Quality trimming reads...', start)
    Detector.quality_trim_reads(detector)
    for i in range(arguments.number_subsamples):
        printtime('Working on subsampling {} of {}'.format(str(i + 1), arguments.number_subsamples), start)
        Detector.subsample_reads(detector)
        Detector.run_jellyfish(detector)
        Detector.write_mer_file(detector)
        Detector.run_bbmap(detector)
        Detector.read_samfile(detector)
    printtime('Writing output and cleaning up temporary files!', start)
    Detector.write_output(detector)
    Detector.cleanup(detector)
    printtime('Contamination detection complete! :D', start)
