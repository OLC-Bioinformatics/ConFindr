#!/usr/bin/env python
import statistics
import csv
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
from Bio.Blast import NCBIXML
from Bio import SeqIO
from io import StringIO
from biotools import mash


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
        self.genus = 'NA'
        self.genus_database = 'NA'
        self.cross = 'NA'


class Detector(object):
    def __init__(self, args):
        self.fastq_directory = args.fastq_directory
        self.samples = dict()
        self.tmpdir = args.output_name + 'tmp/'
        self.outfile = args.output_name + '.csv'
        self.logfile = args.output_name + '.log'
        if not os.path.exists(self.tmpdir):
            os.makedirs(self.tmpdir)
        self.threads = args.threads
        self.kmer_size = args.kmer_size
        self.subsample_depth = args.subsample_depth
        self.kmer_cutoff = args.kmer_cutoff
        self.databasedir = args.databases

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

    def determine_genus(self):
        """
        Figure out which genus each sample is in using MASH, so genus-specific rMLST databases can be used.
        Hopefully makes sensitivity even higher than before.
        :return:
        """
        # Step 0.1: Link the fastq files in each sample to the tmp directory. If they're already there, perhaps
        # due to a previous failed run that didn't delete the tmp directory, pass.
        fastq_files = glob.glob(self.fastq_directory + '/*.f*q*')
        for name in fastq_files:
            try:
                os.link(name, self.tmpdir + name.split('/')[-1])
            except FileExistsError:
                pass
            except OSError:
                shutil.copy(name, self.tmpdir + name.split('/')[-1])
        # Step 1: Run the mashsippr to determine genus. This should create a mash.csv file in self.tmpdir/reports
        cmd = 'python3 -m confindr.mashsippr -s {} -t {} {}'.format(self.tmpdir, self.databasedir, self.tmpdir)
        with open(self.logfile, 'a+') as logfile:
            subprocess.call(cmd, shell=True, stdout=logfile, stderr=logfile)

        # Now read through the csv, and assign a genus to each sample.
        try:
            with open(self.tmpdir + '/reports/mash.csv') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    for sample in self.samples:
                        if row['Strain'] in self.samples[sample].forward_reads:
                            self.samples[sample].genus = row['ReferenceGenus']
        except FileNotFoundError:
            print('WARNING: Could not find the MASH result file. All samples will be processed with'
                  ' non-specific databases..')

    def prepare_genusspecific_databases(self):
        databases_to_create = list()
        genes_to_exclude = dict()
        # Read in the profiles file.
        f = open('{}/profiles.txt'.format(self.databasedir))
        lines = f.readlines()
        f.close()
        # Parse the profiles file to know what genes to exclude for every genus.
        for line in lines:
            line = line.rstrip()
            data = line.split(':')
            genus = data[0]
            genes = data[1].split(',')
            genes_to_exclude[genus] = genes
        # Figure out which databases you need to create, and put them into a list.
        for sample in self.samples:
            if self.samples[sample].genus not in databases_to_create and self.samples[sample].genus != 'NA'\
                    and self.samples[sample].genus in genes_to_exclude:
                databases_to_create.append(self.samples[sample].genus)
        # Now create the db with the aid of SeqIO.
        for db in databases_to_create:
            genus_database = '{}/{}_db.fasta'.format(self.databasedir, db)
            if not os.path.exists(genus_database):
                print('Creating database for {}...'.format(db))
                f = open(genus_database, 'w')
                sequences = SeqIO.parse(os.path.join(self.databasedir, 'rMLST_combined.fasta'), 'fasta')
                for item in sequences:
                    if item.id.split('_')[0] not in genes_to_exclude[db]:
                        f.write('>' + item.id + '\n')
                        f.write(str(item.seq) + '\n')
                f.close()
        for sample in self.samples:
            genus_database = 'databases/{}_db.fasta'.format(self.samples[sample].genus)
            if os.path.isfile(genus_database):
                self.samples[sample].genus_database = genus_database

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
        bad_samples = list()
        for sample in self.samples:
            # For paired reads.
            if self.samples[sample].reverse_reads != 'NA':
                self.samples[sample].forward_rmlst_reads = self.tmpdir + self.samples[sample].forward_reads.split('/')[-1] + '_rmlst'
                self.samples[sample].reverse_rmlst_reads = self.tmpdir + self.samples[sample].reverse_reads.split('/')[-1] + '_rmlst'
                if self.samples[sample].genus_database == 'NA':
                    db = os.path.join(self.databasedir, 'rMLST_combined.fasta')
                else:
                    db = self.samples[sample].genus_database
                cmd = 'bbduk.sh ref={} in1={} in2={} outm={} outm2={} threads={}'.format(db,
                                                                              self.samples[sample].forward_reads,
                                                                              self.samples[sample].reverse_reads,
                                                                              self.samples[sample].forward_rmlst_reads,
                                                                              self.samples[sample].reverse_rmlst_reads,
                                                                                         str(self.threads))
            else:
                # Unpaired reads.
                db = os.path.join(self.databasedir, 'rMLST_combined.fasta')
                self.samples[sample].forward_rmlst_reads = self.tmpdir + self.samples[sample].forward_reads.split('/')[-1] + '_rmlst'
                cmd = 'bbduk.sh ref={} in={} outm={} threads={}'.format(db,
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
                bad_samples.append(sample)
        # If any of the samples are bad, delete them so they don't go for further analysis.
        for sample in bad_samples:
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
                cmd = 'reformat.sh in1={} in2={} out1={} out2={} overwrite samplebasestarget={}'.format(self.samples[sample].trimmed_forward_reads,
                                                                                                        self.samples[sample].trimmed_reverse_reads,
                                                                                                        self.samples[sample].subsampled_forward,
                                                                                                        self.samples[sample].subsampled_reverse,
                                                                                                        str(self.subsample_depth * 35000))
            else:
                self.samples[sample].subsampled_forward = self.samples[sample].trimmed_forward_reads + '_subsampled'
                cmd = 'reformat.sh in={} out={} overwrite samplebasestarget={}'.format(self.samples[sample].trimmed_forward_reads,
                                                                                       self.samples[sample].subsampled_forward,
                                                                                       str(self.subsample_depth * 35000))
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
                cmd = 'jellyfish count -m {} -s 100M --bf-size 100M -C -F 2 {} {} -o {} -t {}'.format(str(self.kmer_size), self.samples[sample].subsampled_forward,
                                                                                                      self.samples[sample].subsampled_reverse,
                                                                                                      self.samples[sample].jellyfish_file,
                                                                                                      str(self.threads))
            else:
                cmd = 'jellyfish count -m {} -s 100M --bf-size 100M -C -F 1 {} -o {} -t {}'.format(str(self.kmer_size), self.samples[sample].subsampled_forward,
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
                    if int(fastas[i].replace('>', '')) >= self.kmer_cutoff:
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
            if self.samples[sample].genus_database == 'NA':
                db = os.path.join(self.databasedir, 'rMLST_combined.fasta')
            else:
                db = self.samples[sample].genus_database
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
                    if "1X" in match.cigarstring and match.query_alignment_length == self.kmer_size:
                        reference = samfile.getrname(match.reference_id)
                        to_blast.append(mer_dict[reference])
                # BLAST our potentially contaminating sequences in parallel to speed things up because BLASTing
                # is painfully slow.
                pool = multiprocessing.Pool(processes=self.threads)
                db_list = [db] * len(to_blast)
                results = pool.starmap(self.present_in_db, zip(to_blast, db_list))
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
        for sample in self.samples:
            if self.samples[sample].genus_database == 'NA':
                db = os.path.join(self.databasedir, 'rMLST_combined.fasta')
            else:
                db = self.samples[sample].genus_database
            db_files = ['.nhr', '.nin', '.nsq']
            db_present = True
            for db_file in db_files:
                if len(glob.glob(db + '*' + db_file)) == 0:
                    db_present = False
            if not db_present:
                print('Making database!')
                cmd = 'makeblastdb -dbtype nucl -in ' + db
                with open(self.logfile, 'a+') as logfile:
                    subprocess.call(cmd, shell=True, stderr=logfile)

    def present_in_db(self, query_sequence, database):
        """
        Checks if a sequence is present in our rMLST database, as some overhangs on reads can be in repetitive regions
        that could cause false positives and should therefore be screened out.
        :param query_sequence: nucleotide sequence, as a string.
        :param database: Database to be used.
        :return: True if sequence is found in rMLST database, False if it isn't.
        """
        # Blast the sequence against our database.
        blastn = NcbiblastnCommandline(db=database, outfmt=5)
        stdout, stderr = blastn(stdin=query_sequence)
        # If there's any full-length result, the sequence is present. No result means not present.
        if stdout:
            for record in NCBIXML.parse(StringIO(stdout)):
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.align_length == self.kmer_size:
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

    def cross_contamination(self):
        # Need to extract taxonomy information from the query_id field. Seems like it'll be a bit of a pain.
        for sample in self.samples:
            genuses_present = list()
            mash.screen('{}/refseq.msh'.format(self.databasedir), self.samples[sample].forward_reads,
                        self.samples[sample].reverse_reads, threads=self.threads, w='', i='0.95')
            screen_output = mash.read_mash_screen('screen.tab')
            for item in screen_output:
                genus = item.query_id.split('/')[-3]
                if genus not in genuses_present:
                    genuses_present.append(genus)
            self.samples[sample].cross = genuses_present
            if len(self.samples[sample].cross) < 2:
                self.samples[sample].cross = ['NA']

    def write_output(self):
        """
        Writes our output. Takes the median number of contaminating SNVs across our subsample reps, and the highest
        number of unique kmers.
        """
        f = open(self.outfile, 'w')
        f.write('Sample,Genus,NumContamSNVs,NumUniqueKmers,CrossContamination,ContamStatus\n')
        for sample in self.samples:
            try:  # This try/except shouldn't be necessary, but it's good insurance I guess.
                snv_median = statistics.median(self.samples[sample].snv_count)
            except statistics.StatisticsError:
                snv_median = 0
            if snv_median > 0 or self.samples[sample].unique_kmers > 50000 or len(self.samples[sample].cross) > 1:
                contam_status = 'Contaminated'
            else:
                contam_status = 'Clean'
            f.write('{},{},{},{},{},{}\n'.format(sample.split('/')[-1], self.samples[sample].genus, str(statistics.median(self.samples[sample].snv_count)),
                                                str(self.samples[sample].unique_kmers), self.samples[sample].cross,
                                                contam_status))
            # Should also write more detailed statistics somewhere - namely, number of SNVs per subsample in addition
            # to the median.
        f.close()

    def cleanup(self):
        """
        Gets rid of our tmp directory and all its files.
        """
        shutil.rmtree(self.tmpdir)


def check_dependencies():
    dependencies = ['bbmap.sh', 'bbduk.sh', 'blastn', 'jellyfish', 'mash']
    for dep in dependencies:
        is_present = shutil.which(dep)
        if is_present is None:
            raise ModuleNotFoundError('ERROR! Could not find executable for: {}!'.format(dep))


if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_directory", help="Folder that contains fastq files you want to check for contamination. "
                                                "Will find any fastq file that contains .fq or .fastq in the filename.")
    parser.add_argument('output_name', help='Base name for output/temporary directories.')
    parser.add_argument('databases', help='Databases folder. Should contain rMLST_combined.fasta, profiles.txt, '
                                          'and refseq.msh as well as RefSeqSketchesDefaults.msh')
    cpu_count = multiprocessing.cpu_count()
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
    arguments = parser.parse_args()
    check_dependencies()
    detector = Detector(arguments)
    Detector.parse_fastq_directory(detector)
    printtime('Finding genus of samples...', start)
    Detector.determine_genus(detector)
    Detector.prepare_genusspecific_databases(detector)
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
    printtime('Looking for cross-species contamination...', start)
    Detector.cross_contamination(detector)
    printtime('Writing output and cleaning up temporary files!', start)
    Detector.write_output(detector)
    Detector.cleanup(detector)
    printtime('Contamination detection complete! :D', start)
