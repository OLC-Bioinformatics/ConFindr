import jellyfish
import shutil
import os
import pysam
import genome_size
import gzip
import bz2
import subprocess
import run_clark
# from Bio.SeqUtils import GC


# TODO: Add option to try to remove reads that have bad kmers in them - maybe not necessary, but could be useful.
class ContamDetect:

    def parse_fastq_directory(self):
        """
        Should be the first thing called on a ContamDetect object.
        :return: List of fastqpairs in nested array [[forward1, reverse1], [forward2, reverse2]]
        """
        # TODO: Make this handle single reads as well as paired.
        import glob
        # Get a list of all fastq files. For some reason, having/not having the slash doesn't seem to matter on the
        # fastqfolder argument. These should be all the common extensions
        fastq_files = glob.glob(self.fastq_folder + "/*.fastq*")
        fastq_files += glob.glob(self.fastq_folder + "/*.fq*")
        fastq_pairs = list()
        for name in fastq_files:
            # If forward and reverse reads are present, put them in a list of paired files.
            # May need to add support for other naming conventions too. Supports both _R1 and _1 type conventions.
            if "R1" in name and os.path.isfile(name.replace("R1", "R2")):
                fastq_pairs.append([name, name.replace("R1", "R2")])
            # Other naming convention support.
            elif "_1" in name and os.path.isfile(name.replace("_1", "_2")):
                fastq_pairs.append([name, name.replace("_1", "_2")])
        return fastq_pairs

    def run_jellyfish(self, fastq, threads):
        """
        Runs jellyfish at kmer length of self.kmer_size. Writes kmer sequences to mer_sequences.fasta
        :param fastq: An array with forward reads at index 0 and reverse reads at index 1.
        :return: integer num_mers, which is number of kmers in the reads at that kmer size.
        """
        # Send files to check if they're compressed. If they are, create uncompressed version that jellyfish can handle.
        to_remove = list()
        to_use = list()
        for j in range(len(fastq)):
            uncompressed = ContamDetect.uncompress_file(fastq[j])
            if 'bz2' in fastq[j]:
                to_use.append(fastq[j].replace('bz2', ''))
                to_remove.append(fastq[j].replace('.bz2', ''))
            elif 'gz' in fastq[j]:
                to_use.append(fastq[j].replace('.gz', ''))
                to_remove.append(fastq[j].replace('.gz', ''))
            else:
                to_use.append(fastq[j])
        # print('Running jellyfish...')
        # Run jellyfish!
        cmd = 'jellyfish count -m ' + str(self.kmer_size) + ' -s 100M --bf-size 100M -t ' + str(threads) + ' -C -F 2 ' +\
              to_use[0] + ' ' + to_use[1]
        # os.system(cmd)
        subprocess.call(cmd, shell=True)
        os.rename('mer_counts.jf', self.output_file + 'tmp/mer_counts.jf')
        # If we had to uncompress files, remove the uncompressed versions.
        if uncompressed:
            for f in to_remove:
                try:
                    # print(f)
                    os.remove(f)
                except:# Needed in case the file has already been removed - figure out the specific error soon.
                    pass

    # TODO: Consider changing this to just calling jellyfish dump, since then you wouldn't have to get the python
    # jellyfish thing installed (and then you could use python3!).
    def write_mer_file(self, jf_file):
        """
        :param jf_file: .jf file created by jellyfish to be made into a fasta file
        :return: The number of unique kmers in said file.
        """
        mf = jellyfish.ReadMerFile(jf_file)
        # print('Writing mer sequences to file...')
        outstr = list()
        i = 1
        for mer, count in mf:
            if count > 5:
                outstr.append('>mer' + str(i) + '_' + str(count) + '\n')
                outstr.append(str(mer) + '\n')
                i += 1
        f = open(self.output_file + 'tmp/mer_sequences.fasta', 'w')
        f.write(''.join(outstr))
        f.close()
        num_mers = i
        return num_mers

    def run_bbmap(self, pair, threads):
        """
        Runs bbmap on mer_sequences.fasta, against mer_sequences.fasta, outputting to test.sam. Important to set
        ambig=all so kmers don't just match with themselves.
        """
        # TODO: Try to figure out how to make this samfile smaller by excluding useless info so parsing is faster.
        if os.path.isdir('ref'):
            shutil.rmtree('ref')
        cmd = 'bbmap.sh ref=' + self.output_file + 'tmp/mer_sequences.fasta in=' + self.output_file + 'tmp/mer_sequences.fasta ambig=all ' \
              'outm=' + self.output_file + 'tmp/' + pair[0].split('/')[-1] + '.sam idfilter=0.96 minid=0.95 subfilter=1 insfilter=0 ' \
                                                     'delfilter=0 indelfilter=0 nodisk threads=' + str(threads)
        # os.system(cmd)
        # print('Running bbmap...')
        with open(self.output_file + 'tmp/junk.txt', 'w') as outjunk:
            subprocess.call(cmd, shell=True, stderr=outjunk)

    @staticmethod
    def uncompress_file(filename):
        """
        If a file is gzipped or bzipped, creates an uncompressed copy of that file in the same folder
        :param filename: Path to file you want to uncompress
        :return: True if the file needed to be uncompressed, otherwise false.
        """
        uncompressed = False
        if ".gz" in filename:
            in_gz = gzip.open(filename, 'rb')
            out = open(filename.replace('.gz', ''), 'w')
            out.write(in_gz.read())
            out.close()
            uncompressed = True
        elif ".bz2" in filename:
            in_bz2 = bz2.BZ2File(filename, 'rb')
            out = open(filename.replace('.bz2', ''), 'w')
            out.write(in_bz2.read())
            out.close()
            uncompressed = True
        return uncompressed

    def read_samfile(self, num_mers, fastq):
        """
        :param num_mers: Number of unique kmers for the sample be looked at. Found by write_mer_file.
        :param fastq: Array with forward read filepath at index 0, reverse read filepath at index 1.
        Parse through the SAM file generated by bbmap to find how often contaminating alleles are present.
        Also calls methods from genome_size.py in order to estimate genome size (good for finding cross-species contam).
        Writes results to user-specified output file.
        """
        i = 1
        # print('Reading sam result file...')
        # Open up the alignment file.
        # tot_ratio = 0.0
        # os.system('grep -v "31=" ' + self.output_file + 'tmp/' + fastq[0].split('/')[-1] + '.sam > test.sam')
        samfile = pysam.AlignmentFile(self.output_file + 'tmp/' + fastq[0].split('/')[-1] + '.sam', 'r')
        # samfile = pysam.AlignmentFile('test.sam', 'r')
        for match in samfile:
            # We're interested in full-length matches with one mismatch. This gets us that.
            if "1X" in match.cigarstring and match.query_alignment_length == self.kmer_size:
                query = match.query_name
                reference = samfile.getrname(match.reference_id)
                # If either the query or reference are singletons, chuck them, because those are (probably) just sequencing
                # error noise.
                # TODO: Try adjusting this cutoff - looks to be (potentially) too low.
                if '_1' not in query and '_1' not in reference and '_2' not in query and '_2' not in reference:
                    query_kcount = float(query.split('_')[-1])
                    ref_kcount = float(reference.split('_')[-1])
                    # Assuming that the contamination isn't terrible (aka 50/50 or something), we can chuck everything
                    # that has relatively equal ratios of kmers, as we're only interested in low-ratio stuff.
                    if query_kcount/ref_kcount < 0.5 or ref_kcount/query_kcount < 0.5:
                        i += 1
        # Try to get estimated genome size.
        # print('Estimating genome size...')
        # Make jellyfish run a histogram.
        genome_size.run_jellyfish_histo(self.output_file)
        # Find total number of mers and peak coverage value so estimated genome size can be calculated.
        peak, total_mers = genome_size.get_peak_kmers(self.output_file + 'tmp/histogram.txt')
        # Calculate the estimated size
        estimated_size = genome_size.get_genome_size(total_mers, peak)
        # Large estimated size means cross species contamination is likely. Run CLARK-light to figure out which species
        # are likely present
        if estimated_size > 10000000 and self.classify:
            printtime('Cross-species contamination suspected! Running CLARK for classification.', self.start)
            run_clark.classify_metagenome('bacteria/', fastq, self.threads)
            clark_results = run_clark.read_clark_output('abundance.csv')
        else:
            clark_results = 'NA'
        # print('Estimating coverage...')
        estimated_coverage = ContamDetect.estimate_coverage(estimated_size, fastq)
        f = open(self.output_file, 'a+')
        # Calculate how often we have potentially contaminating kmers.
        percentage = (100.0 * float(i)/float(num_mers))
        f.write(fastq[0] + ',' + str(percentage) + ',' + str(num_mers) + ',' + str(estimated_size) + ',' +
                str(estimated_coverage) + ',' + clark_results + '\n')
        f.close()
        # print(tot_ratio/float(i))
        # print(i)

    @staticmethod
    def estimate_coverage(estimated_size, pair):
        """
        :param estimated_size: Estimated size of genome, in basepairs. Found using genome_size.get_genome_size
        :param pair: Array with structure [path_to_forward_reads, path_to_reverse_reads].
        :return: Estimated coverage depth of genome, as an integer.
        """
        # Use some shell magic to find how many basepairs in forward fastq file - cat it into paste, which lets cut take
        # only the second column (which has the sequence), and then count the number of characters.
        if ".gz" in pair[0]:
            cmd = 'zcat ' + pair[0] + ' | paste - - - - | cut -f 2 | wc -c'
        elif ".bz2" in pair[0]:
            cmd = 'bzcat ' + pair[0] + ' | paste - - - - | cut -f 2 | wc -c'
        else:
            cmd = 'cat ' + pair[0] + ' | paste - - - - | cut -f 2 | wc -c'
        number_bp = int(subprocess.check_output(cmd, shell=True))
        # Multiply by two, since reverse reads are present.
        number_bp *= 2
        return number_bp/estimated_size

    def __init__(self, args, start):
        self.fastq_folder = args.fastq_folder
        self.output_file = args.output_file
        self.threads = args.threads
        self.kmer_size = args.kmer_size
        self.classify = args.classify
        self.start = start
        if not os.path.isdir(self.output_file + 'tmp'):
            os.makedirs(self.output_file + 'tmp')
        f = open(self.output_file, 'w')
        f.write('File,Percentage,NumUniqueKmers,EstimatedGenomeSize,EstimatedCoverage,CrossContamination\n')
        f.close()


if __name__ == '__main__':
    import argparse
    import time
    import multiprocessing
    from accessoryFunctions.accessoryFunctions import printtime

    # Check the number of CPUs available on the system to be used by bbmap.
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
    arguments = parser.parse_args()
    detector = ContamDetect(arguments, start)
    paired_files = ContamDetect.parse_fastq_directory(detector)
    sample_num = 1
    for pair in paired_files:
        for i in range(len(pair)):
            pair[i] = os.path.abspath(pair[i])
        printtime('Working on sample ' + str(sample_num) + ' of ' + str(len(paired_files)), start)
        # num_mers = 3500000
        printtime('Running jellyfish...', start)
        ContamDetect.run_jellyfish(detector, pair, arguments.threads)
        printtime('Writing mers to file...', start)
        num_mers = ContamDetect.write_mer_file(detector, arguments.output_file + 'tmp/mer_counts.jf')
        printtime('Finding mismatching mers...', start)
        ContamDetect.run_bbmap(detector, pair, arguments.threads)
        # TODO: Get the samfile reading done in parallel maybe?
        printtime('Generating contamination statistics...', start)
        ContamDetect.read_samfile(detector, num_mers, pair)
        sample_num += 1
    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)
    shutil.rmtree(arguments.output_file + 'tmp')
    printtime("Finished contamination detection!", start)

