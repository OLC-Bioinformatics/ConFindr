import jellyfish
import shutil
import os
import pysam
import genome_size


class ContamDetect:

    def parse_fastq_directory(self):
        """
        Should be the first thing called on a ContamDetect object.
        :return: List of fastqpairs in nested array [[forward1, reverse1], [forward2, reverse2]]
        """
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

    def run_jellyfish(self, fastq):
        """
        Runs jellyfish at kmer length of self.kmer_size. Writes kmer sequences to mer_sequences.fasta
        :param fastq: An array with forward reads at index 0 and reverse reads at index 1.
        :return: integer num_mers, which is number of kmers in the reads at that kmer size.
        """
        print('Running jellyfish...')
        cmd = 'jellyfish count -m ' + str(self.kmer_size) + ' -s 100M --bf-size 100M -t 12 -C -F 2 ' + fastq[0] \
               + ' ' + fastq[1]
        os.system(cmd)

        mf = jellyfish.ReadMerFile('mer_counts.jf')
        print('Writing mer sequences to file...')
        outstr = ""
        i = 1
        for mer, count in mf:
            outstr += '>mer' + str(i) + '_' + str(count) + '\n'
            outstr += str(mer) + '\n'
            i += 1
        f = open('mer_sequences.fasta', 'w')
        f.write(outstr)
        f.close()
        num_mers = i
        return num_mers

    @staticmethod
    def run_bbmap(pair):
        """
        Runs bbmap on mer_sequences.fasta, against mer_sequences.fasta, outputting to test.sam. Important to set
        ambig=all so kmers don't just match with themselves.
        """
        if os.path.isdir('ref'):
            shutil.rmtree('ref')
        cmd = 'bbmap.sh ref=mer_sequences.fasta in=mer_sequences.fasta ambig=all outm=tmp/' + pair[0].split('/')[-1] + '.sam subfilter=2'
        os.system(cmd)

    def read_samfile(self, num_mers, fastq):
        i = 1
        print('Reading sam result file...')
        # Open up the alignment file.
        # tot_ratio = 0.0
        samfile = pysam.AlignmentFile('tmp/' + fastq[0].split('/')[-1] + '.sam', 'r')
        for match in samfile:
            # We're interested in full-length matches with one mismatch. This gets us that.
            if "1X" in match.cigarstring and match.query_alignment_length == self.kmer_size:
                query = match.query_name
                reference = samfile.getrname(match.reference_id)
                # If either the query or reference are singletons, chuck them, because those are (probably) just sequencing
                # error noise.
                if '_1' not in query and '_1' not in reference and '_2' not in query and '_2' not in reference:
                    query_kcount = float(query.split('_')[-1])
                    ref_kcount = float(reference.split('_')[-1])
                    # Assuming that the contamination isn't terrible (aka 50/50 or something), we can chuck everything
                    # that has relatively equal ratios of kmers, as we're only interested in low-ratio stuff.
                    if query_kcount/ref_kcount < 0.5 or ref_kcount/query_kcount < 0.5:
                        # if query_kcount < ref_kcount:
                        #     ratio = str(query_kcount/ref_kcount)
                        # else:
                        #     ratio = str(ref_kcount/query_kcount)
                        # tot_ratio += float(ratio)
                        # print(query, reference, ratio)
                        i += 1
        # Try to get estimated genome size.
        print('Estimating genome size...')
        genome_size.run_jellyfish_histo()
        peak, total_mers = genome_size.get_peak_kmers('histogram.txt')
        estimated_size = genome_size.get_genome_size(total_mers, peak)
        f = open(self.output_file, 'a+')
        # Calculate how often we have potentially contaminating kmers.
        percentage = (100.0 * float(i)/float(num_mers))
        f.write(fastq[0] + ',' + str(percentage) + ',' + str(num_mers) + ',' + str(estimated_size) + '\n')
        f.close()
        # print(tot_ratio/float(i))
        # print(i)

    def __init__(self, args):
        self.fastq_folder = args.fastq_folder
        self.output_file = args.output_file
        self.threads = args.threads
        self.kmer_size = args.kmer_size
        if not os.path.isdir('tmp'):
            os.makedirs('tmp')
        f = open(self.output_file, 'w')
        f.write('File,Percentage,NumUniqueKmers,EstimatedGenomeSize\n')
        f.close()


if __name__ == '__main__':
    import argparse
    import time
    import multiprocessing

    # Check the number of CPUs available on the system to be used by bbmap.
    cpu_count = multiprocessing.cpu_count()
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_folder", help="Folder that contains fastq files you want to check for contamination. "
                                             "Will find any fastq file that contains .fq or .fastq in the filename.")
    parser.add_argument("output_file", help="Base name of the output csv you want to create. (.csv extension is added"
                                            "by the program).")
    parser.add_argument("-k", "--kmer_size", type=int, default=31, help="Size of kmer to use. Experimental feature."
                                                                        "Probably don't mess with it.")
    parser.add_argument("-t", "--threads", type=int, default=cpu_count, help="Number of CPUs to run analysis on."
                                                                             " Defaults to number of CPUs on the system.")
    arguments = parser.parse_args()
    detector = ContamDetect(arguments)
    fastq_files = ContamDetect.parse_fastq_directory(detector)
    for pair in fastq_files:
        print(pair)
        # num_mers = 3500000
        num_mers = ContamDetect.run_jellyfish(detector, pair)
        ContamDetect.run_bbmap(pair)
        ContamDetect.read_samfile(detector, num_mers, pair)
    end = time.time()
    # TODO: Delete the tmp folder once you're done with everything.
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)
    print("Finished contamination detection in %d:%02d:%02d " % (h, m, s))

