from contam_object import ContamObject
import subprocess
import argparse
import glob
import os


class Detector(object):
    def __init__(self, args):
        self.fastq_directory = args.fastq_directory
        self.samples = dict()
        self.tmpdir = args.output_name + 'tmp/'
        self.database = args.rmlst_database
        if not os.path.exists(self.tmpdir):
            os.makedirs(self.tmpdir)

    def parse_fastq_directory(self):
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
        # Figure out where bbduk is so that we can use the adapter file.
        cmd = 'which bbduk.sh'
        bbduk_dir = subprocess.check_output(cmd.split()).decode('utf-8')
        bbduk_dir = bbduk_dir.split('/')[:-1]
        bbduk_dir = '/'.join(bbduk_dir)
        for sample in self.samples:
            # This captures paired reads.
            if self.samples[sample].reverse_reads != 'NA':
                self.samples[sample].trimmed_forward_reads = self.tmpdir + self.samples[sample].forward_reads.split('/')[-1] + '_trimmed'
                self.samples[sample].trimmed_reverse_reads = self.tmpdir + self.samples[sample].reverse_reads.split('/')[-1] + '_trimmed'
                cmd = 'bbduk.sh in1={} in2={} out1={} out2={} qtrim=w trimq=20 k=25 minlength=50 forcetrimleft=15' \
                  ' ref={}/resources/adapters.fa overwrite hdist=1 tpe tbo'.format(self.samples[sample].forward_reads,
                                                                                   self.samples[sample].reverse_reads,
                                                                                   self.samples[sample].trimmed_forward_reads,
                                                                                   self.samples[sample].trimmed_reverse_reads,
                                                                                   bbduk_dir)
                subprocess.call(cmd, shell=True)
            else:
                self.samples[sample].trimmed_forward_reads = self.tmpdir + self.samples[sample].forward_reads.split('/')[-1] + '_trimmed'
                cmd = 'bbduk.sh in={} out={} qtrim=w trimq=20 k=25 minlength=50 forcetrimleft=15' \
                      ' ref={}/resources/adapters.fa overwrite hdist=1 tpe tbo'.format(self.samples[sample].forward_reads,
                                                                                       self.samples[sample].trimmed_forward_reads,
                                                                                       bbduk_dir)
                subprocess.call(cmd, shell=True)

    def extract_rmlst_reads(self):
        for sample in self.samples:
            if self.samples[sample].trimmed_reverse_reads != 'NA':
                self.samples[sample].forward_rmlst_reads = self.samples[sample].trimmed_forward_reads + '_rmlst'
                self.samples[sample].reverse_rmlst_reads = self.samples[sample].trimmed_reverse_reads + '_rmlst'
                cmd = 'bbduk.sh ref={} in1={} in2={} outm={} outm2={}'.format(self.database,
                                                                              self.samples[sample].trimmed_forward_reads,
                                                                              self.samples[sample].trimmed_reverse_reads,
                                                                              self.samples[sample].forward_rmlst_reads,
                                                                              self.samples[sample].reverse_rmlst_reads)
                subprocess.call(cmd, shell=True)
            else:
                self.samples[sample].forward_rmlst_reads = self.samples[sample].trimmed_forward_reads + '_rmlst'
                cmd = 'bbduk.sh ref={} in={} outm={} '.format(self.database,
                                                              self.samples[sample].trimmed_forward_reads,
                                                              self.samples[sample].forward_rmlst_reads)
                subprocess.call(cmd, shell=True)

    def subsample_reads(self):
        for sample in self.samples:
            if self.samples[sample].reverse_rmlst_reads != 'NA':
                self.samples[sample].subsampled_forward = self.samples[sample].forward_rmlst_reads + '_subsampled'
                self.samples[sample].subsampled_reverse = self.samples[sample].reverse_rmlst_reads + '_subsampled'
                cmd = 'reformat.sh in1={} in2={} out1={} out2={} overwrite samplebasestarget=700000'.format(self.samples[sample].forward_rmlst_reads,
                                                                                                            self.samples[sample].reverse_rmlst_reads,
                                                                                                            self.samples[sample].subsampled_forward,
                                                                                                            self.samples[sample].subsampled_reverse)
                subprocess.call(cmd, shell=True)
            else:
                self.samples[sample].subsampled_forward = self.samples[sample].forward_rmlst_reads + '_subsampled'
                cmd = 'reformat.sh in={} out={} overwrite samplebasestarget=700000'.format(self.samples[sample].forward_rmlst_reads,
                                                                                           self.samples[sample].subsampled_forward)
                subprocess.call(cmd, shell=True)

    def run_jellyfish(self):
        for sample in self.samples:
            self.samples[sample].jellyfish_file = self.tmpdir + self.samples[sample].forward_reads.split('/')[-1] + '_jellyfish'
            # TODO: Add in threading options. Accomodate non-paired reads.
            cmd = 'jellyfish count -m 31 -s 100M --bf-size 100M -t 12 -C -F 2 {} {} -o {}'.format(self.samples[sample].subsampled_forward,
                                                                                                  self.samples[sample].subsampled_reverse,
                                                                                                  self.samples[sample].jellyfish_file)
            subprocess.call(cmd, shell=True)

    def write_mer_file(self):
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
        for sample in self.samples:
            self.samples[sample].samfile = self.samples[sample].mer_fasta.replace('.fasta', '.sam')
            cmd = 'bbmap.sh ref={} in={} outm={} overwrite subfilter=1 insfilter=0delfilter=0 indelfilter=0 ' \
                  'nodisk'.format(self.samples[sample].mer_fasta, self.samples[sample].mer_fasta,
                                  self.samples[sample].samfile)
            subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_directory", help="Folder that contains fastq files you want to check for contamination. "
                                                "Will find any fastq file that contains .fq or .fastq in the filename.")
    parser.add_argument('output_name', help='Base name for output/temporary directories.')
    parser.add_argument('rmlst_database', help='rMLST database, in fasta format.')
    arguments = parser.parse_args()
    detector = Detector(arguments)
    Detector.parse_fastq_directory(detector)
    Detector.quality_trim_reads(detector)
    Detector.extract_rmlst_reads(detector)
    Detector.subsample_reads(detector)
    Detector.run_jellyfish(detector)
    Detector.write_mer_file(detector)
    Detector.run_bbmap(detector)