#!/usr/bin/env python

import multiprocessing
from Bio import SeqIO
import urllib.request
import subprocess
import argparse
import tempfile
import logging
import glob
import os


def main():
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_folder',
                        type=str,
                        required=True,
                        help='Folder to first store temporary files, and eventually store the created database.')
    parser.add_argument('-i', '--input_folder',
                        type=str,
                        required=True,
                        help='Folder with your input files to try to find core genes. Each gene should be in a '
                             'FASTA file. Expected extension is .fasta')
    parser.add_argument('-g', '--genus',
                        type=str,
                        required=True,
                        help='Name of genus you\'re creating a database for.')
    args = parser.parse_args()

    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)
    # Steps to get this done:
    # 1) Get the RefSeq assembly summary (ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt)
    download_refseq_summary(args.output_folder)
    # 2) From the RefSeq assembly summary, download complete genomes for your genus of interest.
    download_refseq_genomes(args.output_folder, os.path.join(args.output_folder, 'assembly_summary_refseq.txt'), args.genus)
    # 3) BLAST each of the potential genes to be used against the RefSeq genomes of interest. We only want to keep genes
    # that both hit all genomes, and also hit only once per genome.
    find_hits_per_genome(args.input_folder, args.output_folder)
    # 4) BLAST the potential genes we've found against each other to make sure none of them are similar to each other.
    # 5) ???
    # 6) Profit! (but not actually, free and open source, wooooo!)


def download_refseq_summary(output_folder):
    logging.info('Downloading RefSeq summary...')
    urllib.request.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt',
                               os.path.join(output_folder, 'assembly_summary_refseq.txt'))
    assert os.path.isfile(os.path.join(output_folder, 'assembly_summary_refseq.txt'))


def download_refseq_genomes(output_folder, assembly_summary, genus):
    logging.info('Downloading complete RefSeq genomes for {}. Depending on genus, this may take a while...'.format(genus))
    i = 0
    with open(assembly_summary) as f:
        for line in f:
            if not line.startswith('#'):
                x = line.split('\t')
                organism = x[7]
                level = x[11]
                ftp_folder = x[19]
                download_link = ftp_folder + '/' + ftp_folder.split('/')[-1] + '_genomic.fna.gz'
                if genus in organism and 'PHAGE' not in organism.upper() and 'Complete' in level:
                    i += 1
                    output_file = os.path.join(output_folder, 'genome_{}.fasta.gz'.format(i))
                    urllib.request.urlretrieve(download_link, output_file)
                    # System call to gzip since it's faster
                    subprocess.call('gunzip {}'.format(output_file), shell=True)
    logging.info('Done downloading! Got {} genomes.'.format(i))


def find_hits_per_genome(genes_folder, genomes_folder):
    # Make blast DBs for all of our genomes.
    genomes = sorted(glob.glob(os.path.join(genomes_folder, '*.fasta')))
    logging.info('Creating BLAST databases for genomes of interest.')
    for genome in genomes:
        cmd = 'makeblastdb -dbtype nucl -in {}'.format(genome)
        subprocess.call(cmd, shell=True)
    # Now that Blast DBs are created, take the first allele from each gene file (it's assumed alleles are REALLY
    # similar), and BLAST it against each of the genomes.
    genes = sorted(glob.glob(os.path.join(genes_folder, '*.fasta')))
    for gene in genes:
        i = 0
        for sequence in SeqIO.parse(gene, 'fasta'):
            if i == 0:
                with tempfile.TemporaryDirectory() as tmpdir:
                    seqfile = os.path.join(tmpdir, 'sequence.fasta')
                    SeqIO.write([sequence], seqfile, 'fasta')
                    one_hit_in_all_genomes = True
                    for genome in genomes:
                        blast_file = os.path.join(tmpdir, 'blast_out.tsv')
                        cmd = 'blastn -query {seqfile} -db {genome} -out {outfile} -outfmt ' \
                              '"6 qseqid sseqid pident length qlen qstart qend sstart send evalue"'.format(seqfile=seqfile,
                                                                                                           genome=genome,
                                                                                                           outfile=blast_file)
                        subprocess.call(cmd, shell=True)
                        number_hits = 0
                        with open(blast_file) as f:
                            for line in f:
                                blast_result = BlastResult(line.rstrip())
                                if blast_result.percent_identity >= 90 and blast_result.query_coverage >= 90:
                                    number_hits += 1
                        if number_hits != 1:
                            print(genome)
                            one_hit_in_all_genomes = False
                    print('{},{}'.format(gene, one_hit_in_all_genomes))
            i += 1


class BlastResult:
    def __init__(self, blast_tabdelimited_line):
        # With my custom output format, headers are:
        # Index 0: query sequence name
        # Index 1: subject sequence name
        # Index 2: percent identity
        # Index 3: alignment length
        # Index 4: query sequence length
        # Index 5: query start position
        # Index 6: query end position
        # Index 7: subject start position
        # Index 8: subject end position
        # Index 9: evalue
        x = blast_tabdelimited_line.rstrip().split()
        self.query_name = x[0]
        self.subject_name = x[1]
        self.percent_identity = float(x[2])
        self.alignment_length = int(x[3])
        self.query_sequence_length = int(x[4])
        self.query_start_position = int(x[5])
        self.query_end_position = int(x[6])
        self.subject_start_position = int(x[7])
        self.subject_end_position = int(x[8])
        self.evalue = float(x[9])
        # Also need to have amount of query sequence covered as a percentage.
        self.query_coverage = 100.0 * self.alignment_length/self.query_sequence_length


if __name__ == '__main__':
    main()
