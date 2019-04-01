#!/usr/bin/env python

from Bio import SeqIO
import urllib.request
import subprocess
import argparse
import tempfile
import logging
import glob
import csv
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
    parser.add_argument('--desired_number_genes',
                        type=int,
                        default=50,
                        help='Minimum number of genes you want to find.')
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
    potential_genes = get_potential_genes(os.path.join(args.output_folder, 'gene_hit_report.csv'), args.desired_number_genes)
    genomes = sorted(glob.glob(os.path.join(args.output_folder, '*.fasta')))
    confirmed_genes = check_for_similar_genes(potential_genes, genomes)
    for gene in confirmed_genes:
        cmd = 'cat {} >> {}'.format(gene, args.genus + '_db_cgderived.fasta')
        subprocess.call(cmd, shell=True)
    # 5) ???
    # 6) Profit! (but not actually, free and open source, wooooo!)


def check_for_similar_genes(potential_genes, genomes):
    # For each of our potential genes make a blast DB.
    confirmed_genes = list()
    for potential_gene in potential_genes:
        cmd = 'makeblastdb -dbtype nucl -in {}'.format(potential_gene)
        subprocess.call(cmd, shell=True)
    # Then, blast each gene against all other genes, and raise warnings if you find any significant-looking hits.
    for gene1 in potential_genes:
        for gene2 in potential_genes:
            if gene1 != gene2:
                similar_genes_found = False
                with tempfile.TemporaryDirectory() as tmpdir:
                    blast_file = os.path.join(tmpdir, 'blast_out.tsv')
                    cmd = 'blastn -query {seqfile} -db {genome} -out {outfile} -outfmt ' \
                          '"6 qseqid sseqid pident length qlen qstart qend sstart send evalue"'.format(seqfile=gene1,
                                                                                                       genome=gene2,
                                                                                                       outfile=blast_file)
                    subprocess.call(cmd, shell=True)
                    with open(blast_file) as f:
                        for line in f:
                            blast_result = BlastResult(line.rstrip())
                            if blast_result.percent_identity >= 70 or blast_result.query_coverage >= 50:
                                similar_genes_found = True
                if gene1 not in confirmed_genes and similar_genes_found is False:
                    confirmed_genes.append(gene1)
    # Also check that our confirmed genes only hit each genome once, with very loose settings.
    really_confirmed_genes = list()
    for confirmed_gene in confirmed_genes:
        with tempfile.TemporaryDirectory() as tmpdir:
            count = 0
            for contig in SeqIO.parse(confirmed_gene, 'fasta'):
                if count == 0:
                    SeqIO.write([contig], os.path.join(tmpdir, 'sequence.fasta'), 'fasta')
                    count += 1
            only_one_per_genome = True
            for genome in genomes:
                hits = 0
                blast_file = os.path.join(tmpdir, 'blast_out.tsv')
                cmd = 'blastn -query {seqfile} -db {genome} -out {outfile} -outfmt ' \
                      '"6 qseqid sseqid pident length qlen qstart qend sstart send evalue"'.format(seqfile=os.path.join(tmpdir, 'sequence.fasta'),
                                                                                                   genome=genome,
                                                                                                   outfile=blast_file)
                subprocess.call(cmd, shell=True)
                with open(blast_file) as f:
                    for line in f:
                        blast_result = BlastResult(line.rstrip())
                        if blast_result.percent_identity >= 70 or blast_result.query_coverage >= 50:
                            hits += 1
                if hits > 1:
                    only_one_per_genome = False
            if only_one_per_genome is True:
                really_confirmed_genes.append(confirmed_gene)
    return really_confirmed_genes


def get_potential_genes(gene_report, desired_genes):
    proportion_in_genomes = dict()
    potential_genes = list()
    lowest_proportion = 1
    with open(gene_report) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene = row['Gene']
            proportion = float(row['OneHitPerGenome'])
            proportion_in_genomes[gene] = proportion
    sorted_proportions = sorted(proportion_in_genomes.items(), key=lambda kv: kv[1], reverse=True)
    genes_added = 0
    for gene, proportion in sorted_proportions:
        if proportion == 1:
            potential_genes.append(gene)
            genes_added += 1
        elif genes_added < desired_genes:
            potential_genes.append(gene)
            genes_added += 1
            lowest_proportion = proportion
    logging.info('Found {} genes. Lowest proportion found was {}'.format(genes_added, lowest_proportion))
    return potential_genes


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
                    # Make sure files are big enough to be genomes and aren't phage/plasmid/something else.
                    if os.path.getsize(output_file.replace('.gz', '')) < 2000000:
                        os.remove(output_file.replace('.gz', ''))
                        i -= 1
    logging.info('Done downloading! Got {} genomes.'.format(i))


def find_hits_per_genome(genes_folder, genomes_folder):
    # Make blast DBs for all of our genomes.
    genomes = sorted(glob.glob(os.path.join(genomes_folder, '*.fasta')))
    genome_hit_report_file = os.path.join(genomes_folder, 'genome_hit_report.csv')
    gene_report_file = os.path.join(genomes_folder, 'gene_hit_report.csv')
    with open(gene_report_file, 'w') as f:
        f.write('Gene,OneHitPerGenome\n')
    with open(genome_hit_report_file, 'w') as f:
        to_write = 'Gene,'
        for genome in genomes:
            to_write += genome + ','
        to_write = to_write[:-1]
        f.write(to_write + '\n')
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
                    genomes_with_one_hit = 0
                    hits_per_genome = dict()
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
                        hits_per_genome[genome] = number_hits
                        if number_hits == 1:
                            genomes_with_one_hit += 1
                    with open(genome_hit_report_file, 'a+') as f:
                        to_write = gene + ','
                        for genome in genomes:
                            to_write += str(hits_per_genome[genome]) + ','
                        to_write = to_write[:-1]
                        f.write(to_write + '\n')
                    with open(gene_report_file, 'a+') as f:
                        f.write('{},{}\n'.format(gene, genomes_with_one_hit/len(genomes)))
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
