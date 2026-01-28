#!/usr/bin/env python3
"""
create_genus_specific_db.py

Create a genus-specific core-gene derived database using input gene FASTA
files and RefSeq genomes.
"""

# Standard library imports
from glob import glob
from typing import (
    Dict,
    List
)
import argparse
import csv
import logging
import os
import subprocess
import tempfile
import urllib.request

# Third-party imports
from Bio import SeqIO


def main() -> None:
    """
    Command-line entrypoint for creating a genus-specific database.

    This performs the end-to-end workflow:
    1. Download RefSeq assembly summary
    2. Download complete genomes for the requested genus
    3. BLAST query genes against genomes and filter candidates
    4. Remove internally-similar genes and write a combined FASTA

    Returns:
        None
    """
    logging.basicConfig(
        format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o', '--output_folder',
        type=str,
        required=True,
        help=(
            'Folder to first store temporary files, and eventually store the '
            'created database.'
        )
    )
    parser.add_argument(
        '-i', '--input_folder',
        type=str,
        required=True,
        help=(
            'Folder with your input files to try to find core genes. Each '
            'gene should be in a FASTA file. Expected extension is .fasta'
        )
    )
    parser.add_argument(
        '-g', '--genus',
        type=str,
        required=True,
        help='Name of genus you\'re creating a database for.'
    )
    parser.add_argument(
        '--desired_number_genes',
        type=int,
        default=50,
        help='Minimum number of genes you want to find.'
    )
    args = parser.parse_args()

    # Create output folder if it doesn't exist
    os.makedirs(args.output_folder, exist_ok=True)

    # Steps to get this done:
    # 1) Get the RefSeq assembly summary
    # (ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/
    # assembly_summary_refseq.txt)
    download_refseq_summary(
        output_folder=args.output_folder
    )

    # 2) From the RefSeq assembly summary, download complete genomes for your
    # genus of interest.
    download_refseq_genomes(
        output_folder=args.output_folder,
        assembly_summary=os.path.join(
            args.output_folder,
            'assembly_summary_refseq.txt'
        ),
        genus=args.genus
    )

    # 3) BLAST each of the potential genes to be used against the RefSeq
    # genomes of interest. We only want to keep genes that both hit all
    # genomes, and also hit only once per genome.
    find_hits_per_genome(
        genes_folder=args.input_folder,
        genomes_folder=args.output_folder
    )

    # 4) BLAST the potential genes we've found against each other to make sure
    # none of them are similar to each other.
    potential_genes = get_potential_genes(
        gene_report=os.path.join(
            args.output_folder,
            'gene_hit_report.tsv'
        ),
        desired_genes=args.desired_number_genes
    )
    genomes = sorted(glob(os.path.join(args.output_folder, '*.fasta')))

    # Filter out genes that are too similar to each other.
    confirmed_genes = check_for_similar_genes(
        potential_genes=potential_genes,
        genomes=genomes
    )

    # Write out our final database FASTA file.
    for gene in confirmed_genes:
        cmd = f'cat {gene} >> {args.genus}_db_cgderived.fasta'
        subprocess.call(cmd, shell=True)
    # 5) ???
    # 6) Profit! (but not actually, free and open source, wooooo!)


def check_for_similar_genes(
    *,  # Enforce keyword-only arguments
    potential_genes: List[str],
    genomes: List[str]
) -> List[str]:
    """Identify and filter genes that are too similar to other potential genes.

    Args:
        potential_genes: List of file paths to candidate gene FASTA files.
        genomes: List of genome FASTA paths to check hits against.

    Returns:
        A filtered list of confirmed gene file paths.
    """
    # For each of our potential genes make a blast DB.
    confirmed_genes: List[str] = []
    for potential_gene in potential_genes:
        cmd = f'makeblastdb -dbtype nucl -in {potential_gene}'
        subprocess.call(cmd, shell=True)

    # Then, blast each gene against all other genes, and raise warnings if you
    # find any significant-looking hits.
    for gene1 in potential_genes:
        for gene2 in potential_genes:
            if gene1 != gene2:
                similar_genes_found = False
                with tempfile.TemporaryDirectory() as tmpdir:
                    blast_file = os.path.join(tmpdir, 'blast_out.tsv')
                    cmd = (
                        f'blastn -query {gene1} -db {gene2} -out {blast_file} '
                        '-outfmt "6 qseqid sseqid pident length qlen qstart '
                        'qend sstart send evalue"'
                    )
                    subprocess.call(cmd, shell=True)
                    with open(blast_file, encoding='utf-8') as f:
                        for line in f:
                            # Check if the blast hit meets similarity criteria
                            blast_result = BlastResult(line.rstrip())
                            if (
                                blast_result.percent_identity >= 70
                                or blast_result.query_coverage >= 50
                            ):
                                # Set flag
                                similar_genes_found = True

                # If no similar genes were found, add to confirmed genes.
                if (
                    gene1 not in confirmed_genes
                    and similar_genes_found is False
                ):
                    confirmed_genes.append(gene1)

    # Also check that our confirmed genes only hit each genome once, with very
    # loose settings.
    really_confirmed_genes = []
    for confirmed_gene in confirmed_genes:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Initialize count to track first contig only
            count = 0

            # Write out only the first contig from the gene file.
            for contig in SeqIO.parse(confirmed_gene, 'fasta'):
                if count == 0:
                    SeqIO.write(
                        [contig],
                        os.path.join(
                            tmpdir,
                            'sequence.fasta'
                        ),
                        'fasta'
                    )
                    count += 1

            # Blast against each genome and make sure only one hit per genome.
            only_one_per_genome = True
            for genome in genomes:
                # Initialize hit count
                hits = 0

                # Set the path for the blast output file
                blast_file = os.path.join(tmpdir, 'blast_out.tsv')

                # Create a variable for the sequence file path
                seqfile = os.path.join(tmpdir, 'sequence.fasta')

                # Create and run the blast command
                cmd = (
                    f'blastn -query {seqfile} -db {genome} -out {blast_file} '
                    '-outfmt "6 qseqid sseqid pident length qlen qstart qend '
                    'sstart send evalue"'
                )

                # Execute the blast command
                subprocess.call(cmd, shell=True)

                # Parse the blast output and count hits
                with open(blast_file, encoding='utf-8') as f:
                    for line in f:
                        # Create a BlastResult object from the line
                        blast_result = BlastResult(line.rstrip())

                        # Check if the blast hit meets similarity criteria
                        if (
                            blast_result.percent_identity >= 70
                            or blast_result.query_coverage >= 50
                        ):
                            hits += 1

                # If more than one hit found, set flag to False
                if hits > 1:
                    only_one_per_genome = False

            # If only one hit per genome, add to final confirmed genes
            if only_one_per_genome is True:
                really_confirmed_genes.append(confirmed_gene)

    return really_confirmed_genes


def get_potential_genes(
    *,  # Enforce keyword-only arguments
    gene_report: str,
    desired_genes: int
) -> List[str]:
    """
    Select potential genes based on presence across genomes.

    Args:
        gene_report: Path to CSV file with Gene,OneHitPerGenome counts.
        desired_genes: Minimum number of genes to select.

    Returns:
        A list of gene file paths selected as potential core genes.
    """
    # Initialize dictionary to hold proportion of genomes with one hit
    proportion_in_genomes: Dict[str, float] = {}

    # Create list to hold potential genes
    potential_genes: List[str] = []

    # Initialize variable to track lowest proportion
    lowest_proportion: float = 1.0

    # Read the gene report and populate the dictionary
    with open(gene_report, encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Extract gene name and proportion
            gene = row['Gene']
            proportion = float(row['OneHitPerGenome'])
            proportion_in_genomes[gene] = proportion

    # Sort genes by proportion in descending order
    sorted_proportions = sorted(
        proportion_in_genomes.items(),
        key=lambda kv: kv[1],
        reverse=True
    )

    # Initialize count of genes added
    genes_added = 0

    # Select genes based on proportion until desired number is reached
    for gene, proportion in sorted_proportions:
        if proportion == 1:
            potential_genes.append(gene)
            genes_added += 1
        elif genes_added < desired_genes:
            potential_genes.append(gene)
            genes_added += 1
            lowest_proportion = proportion

    # Log the number of genes found and the lowest proportion
    logging.info(
        'Found %s genes. Lowest proportion found was %s',
        genes_added, lowest_proportion
    )

    return potential_genes


def download_refseq_summary(
    *,  # Enforce keyword-only arguments
    output_folder: str
) -> None:
    """
    Download RefSeq assembly_summary_refseq.txt into output_folder.

    Args:
        output_folder: Directory where the assembly summary will be saved.
    """
    logging.info('Downloading RefSeq summary...')

    # Download the assembly summary file
    urllib.request.urlretrieve(
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/'
        'assembly_summary_refseq.txt',
        os.path.join(output_folder, 'assembly_summary_refseq.txt')
    )
    assert os.path.isfile(
        os.path.join(
            output_folder,
            'assembly_summary_refseq.txt'
        )
    )


def download_refseq_genomes(
    *,  # Enforce keyword-only arguments
    output_folder: str,
    assembly_summary: str,
    genus: str
) -> None:
    """
    Download complete RefSeq genomes for a given genus.

    Args:
        output_folder: Directory where downloaded genomes will be placed.
        assembly_summary: Path to the RefSeq assembly_summary_refseq.txt file.
        genus: Genus name to filter organisms by.

    Returns:
        None
    """
    logging.info(
        'Downloading complete RefSeq genomes for %s. Depending on genus, this '
        'may take a while...', genus
    )

    # Initialize genome counter
    i = 0

    # Parse the assembly summary and download genomes
    with open(assembly_summary, encoding='utf-8') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            # Split the line into fields
            x = line.split('\t')
            organism = x[7]
            level = x[11]
            ftp_folder = x[19]
            download_link = (
                ftp_folder + '/' + ftp_folder.split('/')[-1]
                + '_genomic.fna.gz'
            )

            # Check if organism matches genus and is complete
            if (
                genus in organism
                and 'PHAGE' not in organism.upper()
                and 'Complete' in level
            ):
                # Increment genome counter
                i += 1

                # Download and unzip the genome
                output_file = os.path.join(
                    output_folder,
                    f'genome_{i}.fasta.gz'
                )

                # Use urllib to download
                urllib.request.urlretrieve(download_link, output_file)

                # System call to gzip since it's faster
                subprocess.call(f'gunzip {output_file}', shell=True)

                # Make sure files are big enough to be genomes and aren't
                # phage/plasmid/something else.
                if os.path.getsize(output_file.replace('.gz', '')) < 2000000:
                    os.remove(output_file.replace('.gz', ''))
                    # Decrement counter since we didn't keep this genome
                    i -= 1

    logging.info('Done downloading! Got %s genomes.', i)


def find_hits_per_genome(
    *,  # Enforce keyword-only arguments
    genes_folder: str,
    genomes_folder: str
) -> None:
    """
    BLAST each gene against genomes and report single-hit proportions.

    Args:
        genes_folder: Folder containing input gene FASTA files (one allele per
        file).
        genomes_folder: Folder containing RefSeq genomes (downloaded by
        `download_refseq_genomes`).

    Returns:
        None
    """
    # Make blast DBs for all of our genomes.
    genomes = sorted(glob(os.path.join(genomes_folder, '*.fasta')))

    # Set up report file paths
    genome_hit_report_file = os.path.join(
        genomes_folder,
        'genome_hit_report.tsv'
    )
    gene_report_file = os.path.join(genomes_folder, 'gene_hit_report.tsv')

    # Initialize report files
    with open(gene_report_file, 'w', encoding='utf-8') as f:
        f.write('Gene\tOneHitPerGenome\n')

    # Write header for genome hit report
    with open(genome_hit_report_file, 'w', encoding='utf-8') as f:
        to_write = 'Gene\t'
        for genome in genomes:
            to_write += genome + '\t'

        # Remove trailing tab
        to_write = to_write[:-1]
        f.write(to_write + '\n')

    logging.info('Creating BLAST databases for genomes of interest.')
    for genome in genomes:
        cmd = f'makeblastdb -dbtype nucl -in {genome}'

        # Run the command
        subprocess.call(cmd, shell=True)

    # Now that Blast DBs are created, take the first allele from each gene
    # file (it's assumed alleles are REALLY similar), and BLAST it against
    # each of the genomes.
    genes = sorted(glob(os.path.join(genes_folder, '*.fasta')))
    for gene in genes:
        # Initialize allele counter
        i = 0

        for sequence in SeqIO.parse(gene, 'fasta'):
            # Only take the first sequence (allele) from the gene file
            if i == 0:
                with tempfile.TemporaryDirectory() as tmpdir:
                    # Write out the sequence to a temporary FASTA file
                    seqfile = os.path.join(tmpdir, 'sequence.fasta')
                    SeqIO.write([sequence], seqfile, 'fasta')

                    # Initialize count of genomes with one hit and hits
                    # per genome
                    genomes_with_one_hit = 0
                    hits_per_genome = {}

                    # BLAST against each genome
                    for genome in genomes:
                        blast_file = os.path.join(tmpdir, 'blast_out.tsv')

                        # Create and run the blast command
                        cmd = (
                            f'blastn -query {seqfile} -db {genome} '
                            f'-out {blast_file} -outfmt "6 qseqid sseqid '
                            'pident length qlen qstart qend sstart send '
                            'evalue"'
                        )
                        subprocess.call(cmd, shell=True)

                        # Initialize hit counter
                        number_hits = 0

                        # Parse the blast output and count hits
                        with open(blast_file, 'r', encoding='utf-8') as f:
                            for line in f:
                                # Create a BlastResult object from the line
                                blast_result = BlastResult(line.rstrip())

                                # Check if the blast hit meets similarity
                                # criteria
                                if (
                                    blast_result.percent_identity >= 90
                                    and blast_result.query_coverage >= 90
                                ):
                                    number_hits += 1

                        # Store the number of hits for this genome
                        hits_per_genome[genome] = number_hits

                        # If only one hit, increment genomes_with_one_hit
                        if number_hits == 1:
                            genomes_with_one_hit += 1

                    # Write out results to the report files
                    with open(
                        genome_hit_report_file, 'a+', encoding='utf-8'
                    ) as f:
                        # Add a tab to the gene name
                        to_write = gene + '\t'

                        # Iterate through genomes and add hit counts
                        for genome in genomes:
                            to_write += str(hits_per_genome[genome]) + '\t'

                        # Remove trailing tab
                        to_write = to_write[:-1]
                        f.write(to_write + '\n')

                    # Write out to gene report file
                    with open(gene_report_file, 'a+', encoding='utf-8') as f:
                        f.write(
                            f'{gene}\t{genomes_with_one_hit/len(genomes)}\n'
                        )

            # Increment allele counter
            i += 1


class BlastResult:
    """
    Class to hold BLAST result information.
    """
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
        self.query_coverage = (
            100.0 * self.alignment_length/self.query_sequence_length
        )


if __name__ == '__main__':
    main()
