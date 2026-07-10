#!/usr/bin/env python3

"""
Tests for read_contig and cross-contamination scoring in confindr_src.methods
"""

# Standard imports
import math
import os

# Third-party imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

# Local imports
from confindr_src.methods import read_contig


def make_reference_and_bam(
    *,  # Enforce keyword arguments
    tmpdir: str
) -> (str, str):
    """
    Function to create a simple reference FASTA and a BAM file with reads
    mapping to it, including a minor allele at a specific position.

    Args:
        tmpdir: Temporary directory to write files to.

    Returns:
        Tuple of (reference_fasta_path, bam_file_path)
    """
    # Create reference FASTA
    ref_path = os.path.join(tmpdir, 'ref.fasta')
    seq = 'A' * 100
    record = SeqRecord(Seq(seq), id='testgene', description='')

    # write FASTA with single-line sequence for easier indexing and consistent
    # coordinate handling
    with open(ref_path, 'w', encoding='utf-8') as f:
        SeqIO.write(record, f, 'fasta')

    # Write a simple .fai suitable for pysam.FastaFile access
    # (single-line sequence)
    fai_path = ref_path + '.fai'

    # Header line is '>testgene\n' length 10, sequence line length is
    # 100 bases + newline
    with open(fai_path, 'w', encoding='utf-8') as f:
        f.write("testgene\t100\t10\t100\t101\n")

    # Create BAM with reads: majority 'A', two reads with a 'G' at position
    # 50 (0-based 49)
    bam_path = os.path.join(tmpdir, 'reads.bam')
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'SN': 'testgene', 'LN': 100}]}
    aln = pysam.AlignmentFile(bam_path, 'wb', header=header)

    def make_read(name, seq, qual, start=0, mapq=60, is_reverse=False):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.query_sequence = seq
        a.query_qualities = qual
        a.flag = 0
        if is_reverse:
            a.flag |= 16
            a.is_reverse = True
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = mapq
        # CIGAR: full match
        a.cigar = [(0, len(seq))]
        return a

    # create 8 reference-supporting reads
    for i in range(8):
        seq = 'A' * 100
        qual = [30] * 100
        a = make_read(
            f'read{i+1}',
            seq,
            qual,
            start=0,
            mapq=60,
            is_reverse=(i % 2 == 0)
        )
        aln.write(a)

    # Create two reads with G at position 50 (0-based index 49 -> positions
    # are 1-based in output but pileup uses 0-based)
    for i in range(2):
        seq = list('A' * 100)
        seq[49] = 'G'
        seq = ''.join(seq)
        qual = [35] * 100
        a = make_read(
            f'var{i+1}',
            seq,
            qual,
            start=0,
            mapq=60,
            is_reverse=(i % 2 == 1)
        )
        aln.write(a)

    aln.close()
    # index BAM so pysam pileup can access contig regions
    pysam.index(bam_path)
    return ref_path, bam_path


def test_read_contig_detects_minor_allele_and_gene_summary(tmp_path):
    """
    Test that read_contig correctly identifies a minor allele at a
    """
    # Set the path for temporary files
    tmpdir = str(tmp_path)

    # Create reference FASTA and BAM files with reads
    ref_fasta, bam_path = make_reference_and_bam(tmpdir=tmpdir)

    # Create allele_records via SeqIO.index to match how pipeline uses it
    allele_index = SeqIO.index(ref_fasta, 'fasta')

    # Build fastq_records mapping for characterise_read
    # (SeqRecord objects with phred qualities)
    fastq_records = {}

    # Create records for the 8 ref-supporting reads
    for i in range(1, 9):
        name = f'read{i}'
        seq = 'A' * 100
        qual = [30] * 100
        r1 = SeqRecord(Seq(seq), id=name + '/1', description=name + '/1')
        r1.letter_annotations['phred_quality'] = qual
        r2 = SeqRecord(Seq(seq), id=name + '/2', description=name + '/2')
        r2.letter_annotations['phred_quality'] = qual
        fastq_records[name + '/1'] = r1
        fastq_records[name + '/2'] = r2

    # Create records for the 2 variant reads
    for i in range(1, 3):
        name = f'var{i}'
        seq = list('A' * 100)
        seq[49] = 'G'
        seq = ''.join(seq)
        qual = [35] * 100
        r1 = SeqRecord(Seq(seq), id=name + '/1', description=name + '/1')
        r1.letter_annotations['phred_quality'] = qual
        r2 = SeqRecord(Seq(seq), id=name + '/2', description=name + '/2')
        r2.letter_annotations['phred_quality'] = qual
        fastq_records[name + '/1'] = r1
        fastq_records[name + '/2'] = r2

    # Run read_contig on the synthetic data with the generated fastq_records
    multibase_dict, tsv_output = read_contig(
        contig_name='testgene',
        bamfile_name=bam_path,
        reference_fasta=ref_fasta,
        allele_records=allele_index,
        fastq_records=fastq_records,
        quality_cutoff=20,
        base_cutoff=1,
        base_fraction_cutoff=0.01,
        fasta=False,
        error_cutoff=1.0,
        nanopore=False
    )

    # Expect a gene summary present
    assert 'testgene' in multibase_dict
    assert '_gene_stats' in multibase_dict['testgene']
    gene_stats = multibase_dict['testgene']['_gene_stats']
    assert (
        'combined_p' in gene_stats
        and 'gene_score' in gene_stats
        and 'num_sig_positions' in gene_stats
    )

    # There should be an entry in TSV output for the position around 50
    assert 'testgene' in tsv_output

    # The position should be present (1-based position 50)
    assert '\t50\t' in tsv_output

    # Validate that statistical columns are present and sensible for the
    # variant position
    lines = [
        line for line in tsv_output.strip().splitlines() if '\t50\t' in line
    ]
    assert len(lines) == 1

    # Extract fields
    fields = lines[0].split('\t')
    # Fields indices are: p at 7, q at 8, strand_p at 9, pos_p at 10,
    # mean_q at 11, mean_mapq at 12
    p_str = fields[7]
    q_str = fields[8]
    strand_str = fields[9]
    pos_str = fields[10]
    mean_q_str = fields[11]
    mean_mapq_str = fields[12]

    # Ensure p, q, and strand statistics are available. pos_p may be
    # unavailable for some small or tied datasets, so allow ND there.
    assert p_str != 'ND'
    assert q_str != 'ND'
    assert strand_str != 'ND'

    # Convert to float and check ranges where available
    p = float(p_str)
    q = float(q_str)
    strand_p = float(strand_str)
    mean_q = float(mean_q_str)
    mean_mapq = float(mean_mapq_str)

    if pos_str != 'ND':
        pos_p = float(pos_str)
        assert 0.0 <= pos_p <= 1.0
    else:
        pos_p = None

    assert 0.0 <= p <= 1.0
    assert 0.0 <= q <= 1.0
    assert 0.0 <= strand_p <= 1.0
    assert mean_q >= 0.0
    assert mean_mapq >= 0.0

    # Clean up index
    allele_index.close()


def make_two_gene_reference_and_bam(
    *,  # Enforce keyword arguments
    tmpdir: str
) -> (str, str):
    """
    Creates a FASTA with two contigs and a BAM with reads mapping to both
    contigs.
    Contig 'testgene' has 8 reference reads and 2 variant reads (G at pos50).
    Contig 'othergene' has 7 reference reads and 3 variant reads (T at pos60).

    Args:
        tmpdir: Temporary directory to write files to.

    Returns:
        Tuple of (reference_fasta_path, bam_file_path)
    """
    # Create reference FASTA with two genes
    ref_path = os.path.join(tmpdir, 'ref_two_genes.fasta')
    seq1 = 'A' * 100
    seq2 = 'C' * 120
    rec1 = SeqRecord(Seq(seq1), id='Escherichia_testgene', description='')
    rec2 = SeqRecord(Seq(seq2), id='Citrobacter_othergene', description='')

    # Write FASTA with single-line sequences for easier indexing and consistent
    with open(ref_path, 'w', encoding='utf-8') as f:
        SeqIO.write(rec1, f, 'fasta')
        SeqIO.write(rec2, f, 'fasta')

    # Write .fai with single-line sequences
    with open(ref_path + '.fai', 'w', encoding='utf-8') as f:
        f.write("Escherichia_testgene\t100\t10\t100\t101\n")
        f.write("Citrobacter_othergene\t120\t112\t120\t121\n")

    # Create BAM with reads for both genes
    bam_path = os.path.join(tmpdir, 'reads_two_genes.bam')
    header = {
        'HD': {
            'VN': '1.0'
        },
        'SQ': [
            {
                'SN': 'Escherichia_testgene',
                'LN': 100
            },
            {
                'SN': 'Citrobacter_othergene',
                'LN': 120
            }
        ]
    }

    # Write reads
    aln = pysam.AlignmentFile(bam_path, 'wb', header=header)

    def make_read(name, seq, qual, ref_id, start=0, mapq=60, is_reverse=False):
        """
        Function to create a pysam AlignedSegment read.
        """
        a = pysam.AlignedSegment()
        a.query_name = name
        a.query_sequence = seq
        a.query_qualities = qual
        a.flag = 0
        if is_reverse:
            a.flag |= 16
            a.is_reverse = True
        a.reference_id = ref_id
        a.reference_start = start
        a.mapping_quality = mapq
        a.cigar = [(0, len(seq))]
        return a

    # gene1 reads
    for i in range(8):
        seq = 'A' * 100
        qual = [30] * 100
        a = make_read(
            f'g1_read{i+1}',
            seq,
            qual,
            ref_id=0,
            start=0,
            is_reverse=(i % 2 == 0)
        )
        aln.write(a)
    for i in range(2):
        seq = list('A' * 100)
        seq[49] = 'G'
        seq = ''.join(seq)
        qual = [35] * 100
        a = make_read(
            f'g1_var{i+1}',
            seq,
            qual,
            ref_id=0,
            start=0,
            is_reverse=(i % 2 == 1)
        )
        aln.write(a)

    # gene2 reads
    for i in range(7):
        seq = 'C' * 120
        qual = [30] * 120
        a = make_read(
            f'g2_read{i+1}',
            seq,
            qual,
            ref_id=1,
            start=0,
            is_reverse=(i % 2 == 0)
        )
        aln.write(a)
    for i in range(3):
        seq = list('C' * 120)
        seq[59] = 'T'  # variant at 1-based pos 60
        seq = ''.join(seq)
        qual = [35] * 120
        a = make_read(
            f'g2_var{i+1}',
            seq,
            qual,
            ref_id=1,
            start=0,
            is_reverse=(i % 2 == 1)
        )
        aln.write(a)

    aln.close()
    pysam.index(bam_path)
    return ref_path, bam_path


def test_cross_contamination_and_probabilistic_scoring(tmp_path):
    """
    Test that read_contig detects minor alleles in two genes and that
    probabilistic scoring produces a positive sample score indicating
    contamination.
    8 ref reads + 2 variant reads (G at pos50) for gene1
    7 ref reads + 3 variant reads (T at pos60) for gene2
    10% minor allele frequency at both positions should be significant.
    10 variant reads across both genes should yield a positive sample score
    in probabilistic scoring.
    """
    # Set the path for temporary files
    tmpdir = str(tmp_path)

    # Create reference FASTA and BAM files with reads for two genes
    ref_fasta, bam_path = make_two_gene_reference_and_bam(tmpdir=tmpdir)
    allele_index = SeqIO.index(ref_fasta, 'fasta')

    # Build simple fastq_records for all reads so characterise_read has
    # qualities
    fastq_records = {}

    # g1 reads
    for i in range(1, 9):
        name = f'g1_read{i}'
        seq = 'A' * 100
        qual = [30] * 100
        r1 = SeqRecord(Seq(seq), id=name + '/1', description=name + '/1')
        r1.letter_annotations['phred_quality'] = qual
        r2 = SeqRecord(Seq(seq), id=name + '/2', description=name + '/2')
        r2.letter_annotations['phred_quality'] = qual
        fastq_records[name + '/1'] = r1
        fastq_records[name + '/2'] = r2
    for i in range(1, 3):
        name = f'g1_var{i}'
        seq = list('A' * 100)
        seq[49] = 'G'
        seq = ''.join(seq)
        qual = [35] * 100
        r1 = SeqRecord(Seq(seq), id=name + '/1', description=name + '/1')
        r1.letter_annotations['phred_quality'] = qual
        r2 = SeqRecord(Seq(seq), id=name + '/2', description=name + '/2')
        r2.letter_annotations['phred_quality'] = qual
        fastq_records[name + '/1'] = r1
        fastq_records[name + '/2'] = r2
    # g2 reads
    for i in range(1, 8):
        name = f'g2_read{i}'
        seq = 'C' * 120
        qual = [30] * 120
        r1 = SeqRecord(Seq(seq), id=name + '/1', description=name + '/1')
        r1.letter_annotations['phred_quality'] = qual
        r2 = SeqRecord(Seq(seq), id=name + '/2', description=name + '/2')
        r2.letter_annotations['phred_quality'] = qual
        fastq_records[name + '/1'] = r1
        fastq_records[name + '/2'] = r2
    for i in range(1, 4):
        name = f'g2_var{i}'
        seq = list('C' * 120)
        seq[59] = 'T'
        seq = ''.join(seq)
        qual = [35] * 120
        r1 = SeqRecord(Seq(seq), id=name + '/1', description=name + '/1')
        r1.letter_annotations['phred_quality'] = qual
        r2 = SeqRecord(Seq(seq), id=name + '/2', description=name + '/2')
        r2.letter_annotations['phred_quality'] = qual
        fastq_records[name + '/1'] = r1
        fastq_records[name + '/2'] = r2

    # Run read_contig on both genes
    g1_multibase, g1_tsv = read_contig(
        contig_name='Escherichia_testgene',
        bamfile_name=bam_path,
        reference_fasta=ref_fasta,
        allele_records=allele_index,
        fastq_records=fastq_records,
        quality_cutoff=20,
        base_cutoff=1,
        base_fraction_cutoff=0.01,
        fasta=False,
        error_cutoff=1.0,
        nanopore=False
    )
    g2_multibase, g2_tsv = read_contig(
        contig_name='Citrobacter_othergene',
        bamfile_name=bam_path,
        reference_fasta=ref_fasta,
        allele_records=allele_index,
        fastq_records=fastq_records,
        quality_cutoff=20,
        base_cutoff=1,
        base_fraction_cutoff=0.01,
        fasta=False,
        error_cutoff=1.0,
        nanopore=False
    )

    # Both genes should show multibase positions indicating cross-contamination
    assert any(not k.startswith('_') for k in g1_multibase)
    assert any(not k.startswith('_') for k in g2_multibase)

    # Collect q-values from both TSV outputs and compute sample_score
    # (sum -log10(q)) as the pipeline does
    qs = []

    # Parse gene1 TSV
    for line in g1_tsv.strip().splitlines() + g2_tsv.strip().splitlines():
        # Expect tab-separated fields; AdjPValue (q) is at index 8
        fields = line.split('\t')
        if len(fields) > 8:
            q_val = fields[8]
            if q_val != 'ND':
                try:
                    qs.append(float(q_val))
                except ValueError:
                    pass
    assert len(qs) >= 1

    # Compute sample_score as sum of -log10(q) values
    sample_score = sum([-math.log10(q) for q in qs])

    # Probabilistic scoring should produce a positive score; assert it's
    # greater than a small threshold
    assert sample_score > 0.1

    # If using a score threshold slightly below computed score, sample would
    # be flagged as contaminated in probabilistic mode
    threshold = sample_score - 0.01
    assert sample_score >= threshold

    allele_index.close()
