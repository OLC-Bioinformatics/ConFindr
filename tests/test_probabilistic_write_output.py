#!/usr/bin/env python3

"""
Tests for probabilistic write_output function in confindr_src.methods
"""

# Standard imports
import os

# Local imports
from confindr_src.methods import write_output


def test_write_output_probabilistic_true(tmp_path):
    """
    Test write_output with probabilistic contamination detection set to True
    """
    # Prepare output path
    out = tmp_path / 'confindr_report.tsv'

    # Call write_output with parameters that should yield contamination
    # detected
    write_output(
        output_report=str(out),
        sample_name='TestSample',
        multi_positions=0,
        genus='Escherichia',
        total_gene_length=1000,
        database_download_date='ND',
        sample_score=5.0,
        use_probabilistic=True,
        score_threshold=4.0
    )

    assert out.exists()

    # Read and verify contents
    lines = out.read_text().splitlines()
    assert len(lines) == 2

    # Check header and fields
    header = lines[0]
    assert 'Sample' in header
    assert 'ContamStatus' in header

    # Check contamination status and score
    fields = lines[1].split('\t')
    # ContamStatus is the 4th column (0-based index 3)
    assert fields[3] == 'True'
    # Score should be formatted to 3 decimal places
    assert fields[-1] == '5.000'


def test_write_output_probabilistic_false(tmp_path):
    """
    Test write_output with probabilistic contamination detection set to False
    """
    out = tmp_path / 'confindr_report.tsv'
    write_output(
        output_report=str(out),
        sample_name='TestSample',
        multi_positions=10,
        genus='Escherichia',
        total_gene_length=1000,
        database_download_date='ND',
        sample_score=2.0,
        use_probabilistic=True,
        score_threshold=4.0
    )

    # Read and verify contents
    lines = out.read_text().splitlines()
    fields = lines[1].split('\t')
    assert fields[3] == 'False'
    assert fields[-1] == '2.000'


def test_write_output_subreplicate_counts_header_and_values(tmp_path):
    """
    Test write_output with subreplicate counts and summary columns.
    """
    out = tmp_path / 'confindr_report.tsv'
    subreplicate_counts = [0, 3, 5]

    write_output(
        output_report=str(out),
        sample_name='TestSample',
        multi_positions=0,
        genus='Escherichia',
        total_gene_length=1000,
        database_download_date='ND',
        snp_cutoff=3,
        pysam_pass=True,
        sample_score=3.0,
        use_probabilistic=False,
        score_threshold=4.0,
        subreplicate_counts=subreplicate_counts
    )

    lines = out.read_text().splitlines()
    assert len(lines) == 2
    header = lines[0]
    assert 'MeanContamSNVs' in header
    assert 'Subreplicate_1' in header
    assert 'Subreplicate_3' in header

    fields = lines[1].split('\t')
    assert fields[0] == 'TestSample'
    assert fields[1] == 'Escherichia'
    assert fields[2] == '2.67'
    assert fields[3] == '3.00'
    assert fields[4] == '2.05'
    assert fields[5] == 'False'
    assert fields[6] == '3'
    assert fields[7:] == ['0', '3', '5', '1000', 'ND', '3.000']


def test_write_output_pysam_failure_returns_nd(tmp_path):
    """
    Test write_output handles a failed pysam pass condition.
    """
    out = tmp_path / 'confindr_report.tsv'
    write_output(
        output_report=str(out),
        sample_name='TestSample',
        multi_positions=10,
        genus='Escherichia',
        total_gene_length=1000,
        database_download_date='ND',
        pysam_pass=False,
        sample_score=23.0,
        use_probabilistic=True,
        score_threshold=4.0
    )

    fields = out.read_text().splitlines()[1].split('\t')
    assert fields[3] == 'Pysam SamtoolsError'
    assert fields[2] == 'ND'
    assert fields[-1] == '23.000'


def test_write_output_probabilistic_without_threshold_uses_snp_cutoff(tmp_path):
    """
    Test that probabilistic mode without a threshold falls back to deterministic
    SNV cutoff behavior.
    """
    out = tmp_path / 'confindr_report.tsv'
    write_output(
        output_report=str(out),
        sample_name='TestSample',
        multi_positions=4,
        genus='Escherichia',
        total_gene_length=1000,
        database_download_date='ND',
        snp_cutoff=3,
        sample_score=1.0,
        use_probabilistic=True,
        score_threshold=None
    )

    fields = out.read_text().splitlines()[1].split('\t')
    assert fields[3] == 'True'
    assert fields[2] == '4'
