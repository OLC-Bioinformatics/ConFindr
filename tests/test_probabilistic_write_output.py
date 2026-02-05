#!/usr/bin/env python3

"""
Tests for probabilistic write_output function in confindr_src.methods
"""

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
    assert 'Sample' in header or 'SampleName' in header

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
