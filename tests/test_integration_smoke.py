#!/usr/bin/env python3

"""
Tests for integration of ConFindr components, including output writing and
CLI argument parsing for probabilistic scoring.
"""

# Standard imports
import subprocess

# Local imports
from confindr_src.methods import write_output


def test_write_output_creates_tsv(tmp_path):
    """
    Test that write_output creates a TSV report with expected headers.
    """
    # Define output path
    out = tmp_path / 'confindr_report.tsv'

    # Call write_output with test parameters
    write_output(
        output_report=str(out),
        sample_name='TestSample',
        multi_positions=2,
        genus='Fakella',
        total_gene_length=1000,
        database_download_date='ND'
    )
    assert out.exists()

    # Read first line of output TSV to check headers
    first_line = out.read_text().splitlines()[0]

    # Ensure TSV header likely contains key fields
    assert 'Sample' in first_line or 'SampleName' in first_line
    assert 'ContamStatus' in first_line or 'Contamination' in first_line


def test_confindr_help_includes_prob_scoring_and_downsampling():
    """
    Test that the ConFindr CLI help includes probabilistic scoring and
    downsampling/subreplicate options.
    """
    result = subprocess.run(
        [
            'python3',
            'confindr_src/confindr.py',
            '-h'
        ],
        capture_output=True,
        text=True,
        check=True
    )

    assert result.returncode == 0
    assert '--use-prob-scoring' in result.stdout
    assert '--score-threshold' in result.stdout
    assert '--downsample_depth' in result.stdout
    assert '--subreplicates' in result.stdout
    assert '--subreplicate-seed' in result.stdout
    assert '--subreplicate-consensus' in result.stdout
