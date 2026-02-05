#!/usr/bin/env python3

"""
Test suite for ConFindr methods.
"""

# Standard imports
import csv
import os
import subprocess
import shutil

# Third-party imports
from Bio import SeqIO
import pytest

# Local imports
from confindr_src.methods import (
    dependency_check,
    find_paired_reads,
    find_unpaired_reads,
    run_cmd,
    number_of_bases_above_threshold,
    check_valid_base_fraction,
    find_total_sequence_length,
    write_output,
    base_dict_to_string,
    check_acceptable_xmx,
    load_fastq_records
)

# Ensure that the parent directory is in the system path for imports
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)


def test_integration():
    """
    Integration test for ConFindr using test samples.
    This test runs ConFindr on a set of test samples and verifies the
    contamination status and genus calls against expected results.
    7 test samples are used, with known contamination statuses and genera.
    3 samples are contaminated, 4 are not.
    1 sample is cross-contaminated with multiple genera.
    2 samples are Escherichia, 2 are Salmonella, 2 are Listeria, and
    1 is Salmonella:Citrobacter.
    1 sample is cross-contaminated with Escherichia, Salmonella, and Listeria.
    """
    # Expected results
    correct_contamination_calls = {
        'SRX5084910_SRR8268082': 'True',
        'SRX5084911_SRR8268081': 'False',
        'SRX5084914_SRR8268078': 'True',
        'SRX5084915_SRR8268077': 'False',
        'SRX5084941_SRR8268051': 'True',
        'SRX5084940_SRR8268052': 'False',
        'SRX5084995_SRR8267997': 'True'
    }
    correct_genera = {
        'SRX5084910_SRR8268082': 'Escherichia',
        'SRX5084911_SRR8268081': 'Escherichia',
        'SRX5084914_SRR8268078': 'Salmonella',
        'SRX5084915_SRR8268077': 'Salmonella',
        'SRX5084940_SRR8268052': 'Listeria',
        'SRX5084941_SRR8268051': 'Listeria',
        'SRX5084995_SRR8267997': 'Salmonella:Citrobacter'
    }

    # Run ConFindr
    cmd = [
        "confindr.py",
        "-i", "tests/test_samples",
        "-o", "confindr_integration_output",
        "-d", "databases",
        "-k"
    ]

    # Execute the command
    subprocess.call(" ".join(cmd), shell=True)

    # Define output path
    out_path = 'confindr_integration_output/confindr_report.tsv'

    # Ensure output file exists
    if not os.path.exists(out_path):
        # Integration runs are network- and tool-dependent; skip when output
        # is not produced in this environment
        pytest.skip(
            'Integration run did not produce output (requires '
            'databases/downloads)'
        )

    # Verify results
    with open(out_path, encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            # Get sample name
            sample = row['Sample']

            # Check contamination status and genus
            if 'cross_contaminated' not in sample:
                assert row['ContamStatus'] == correct_contamination_calls[
                    sample
                ]
                assert row['Genus'] == correct_genera[sample]
            else:
                assert row['ContamStatus'] == correct_contamination_calls[
                    sample
                ]
                genera = row['Genus'].split(':')
                # Check that all three genera are reported in the cross-
                # contaminated sample
                assert (
                    'Salmonella' in genera
                    and 'Escherichia' in genera
                    and 'Listeria' in genera
                )
    # Clean up output directories
    shutil.rmtree('confindr_integration_output')
    shutil.rmtree('databases')


def test_present_dependency():
    """
    Test that a known present dependency is detected.
    """
    assert dependency_check(dependency='ls') is True


def test_nonexistent_dependency():
    """
    Test that a known absent dependency is not detected.
    """
    assert dependency_check(dependency='fake_dependency') is False


def test_r1_fastqs():
    """
    Test finding R1/R2 paired FASTQs.
    """
    assert find_paired_reads(
        fastq_directory='tests/fake_fastqs/'
    ) == [
        [
            'tests/fake_fastqs/test_R1.fastq.gz',
            'tests/fake_fastqs/test_R2.fastq.gz'
        ]
    ]


def test_1_fastqs():
    """
    Test finding _1/_2 paired FASTQs.
    """
    assert find_paired_reads(
        fastq_directory='tests/fake_fastqs/',
        forward_id='_1',
        reverse_id='_2'
    ) == [
        [
            'tests/fake_fastqs/test_1.fastq.gz',
            'tests/fake_fastqs/test_2.fastq.gz'
        ]
    ]


def test_empty_fastqs():
    """
    Test finding no paired FASTQs with non-matching identifiers.
    """
    assert not find_paired_reads(
        fastq_directory='tests/fake_fastqs/',
        forward_id='_asdf',
        reverse_id='_fdsa'
    )


def test_unpaired_fastq():
    """
    Test finding unpaired FASTQs.
    """
    assert ['tests/fake_fastqs/test_alone.fastq.gz'] == find_unpaired_reads(
        fastq_directory='tests/fake_fastqs'
    )[2]


def test_run_cmd_success():
    """
    Test running a successful command.
    """
    cmd = 'echo asdf'
    out, err = run_cmd(cmd=cmd)
    assert out == 'asdf\n'
    assert err == ''


def test_run_cmd_failure_exit_code():
    """
    Test running a failing command with non-zero exit code.
    """
    with pytest.raises(subprocess.CalledProcessError):
        run_cmd(cmd='garbagecommandthatdoesnotwork')


def test_two_hq_bases_above_threshold():
    """
    Test with two high-quality bases above threshold.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 80, 'A': 20}
    ) == 2


def test_just_one_hq_bases_above_threshold():
    """
    Test with just one high-quality base above threshold.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 99, 'A': 1}
    ) == 1


def test_two_hq_bases_above_threshold_custom_params():
    """
    Test with two high-quality bases above threshold with custom parameters.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 99, 'A': 1},
        base_count_cutoff=1
    ) == 2


def test_just_one_hq_base_above_threshold_custom_params():
    """
    Test with just one high-quality base above threshold with custom parameters
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 96, 'A': 4},
        base_count_cutoff=5
    ) == 1


def test_three_hq_bases_above_threshold():
    """
    Test with three high-quality bases above threshold.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 90, 'A': 10, 'T': 10}
    ) == 3


def test_two_out_of_three_hq_bases_above_threshold():
    """
    Test with two out of three high-quality bases above threshold.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 90, 'A': 9, 'T': 1}
    ) == 2


def test_two_hq_bases_above_fraction_threshold():
    """
    Test with two high-quality bases above fraction threshold.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 80, 'A': 20},
        base_fraction_cutoff=0.05
    ) == 2


def test_two_hq_bases_above_fraction_threshold_low_coverage():
    """
    Test with two high-quality bases above fraction threshold with low coverage
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 9, 'A': 1},
        base_fraction_cutoff=0.05
    ) == 1


def test_two_hq_bases_above_fraction_threshold_low_coverage_one_base_counts():
    """
    Test with two high-quality bases above fraction threshold with low coverage
    and custom base count cutoff.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 9, 'A': 1},
        base_count_cutoff=1,
        base_fraction_cutoff=0.05
    ) == 2


def test_just_one_hq_bases_above_fraction_threshold():
    """
    Test with just one high-quality base above fraction threshold.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 99, 'A': 1},
        base_fraction_cutoff=0.05
    ) == 1


def test_three_hq_bases_above_fraction_threshold():
    """
    Test with three high-quality bases above fraction threshold.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 90, 'A': 10, 'T': 10},
        base_fraction_cutoff=0.05
    ) == 3


def test_two_out_of_three_hq_bases_above_fraction_threshold():
    """
    Test with two out of three high-quality bases above fraction threshold.
    """
    assert number_of_bases_above_threshold(
        high_quality_base_count={'G': 90, 'A': 9, 'T': 1},
        base_fraction_cutoff=0.05
    ) == 2


def test_valid_base_fraction_none():
    """
    Test valid base fraction with None input.
    """
    assert check_valid_base_fraction(base_fraction=None) is True


def test_valid_base_fraction_zero():
    """
    Test valid base fraction with zero input.
    """
    assert check_valid_base_fraction(base_fraction=0.0) is True


def test_valid_base_fraction_one():
    """
    Test valid base fraction with one input.
    """
    assert check_valid_base_fraction(base_fraction=1.0) is True


def test_valid_base_fraction_between_zero_one():
    """
    Test valid base fraction with a value between zero and one.
    """
    assert check_valid_base_fraction(base_fraction=0.2) is True


def test_invalid_base_fraction():
    """
    Test invalid base fraction with a value greater than one.
    """
    assert check_valid_base_fraction(base_fraction=1.2) is False


def test_total_length_fasta():
    """
    Test total sequence length calculation from a FASTA file.
    """
    assert find_total_sequence_length(fasta_file='tests/rmlst.fasta') == 20862


def test_write_output_creates_file_if_does_not_exist():
    """
    Test that write_output creates a new file if it does not exist.
    """
    write_output(
        output_report='tests/confindr_report.tsv',
        sample_name='Test',
        multi_positions=55,
        genus='Fakella',
        total_gene_length=888,
        database_download_date='ND'
    )
    assert os.path.isfile('tests/confindr_report.tsv') is True


def test_write_output_appends_if_file_does_exist():
    """
    Test that write_output appends to an existing file.
    """
    write_output(
        output_report='tests/confindr_report.tsv',
        sample_name='Test',
        multi_positions=55,
        genus='Fakella',
        total_gene_length=888,
        database_download_date='ND'
    )
    with open('tests/confindr_report.tsv', encoding='utf-8') as f:
        lines = f.readlines()
    assert len(lines) > 2


def test_base_dict_to_string_two_base_descending():
    """
    Test base_dict_to_string with two bases in descending order.
    """
    assert base_dict_to_string(base_dict={'A': 18, 'C': 3}) == 'A:18;C:3'


def test_base_dict_to_string_two_base_ascending():
    """
    Test base_dict_to_string with two bases in ascending order.
    """
    assert base_dict_to_string(base_dict={'A': 8, 'C': 33}) == 'C:33;A:8'


def test_base_dict_to_string_three_bases():
    """
    Test base_dict_to_string with three bases.
    """
    assert base_dict_to_string(
        base_dict={'A': 5, 'T': 88, 'C': 33}
    ) == 'T:88;C:33;A:5'


def test_valid_xmx_string_gigabytes():
    """
    Test valid Xmx string in gigabytes.
    """
    assert check_acceptable_xmx(xmx_string='20g') is True
    assert check_acceptable_xmx(xmx_string='20G') is True


def test_valid_xmx_string_megabytes():
    """
    Test valid Xmx string in megabytes.
    """
    assert check_acceptable_xmx(xmx_string='20m') is True
    assert check_acceptable_xmx(xmx_string='20M') is True


def test_valid_xmx_string_kilobytes():
    """
    Test valid Xmx string in kilobytes.
    """
    assert check_acceptable_xmx(xmx_string='550k') is True
    assert check_acceptable_xmx(xmx_string='550K') is True


def test_invalid_xmx_bad_suffix():
    """
    Test invalid Xmx string with bad suffix.
    """
    assert check_acceptable_xmx(xmx_string='600u') is False


def test_invalid_xmx_float():
    """
    Test invalid Xmx string with float value.
    """
    assert check_acceptable_xmx(xmx_string='2.2G') is False


def test_invalid_xmx_not_an_integer():
    """
    Test invalid Xmx string that is not an integer.
    """
    assert check_acceptable_xmx(xmx_string='asdfK') is False


# FASTQ headers can be in different formats, e.g. Casava 1.8, deposited in SRA,
# pre-Casava, etc., and may also be split across multiple lines.
# A Pytest fixture is defined first, and then the unit tests for each different
# FASTQ header format afterwards.
@pytest.fixture(name='parse_fastq_header', scope='function')
def _parse_fastq_header(
    request: pytest.FixtureRequest
) -> tuple[list[str], list[str], list[str], list[str]]:
    """
    Fixture to compare the FASTQ headers obtained from two different methods:
    `load_fastq_records` and `SeqIO.parse`. It ensures that the headers are
    parsed consistently by both methods.

    Args:
        request: A pytest request object that provides access to the parameters
        passed to the fixture.

    Returns:
        A tuple containing four lists:
        - r1_load_fastq_names: Sorted list of FASTQ header names from
        `load_fastq_records` for forward reads.
        - r2_load_fastq_names: Sorted list of FASTQ header names from
        `load_fastq_records` for reverse reads.
        - r1_characterise_read_names: Sorted list of FASTQ header names from
        `SeqIO.parse` for forward reads.
        - r2_characterise_read_names: Sorted list of FASTQ header names from
        `SeqIO.parse` for reverse reads.
    """
    # Get the R1 and R2 FASTQ file paths from the test parameters
    r1_path, r2_path = request.param

    # Obtain the parts of the FASTQ headers from load_fastq_records()
    r1_load_fastq = load_fastq_records(
        gz=r1_path,
        paired=True,
        forward=True
    )
    r2_load_fastq = load_fastq_records(
        gz=r2_path,
        paired=True,
        forward=False
    )

    # Sort the header names obtained from load_fastq_records() for comparison
    r1_load_fastq_names = sorted(list(r1_load_fastq.keys()))
    r2_load_fastq_names = sorted(list(r2_load_fastq.keys()))

    # Obtain the parts of the FASTQ headers from characterise_read()
    r1_characterise_read = SeqIO.to_dict(
        SeqIO.parse(
            r1_path, 'fastq'
        )
    )

    r2_characterise_read = SeqIO.to_dict(
        SeqIO.parse(
            r2_path, 'fastq'
        )
    )

    # Initialise lists to store the header names obtained from
    # characterise_read()
    r1_characterise_read_names = []
    r2_characterise_read_names = []

    # Similar logic as in characterise_read() to extract the read names
    for record in r1_characterise_read.values():
        if record.description.split(' ')[0].endswith('/1'):
            r1_characterise_read_names.append(
                record.description.split(' ')[0]
            )
        else:
            r1_characterise_read_names.append(
                record.description.split(' ')[0] + '/1'
            )
    for record in r2_characterise_read.values():
        if record.description.split(' ')[0].endswith('/2'):
            r2_characterise_read_names.append(
                record.description.split(' ')[0]
            )
        else:
            r2_characterise_read_names.append(
                record.description.split(' ')[0] + '/2'
            )

    # Sort the header names obtained from characterise_read() for comparison
    r1_characterise_read_names = sorted(r1_characterise_read_names)
    r2_characterise_read_names = sorted(r2_characterise_read_names)

    return (
        r1_load_fastq_names,
        r2_load_fastq_names,
        r1_characterise_read_names,
        r2_characterise_read_names
    )


# Miseq Casava
@pytest.mark.parametrize('parse_fastq_header', [
    (
        'tests/real_fastqs/miseq_casava_R1.fastq',
        'tests/real_fastqs/miseq_casava_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_miseq_casava(parse_fastq_header):
    """
    Test that the headers are parsed correctly for the MiSeq Casava-formatted
    FASTQ files
    """
    (
        r1_load_fastq_names,
        r2_load_fastq_names,
        r1_characterise_read_names,
        r2_characterise_read_names
    ) = parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [
        name.split('/1')[0] for name in r1_load_fastq_names
    ] == [
        name.split('/2')[0] for name in r2_load_fastq_names
    ]


# Miseq Casava SRA
@pytest.mark.parametrize('parse_fastq_header', [
    (
        'tests/real_fastqs/miseq_casava_sra_R1.fastq',
        'tests/real_fastqs/miseq_casava_sra_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_miseq_casava_sra(parse_fastq_header):
    """
    Test that the headers are parsed correctly for the SRA version of the MiSeq
    Casava-formatted FASTQ files
    """
    (
        r1_load_fastq_names,
        r2_load_fastq_names,
        r1_characterise_read_names,
        r2_characterise_read_names
    ) = parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [
        name.split('/1')[0] for name in r1_load_fastq_names
    ] == [
        name.split('/2')[0] for name in r2_load_fastq_names
    ]


# Miseq Casava multilane
@pytest.mark.parametrize('parse_fastq_header', [
    (
        'tests/real_fastqs/miseq_casava_multilane_R1.fastq',
        'tests/real_fastqs/miseq_casava_multilane_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_miseq_casava_multilane(parse_fastq_header):
    """
    Test that the headers are parsed correctly for the multilane version of the
    MiSeq Casava-formatted FASTQ files
    """
    (
        r1_load_fastq_names,
        r2_load_fastq_names,
        r1_characterise_read_names,
        r2_characterise_read_names
    ) = parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [
        name.split('/1')[0] for name in r1_load_fastq_names
    ] == [
        name.split('/2')[0] for name in r2_load_fastq_names
    ]


# Hiseq pre-Casava
@pytest.mark.parametrize('parse_fastq_header', [
    (
        'tests/real_fastqs/hiseq_precasava_R1.fastq',
        'tests/real_fastqs/hiseq_precasava_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_hiseq_precasava(parse_fastq_header):
    """
    Test that the headers are parsed correctly for the HiSeq pre-Casava-
    formatted FASTQ files
    """
    (
        r1_load_fastq_names,
        r2_load_fastq_names,
        r1_characterise_read_names,
        r2_characterise_read_names
    ) = parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [
        name.split('/1')[0] for name in r1_load_fastq_names
    ] == [
        name.split('/2')[0] for name in r2_load_fastq_names
    ]


# Hiseq pre-Casava SRA
@pytest.mark.parametrize('parse_fastq_header', [
    (
        'tests/real_fastqs/hiseq_precasava_sra_R1.fastq',
        'tests/real_fastqs/hiseq_precasava_sra_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_hiseq_precasava_sra(parse_fastq_header):
    """
    Test that the headers are parsed correctly for the SRA version of the HiSeq
    pre-Casava-formatted FASTQ files
    """
    (
        r1_load_fastq_names,
        r2_load_fastq_names,
        r1_characterise_read_names,
        r2_characterise_read_names
    ) = parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [
        name.split('/1')[0] for name in r1_load_fastq_names
    ] == [
        name.split('/2')[0] for name in r2_load_fastq_names
    ]


# Hiseq pre-Casava multilane
@pytest.mark.parametrize('parse_fastq_header', [
    (
        'tests/real_fastqs/hiseq_precasava_multilane_R1.fastq',
        'tests/real_fastqs/hiseq_precasava_multilane_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_hiseq_precasava_multilane(parse_fastq_header):
    """
    Test that the headers are parsed correctly for the multilane version of the
    HiSeq pre-Casava-formatted FASTQ files
    """
    (
        r1_load_fastq_names,
        r2_load_fastq_names,
        r1_characterise_read_names,
        r2_characterise_read_names
    ) = parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [
        name.split('/1')[0] for name in r1_load_fastq_names
    ] == [
        name.split('/2')[0] for name in r2_load_fastq_names
    ]
