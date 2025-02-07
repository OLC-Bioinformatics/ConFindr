from confindr_src.methods import *
from Bio import SeqIO
import subprocess
import pytest
import shutil
import csv
import os

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)

#---------------------------
# Test ConFindr integration
#---------------------------

def test_integration():
    correct_contamination_calls = {'SRX5084910_SRR8268082': 'True',
                                   'SRX5084911_SRR8268081': 'False',
                                   'SRX5084914_SRR8268078': 'True',
                                   'SRX5084915_SRR8268077': 'False',
                                   'SRX5084941_SRR8268051': 'True',
                                   'SRX5084940_SRR8268052': 'False',
                                   'SRX5084995_SRR8267997': 'True'
                                   }
    correct_genera = {'SRX5084910_SRR8268082': 'Escherichia',
                      'SRX5084911_SRR8268081': 'Escherichia',
                      'SRX5084914_SRR8268078': 'Salmonella',
                      'SRX5084915_SRR8268077': 'Salmonella',
                      'SRX5084940_SRR8268052': 'Listeria',
                      'SRX5084941_SRR8268051': 'Listeria',
                      'SRX5084995_SRR8267997': 'Salmonella:Citrobacter'}
    
    cmd = [
        "confindr.py",
        "-i", "tests/test_samples",
        "-o", "confindr_integration_output",
        "-d", "databases",
        "-k",
        "-Xmx", "6g"
    ]
    subprocess.call(" ".join(cmd), shell=True)
    with open('confindr_integration_output/confindr_report.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sample = row['Sample']
            if 'cross_contaminated' not in sample:
                assert row['ContamStatus'] == correct_contamination_calls[sample]
                assert row['Genus'] == correct_genera[sample]
            else:
                assert row['ContamStatus'] == correct_contamination_calls[sample]
                genera = row['Genus'].split(':')
                assert 'Salmonella' in genera and 'Escherichia' in genera and 'Listeria' in genera
    shutil.rmtree('confindr_integration_output')
    shutil.rmtree('databases')

#-------------------
# Test dependencies
#-------------------

def test_present_dependency():
    assert dependency_check('ls') is True


def test_nonexistent_dependency():
    assert dependency_check('fake_dependency') is False

#------------------------------------
# Test filename pattern recognition
#------------------------------------

def test_r1_fastqs():
    assert find_paired_reads('tests/fake_fastqs/') == [['tests/fake_fastqs/test_R1.fastq.gz',
                                                        'tests/fake_fastqs/test_R2.fastq.gz']]

def test_1_fastqs():
    assert find_paired_reads('tests/fake_fastqs/', forward_id='_1',
                             reverse_id='_2') == [['tests/fake_fastqs/test_1.fastq.gz',
                                                   'tests/fake_fastqs/test_2.fastq.gz']]

def test_empty_fastqs():
    assert find_paired_reads('tests/fake_fastqs/', forward_id='_asdf', reverse_id='_fdsa') == []


def test_unpaired_fastq():
    assert ['tests/fake_fastqs/test_alone.fastq.gz'] == find_unpaired_reads('tests/fake_fastqs')[2]

#-------------------------
# Test subprocess calling
#-------------------------

def test_run_cmd_success():
    cmd = 'echo asdf'
    out, err = run_cmd(cmd)
    assert out == 'asdf\n'
    assert err == ''

def test_run_cmd_failure_exit_code():
    with pytest.raises(subprocess.CalledProcessError):
        run_cmd('garbagecommandthatdoesnotwork')

#----------------------------------
# Test base count cutoff reporting
#----------------------------------

def test_two_hq_bases_above_threshold():
    assert number_of_bases_above_threshold({'G': 80, 'A': 20}) == 2


def test_just_one_hq_bases_above_threshold():
    assert number_of_bases_above_threshold({'G': 99, 'A': 1}) == 1


def test_two_hq_bases_above_threshold_custom_params():
    assert number_of_bases_above_threshold({'G': 99, 'A': 1}, base_count_cutoff=1) == 2


def test_just_one_hq_base_above_threshold_custom_params():
    assert number_of_bases_above_threshold({'G': 96, 'A': 4}, base_count_cutoff=5) == 1


def test_three_hq_bases_above_threshold():
    assert number_of_bases_above_threshold({'G': 90, 'A': 10, 'T': 10}) == 3


def test_two_out_of_three_hq_bases_above_threshold():
    assert number_of_bases_above_threshold({'G': 90, 'A': 9, 'T': 1}) == 2

#-------------------------------------
# Test base fraction cutoff reporting
#-------------------------------------

def test_two_hq_bases_above_fraction_threshold():
    assert number_of_bases_above_threshold({'G': 80, 'A': 20}, base_fraction_cutoff=0.05) == 2


def test_two_hq_bases_above_fraction_threshold_low_coverage():
    assert number_of_bases_above_threshold({'G': 9, 'A': 1}, base_fraction_cutoff=0.05) == 1


def test_two_hq_bases_above_fraction_threshold_low_coverage_one_base_counts():
    assert number_of_bases_above_threshold({'G': 9, 'A': 1}, base_count_cutoff=1, base_fraction_cutoff=0.05) == 2


def test_just_one_hq_bases_above_fraction_threshold():
    assert number_of_bases_above_threshold({'G': 99, 'A': 1}, base_fraction_cutoff=0.05) == 1


def test_three_hq_bases_above_fraction_threshold():
    assert number_of_bases_above_threshold({'G': 90, 'A': 10, 'T': 10}, base_fraction_cutoff=0.05) == 3


def test_two_out_of_three_hq_bases_above_fraction_threshold():
    assert number_of_bases_above_threshold({'G': 90, 'A': 9, 'T': 1}, base_fraction_cutoff=0.05) == 2


def test_valid_base_fraction_none():
    assert check_valid_base_fraction(None) is True


def test_valid_base_fraction_zero():
    assert check_valid_base_fraction(0.0) is True


def test_valid_base_fraction_one():
    assert check_valid_base_fraction(1.0) is True


def test_valid_base_fraction_between_zero_one():
    assert check_valid_base_fraction(0.2) is True


def test_invalid_base_fraction():
    assert check_valid_base_fraction(1.2) is False

#----------------------------------------
# Test total sequence length calculation
#----------------------------------------

def test_total_length_fasta():
    assert find_total_sequence_length('tests/rmlst.fasta') == 20862

#----------------------
# Test report creation
#----------------------

def test_write_output_creates_file_if_does_not_exist():
    write_output(output_report='tests/confindr_report.csv',
                 sample_name='Test',
                 multi_positions=55,
                 genus='Fakella',
                 total_gene_length=888,
                 database_download_date='ND')
    assert os.path.isfile('tests/confindr_report.csv') is True
    # os.remove('tests/confindr_report.csv')

def test_write_output_appends_if_file_does_exist():
    write_output(output_report='tests/confindr_report.csv',
                 sample_name='Test',
                 multi_positions=55,
                 genus='Fakella',
                 total_gene_length=888,
                 database_download_date='ND')
    with open('tests/confindr_report.csv') as f:
        lines = f.readlines()
    assert len(lines) > 2

#----------------------------------------------
# Test converting base dictionaries to strings
#----------------------------------------------

def test_base_dict_to_string_two_base_descending():
    assert base_dict_to_string({'A': 18, 'C': 3}) == 'A:18;C:3'


def test_base_dict_to_string_two_base_ascending():
    assert base_dict_to_string({'A': 8, 'C': 33}) == 'C:33;A:8'


def test_base_dict_to_string_three_bases():
    assert base_dict_to_string({'A': 5, 'T': 88, 'C': 33}) == 'T:88;C:33;A:5'

#----------------------------------------------
# Check that BBTools memory settings are valid
#----------------------------------------------

def test_valid_xmx_string_gigabytes():
    assert check_acceptable_xmx('20g') is True
    assert check_acceptable_xmx('20G') is True


def test_valid_xmx_string_megabytes():
    assert check_acceptable_xmx('20m') is True
    assert check_acceptable_xmx('20M') is True


def test_valid_xmx_string_kilobytes():
    assert check_acceptable_xmx('550k') is True
    assert check_acceptable_xmx('550K') is True


def test_invalid_xmx_bad_suffix():
    assert check_acceptable_xmx('600u') is False


def test_invalid_xmx_float():
    assert check_acceptable_xmx('2.2G') is False


def test_invalid_xmx_not_an_integer():
    assert check_acceptable_xmx('asdfK') is False

#---------------------------
# Test FASTQ header parsing
#---------------------------
# FASTQ headers can be in different formats, e.g. Casava 1.8, deposited in SRA,
# pre-Casava, etc., and may also be split across multiple lines.
# A Pytest fixture is defined first, and then the unit tests for each different
# FASTQ header format afterwards.

# TODO: Make this fixture more representative of what's actually happening in
# characterise_read() in the future.
# Right now, this is the best we can do without actually calling
# characterise_read(), which currently requires a pileup file to have been 
# generated.
@pytest.fixture
def test_parse_fastq_header(request):
    """
    Fixture to compare the FASTQ headers obtained from two different methods:
    `load_fastq_records` and `SeqIO.parse`. It ensures that the headers are 
    parsed consistently by both methods.

    :param request: A pytest request object that provides access to the parameters 
    passed to the fixture.
    :return: A tuple containing four lists:
        - r1_load_fastq_names: Sorted list of FASTQ header names from `load_fastq_records` for forward reads.
        - r2_load_fastq_names: Sorted list of FASTQ header names from `load_fastq_records` for reverse reads.
        - r1_characterise_read_names: Sorted list of FASTQ header names from `SeqIO.parse` for forward reads.
        - r2_characterise_read_names: Sorted list of FASTQ header names from `SeqIO.parse` for reverse reads.
    """
    r1_path, r2_path = request.param
    # Obtain the parts of the FASTQ headers from load_fastq_records()
    r1_load_fastq = load_fastq_records(
        gz = r1_path,
        paired = True,
        forward = True
    )
    r2_load_fastq = load_fastq_records(
        gz = r2_path,
        paired = True,
        forward = False
    )
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
    r1_characterise_read_names = list()
    r2_characterise_read_names = list()
    # Similar logic as in characterise_read() to extract the read names
    for record in r1_characterise_read.values():
        if record.description.split(' ')[0].endswith('/1'):
            r1_characterise_read_names.append(record.description.split(' ')[0])
        else:
            r1_characterise_read_names.append(record.description.split(' ')[0] + '/1')
    for record in r2_characterise_read.values():
        if record.description.split(' ')[0].endswith('/2'):
            r2_characterise_read_names.append(record.description.split(' ')[0])
        else:
            r2_characterise_read_names.append(record.description.split(' ')[0] + '/2')    
    r1_characterise_read_names = sorted(r1_characterise_read_names)
    r2_characterise_read_names = sorted(r2_characterise_read_names)

    return r1_load_fastq_names, r2_load_fastq_names, r1_characterise_read_names, r2_characterise_read_names

# Miseq Casava
@pytest.mark.parametrize('test_parse_fastq_header', [
    (
        'tests/real_fastqs/miseq_casava_R1.fastq',
        'tests/real_fastqs/miseq_casava_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_miseq_casava(test_parse_fastq_header):
    """
    Test that the headers are parsed correctly for the MiSeq Casava-formatted
    FASTQ files
    """
    r1_load_fastq_names, r2_load_fastq_names, r1_characterise_read_names, r2_characterise_read_names = test_parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [name.split('/1')[0] for name in r1_load_fastq_names] == [name.split('/2')[0] for name in r2_load_fastq_names]

# Miseq Casava SRA
@pytest.mark.parametrize('test_parse_fastq_header', [
    (
        'tests/real_fastqs/miseq_casava_sra_R1.fastq',
        'tests/real_fastqs/miseq_casava_sra_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_miseq_casava_sra(test_parse_fastq_header):
    """
    Test that the headers are parsed correctly for the SRA version of the MiSeq
    Casava-formatted FASTQ files
    """
    r1_load_fastq_names, r2_load_fastq_names, r1_characterise_read_names, r2_characterise_read_names = test_parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [name.split('/1')[0] for name in r1_load_fastq_names] == [name.split('/2')[0] for name in r2_load_fastq_names]

# Miseq Casava multilane
@pytest.mark.parametrize('test_parse_fastq_header', [
    (
        'tests/real_fastqs/miseq_casava_multilane_R1.fastq',
        'tests/real_fastqs/miseq_casava_multilane_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_miseq_casava_multilane(test_parse_fastq_header):
    """
    Test that the headers are parsed correctly for the multilane version of the
    MiSeq Casava-formatted FASTQ files
    """
    r1_load_fastq_names, r2_load_fastq_names, r1_characterise_read_names, r2_characterise_read_names = test_parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [name.split('/1')[0] for name in r1_load_fastq_names] == [name.split('/2')[0] for name in r2_load_fastq_names]

# Hiseq pre-Casava
@pytest.mark.parametrize('test_parse_fastq_header', [
    (
        'tests/real_fastqs/hiseq_precasava_R1.fastq',
        'tests/real_fastqs/hiseq_precasava_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_hiseq_precasava(test_parse_fastq_header):
    """
    Test that the headers are parsed correctly for the HiSeq pre-Casava-
    formatted FASTQ files
    """
    r1_load_fastq_names, r2_load_fastq_names, r1_characterise_read_names, r2_characterise_read_names = test_parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [name.split('/1')[0] for name in r1_load_fastq_names] == [name.split('/2')[0] for name in r2_load_fastq_names]

# Hiseq pre-Casava SRA
@pytest.mark.parametrize('test_parse_fastq_header', [
    (
        'tests/real_fastqs/hiseq_precasava_sra_R1.fastq',
        'tests/real_fastqs/hiseq_precasava_sra_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_hiseq_precasava_sra(test_parse_fastq_header):
    """
    Test that the headers are parsed correctly for the SRA version of the HiSeq
    pre-Casava-formatted FASTQ files
    """
    r1_load_fastq_names, r2_load_fastq_names, r1_characterise_read_names, r2_characterise_read_names = test_parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [name.split('/1')[0] for name in r1_load_fastq_names] == [name.split('/2')[0] for name in r2_load_fastq_names]

# Hiseq pre-Casava multilane
@pytest.mark.parametrize('test_parse_fastq_header', [
    (
        'tests/real_fastqs/hiseq_precasava_multilane_R1.fastq',
        'tests/real_fastqs/hiseq_precasava_multilane_R2.fastq'
    )
], indirect=True)
def test_parse_fastq_header_hiseq_precasava_multilane(test_parse_fastq_header):
    """
    Test that the headers are parsed correctly for the multilane version of the
    HiSeq pre-Casava-formatted FASTQ files
    """
    r1_load_fastq_names, r2_load_fastq_names, r1_characterise_read_names, r2_characterise_read_names = test_parse_fastq_header
    assert r1_load_fastq_names == r1_characterise_read_names
    assert r2_load_fastq_names == r2_characterise_read_names
    assert [name.split('/1')[0] for name in r1_load_fastq_names] == [name.split('/2')[0] for name in r2_load_fastq_names]