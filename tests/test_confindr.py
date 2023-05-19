from confindr_src.methods import *
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
    subprocess.call("confindr.py -i tests/test_samples -o confindr_integration_output -d databases -k -fid '_1' -rid '_2'",
                     shell=True)
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

#-------------------------------
# Test file pattern recognition
#-------------------------------

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
