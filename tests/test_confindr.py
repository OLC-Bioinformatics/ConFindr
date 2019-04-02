import os
import subprocess
import pytest

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)
from confindr_src.confindr import *
from Bio import SeqIO
import subprocess
import shutil
import csv


def test_integration():
    correct_contamination_calls = {'NC_002695_50_HS25_clean': 'False',
                                   'NC_002695_NC_000913_0.05_50_MSv1': 'True',
                                   'NC_002973_50_MSv1_clean': 'False',
                                   'NC_002973_NC_012488_0.05_50_HS25': 'True',
                                   'NC_003198_50_MSv1_clean': 'False',
                                   'NC_003198_NC_003197_0.1_50_HS25': 'True',
                                   'cross_contaminated': 'True'}
    correct_genera = {'NC_002695_50_HS25_clean': 'Escherichia',
                      'NC_002695_NC_000913_0.05_50_MSv1': 'Escherichia',
                      'NC_002973_50_MSv1_clean': 'Listeria',
                      'NC_002973_NC_012488_0.05_50_HS25': 'Listeria',
                      'NC_003198_50_MSv1_clean': 'Salmonella',
                      'NC_003198_NC_003197_0.1_50_HS25': 'Salmonella'}
    subprocess.call('confindr.py -i confindr_integration_tests -o confindr_integration_output -d databases', shell=True)
    with open('confindr_integration_output/confindr_report.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sample = row['Sample']
            if sample != 'cross_contaminated':
                assert row['ContamStatus'] == correct_contamination_calls[sample]
                assert row['Genus'] == correct_genera[sample]
            else:
                assert row['ContamStatus'] == correct_contamination_calls[sample]
                genera = row['Genus'].split(':')
                assert 'Salmonella' in genera and 'Escherichia' in genera and 'Listeria' in genera
    shutil.rmtree('confindr_integration_output')
    shutil.rmtree('databases')


def test_present_dependency():
    assert dependency_check('ls') is True


def test_nonexistent_dependency():
    assert dependency_check('fake_dependency') is False


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


def test_correct_num_multipositions():
    contig_names = list()
    for contig in SeqIO.parse('tests/rmlst.fasta', 'fasta'):
        contig_names.append(contig.id)
    multi_positions = 0
    for contig_name in contig_names:
        multi_position_dict, to_write = read_contig(contig_name=contig_name,
                                                    bamfile_name='tests/contamination.bam',
                                                    reference_fasta='tests/rmlst.fasta',
                                                    quality_cutoff=20,
                                                    base_cutoff=2)
        multi_positions += sum([len(snp_positions) for snp_positions in multi_position_dict.values()])
    assert multi_positions == 24


def test_correct_percent_contam():
    percent_contam, stddev = estimate_percent_contamination('tests/example_contamination.csv')
    assert percent_contam == '18.20'
    assert stddev == '5.89'


def test_run_cmd_success():
    cmd = 'echo asdf'
    out, err = run_cmd(cmd)
    assert out == 'asdf\n'
    assert err == ''


def test_run_cmd_failure_exit_code():
    with pytest.raises(subprocess.CalledProcessError):
        run_cmd('garbagecommandthatdoesnotwork')


# test base_count_cutoff

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

# test base_fraction_cutoff

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


def test_total_length_fasta():
    assert find_total_sequence_length('tests/rmlst.fasta') == 20862


def test_write_output_creates_file_if_does_not_exist():
    write_output(output_report='tests/file_that_does_not_exist.csv',
                 sample_name='Test',
                 multi_positions=55,
                 genus='Fakella',
                 percent_contam=22.2,
                 contam_stddev=1.1,
                 total_gene_length=888,
                 database_download_date='NA')
    assert os.path.isfile('tests/file_that_does_not_exist.csv') is True
    os.remove('tests/file_that_does_not_exist.csv')


def test_write_output_appends_if_file_does_exist():
    write_output(output_report='tests/confindr_report.csv',
                 sample_name='Test',
                 multi_positions=55,
                 genus='Fakella',
                 percent_contam=22.2,
                 contam_stddev=1.1,
                 total_gene_length=888,
                 database_download_date='NA')
    with open('tests/confindr_report.csv') as f:
        lines = f.readlines()
    assert len(lines) > 2


def test_base_dict_to_string_two_base_descending():
    assert base_dict_to_string({'A': 18, 'C': 3}) == 'A:18;C:3'


def test_base_dict_to_string_two_base_ascending():
    assert base_dict_to_string({'A': 8, 'C': 33}) == 'C:33;A:8'


def test_base_dict_to_string_three_bases():
    assert base_dict_to_string({'A': 5, 'T': 88, 'C': 33}) == 'T:88;C:33;A:5'


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

