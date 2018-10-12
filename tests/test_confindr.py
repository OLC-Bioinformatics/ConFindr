import os
import subprocess
import pytest

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)
from confindr.confindr import *
from Bio import SeqIO

# TODO: Make these far more useful than they currently are.


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
        multi_position_dict = read_contig(contig_name=contig_name,
                                          bamfile_name='tests/contamination.bam',
                                          reference_fasta='tests/rmlst.fasta',
                                          report_file='tests/dummy_report',
                                          quality_cutoff=20,
                                          base_cutoff=2)
        multi_positions += len(multi_position_dict)
    os.remove('tests/dummy_report')
    assert multi_positions == 15


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

def test_just_one_hq_bases_above_fraction_threshold():
    assert number_of_bases_above_threshold({'G': 99, 'A': 1}, base_fraction_cutoff=0.05) == 1

def test_three_hq_bases_above_fraction_threshold():
    assert number_of_bases_above_threshold({'G': 90, 'A': 10, 'T': 10}, base_fraction_cutoff=0.05) == 3

def test_two_out_of_three_hq_bases_above_fraction_threshold():
    assert number_of_bases_above_threshold({'G': 90, 'A': 9, 'T': 1}, base_fraction_cutoff=0.05) == 2
