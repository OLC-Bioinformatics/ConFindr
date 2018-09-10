import shutil
import os
import sys

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)
from confindr.confindr import *

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
    assert 'tests/fake_fastqs/test_alone.fastq.gz' in find_unpaired_reads('tests/fake_fastqs')


def test_rmlst_bait():
    pair = ['tests/mashsippr/O157_R1.fastq.gz', 'tests/mashsippr/O157_R2.fastq.gz']
    actual_result = 'AAAAAAACAGCAAATCCGGTGGTCGTAACAACAATGGCCGTATCACCACTCGTCATATCGGTGGTGGCCA' \
                    'CAAGCAGGCTTACCGTATTGTTGACTTCAAACGCAACAAAGACGGTATCCCGGCAGTTGTTGAACGTCTT' \
                    'GAGTACGATCCGAACCGTTCCGCGAACATCGCGCTGGTTCTGTACAAAGACGGTGAACGCCGTTACATCC' \
                    'TGGCCCCTAAAGGCCTGAAAGCTGGCGACCAGATTCAGTC'
    extract_rmlst_genes(pair, 'tests/bait_fasta.fasta', 'asdf_R1.fasta', 'asdf_R2.fasta')
    thing = SeqIO.read('asdf_R1.fasta', 'fasta')
    assert str(thing.seq) == actual_result
    os.remove('asdf_R1.fasta')
    os.remove('asdf_R2.fasta')



