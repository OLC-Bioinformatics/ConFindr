#!/usr/bin/env python3

"""
Tests for helper utilities and wrapper command construction in ConFindr.
"""

# Standard imports
import argparse
import builtins
import gzip
import json
import os
import subprocess
import tarfile
from pathlib import Path
from types import SimpleNamespace

# Third-party imports
import pytest

# Local imports
import confindr_src.create_genus_specific_db as cgdb
import confindr_src.database_setup as dbsetup
import confindr_src.methods as methods
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from confindr_src.methods import (
    _format_seconds,
    base_dict_to_string,
    count_fastq_reads,
    dependency_check,
    downsample_reads,
    estimate_genome_size,
    find_genus_specific_allele_list,
    find_total_sequence_length,
    find_unpaired_reads,
    load_fastq_records,
    setup_allelespecific_database,
    _valid_downsample_depth,
    check_valid_base_fraction,
    check_acceptable_xmx,
    recommend_xmx,
    write_to_logfile,
    find_rmlst_type,
    count_multibase_positions,
)
from confindr_src.wrappers import bbtools, mash


def test_base_dict_to_string_orders_by_count():
    base_dict = {'A': 2, 'C': 10, 'G': 5}
    assert base_dict_to_string(base_dict=base_dict) == 'C:10;G:5;A:2'


def test_find_total_sequence_length(tmp_path):
    fasta = tmp_path / 'test.fasta'
    fasta.write_text('>seq1\nAAAA\n>seq2\nTTTTT\n', encoding='utf-8')
    assert find_total_sequence_length(fasta_file=str(fasta)) == 9


def write_fastq(path, records):
    with open(path, 'w', encoding='utf-8') as handle:
        for i, (seq, qual) in enumerate(records, start=1):
            handle.write(f'@read{i}\n{seq}\n+\n{qual}\n')


def test_load_fastq_records_paired_forward(tmp_path):
    fastq = tmp_path / 'reads_R1.fastq'
    write_fastq(fastq, [('ACTG', 'IIII'), ('TGCA', 'IIII')])
    records = load_fastq_records(gz=str(fastq), paired=True, forward=True)
    assert 'read1/1' in records
    assert 'read2/1' in records
    assert records['read1/1'].seq == 'ACTG'


def test_load_fastq_records_paired_reverse(tmp_path):
    fastq = tmp_path / 'reads_R2.fastq'
    write_fastq(fastq, [('ACTG', 'IIII')])
    records = load_fastq_records(gz=str(fastq), paired=True, forward=False)
    assert 'read1/2' in records
    assert records['read1/2'].seq == 'ACTG'


def test_load_fastq_records_gz_path(tmp_path):
    fastq = tmp_path / 'reads_R1.fastq.gz'
    with gzip.open(fastq, 'wt', encoding='utf-8') as handle:
        handle.write('@read1\nACTG\n+\nIIII\n')
    records = load_fastq_records(gz=str(fastq), paired=True, forward=True)
    assert 'read1/1' in records
    assert records['read1/1'].seq == 'ACTG'


@pytest.mark.parametrize('value,expected', [
    ('10', 10),
    ('50', 50),
    ('100', 100),
])
def test_valid_downsample_depth(value, expected):
    assert _valid_downsample_depth(value) == expected


@pytest.mark.parametrize('value', ['-1', 'abc', '1.5'])
def test_invalid_downsample_depth(value):
    with pytest.raises((argparse.ArgumentTypeError, ValueError, TypeError)):
        _valid_downsample_depth(value)


def test_check_valid_base_fraction():
    assert check_valid_base_fraction(base_fraction=0.0) is True
    assert check_valid_base_fraction(base_fraction=0.5) is True
    assert check_valid_base_fraction(base_fraction=1.0) is True
    assert check_valid_base_fraction(base_fraction=-0.1) is False
    assert check_valid_base_fraction(base_fraction=1.1) is False


def test_check_acceptable_xmx():
    assert check_acceptable_xmx(xmx_string='4g') is True
    assert check_acceptable_xmx(xmx_string='512m') is True
    assert not check_acceptable_xmx(xmx_string='5x')
    assert not check_acceptable_xmx(xmx_string='4.5g')


@pytest.mark.parametrize('available_bytes,expected', [
    (8 * 1024**3, '6g'),
    (16 * 1024**3, '12g'),
])
def test_recommend_xmx(monkeypatch, available_bytes, expected):
    class FakeVmem:
        def __init__(self, available):
            self.available = available
    monkeypatch.setattr('confindr_src.methods.psutil.virtual_memory', lambda: FakeVmem(available_bytes))
    assert recommend_xmx(fraction=0.75, round_gb=True) == expected


def test_count_multibase_positions_excludes_meta():
    multibase_dicts = [
        {'gene1': {10: {'paired': {}}, '_gene_stats': {}}},
        {'gene2': {20: {'paired': {}}, 30: {'paired': {}}, '_gene_stats': {}}}
    ]
    assert count_multibase_positions(multibase_dict_list=multibase_dicts) == 3


def test_find_rmlst_type_writes_report(tmp_path):
    kma_report = tmp_path / 'kma.tsv'
    rmlst_report = tmp_path / 'rmlst.tsv'
    kma_report.write_text(
        '#Template\tScore\n'
        'abcZ_1\t10\n'
        'abcZ_2\t5\n'
        'adk_4\t8\n'
        'adk_3\t12\n',
        encoding='utf-8'
    )
    result = find_rmlst_type(
        kma_report=str(kma_report),
        rmlst_report=str(rmlst_report)
    )
    assert result == ['abcZ_1', 'adk_3']
    assert rmlst_report.read_text().splitlines()[1:] == ['abcZ\t1', 'adk\t3']


def test_format_seconds_and_quality():
    assert _format_seconds(s=3661) == '1:01:01'
    assert _format_seconds(s=59) == '00:59'
    assert _format_seconds(s='not-a-number') == 'N/A'
    assert _format_seconds(s=float('nan')) == 'N/A'


def test_format_seconds_unroundable_input_returns_na():
    class BadNumber:
        def __round__(self):
            raise TypeError('cannot round')

    assert methods._format_seconds(s=BadNumber()) == 'N/A'


def test_estimate_genome_size_variants():
    assert estimate_genome_size(genus='Escherichia') == 4600000
    assert estimate_genome_size(genus='Chlamydophila') == 1000000
    assert estimate_genome_size(genus='') == 4000000
    assert estimate_genome_size(genus='   ') == 4000000


def test_estimate_mean_read_length_handles_eof_and_sample_reads(tmp_path):
    fastq = tmp_path / 'reads.fastq'
    fastq.write_text('@r1\n', encoding='utf-8')
    assert methods.estimate_mean_read_length(fastq_path=str(fastq), sample_reads=100) == 50

    fastq2 = tmp_path / 'reads2.fastq'
    fastq2.write_text('@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n', encoding='utf-8')
    assert methods.estimate_mean_read_length(fastq_path=str(fastq2), sample_reads=1) == 50


def test_parse_bam_bytes_contig_name_decode_error(monkeypatch):
    class FakeAlignmentFile:
        def __init__(self, *args, **kwargs):
            pass
        def pileup(self, *args, **kwargs):
            return []
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    bamfile, pileup = methods.parse_bam(
        bamfile_name='dummy.bam',
        contig_name=b'\xff',
        pysam_fasta=None
    )
    assert pileup == []


def test_init_fastq_index_and_get_fastq_record(monkeypatch, tmp_path):
    fwd = tmp_path / 'reads_R1.fastq'
    rev = tmp_path / 'reads_R2.fastq'
    write_fastq(fwd, [('ACGT', 'IIII')])
    write_fastq(rev, [('TGCA', 'IIII')])
    tmpdir = tmp_path / 'idx'
    tmpdir.mkdir()

    def fake_index_db(db, path, fmt):
        assert fmt == 'fastq'
        if db.endswith('.fwd.idx'):
            return {'read1': 'forward'}
        return {'read1': 'reverse'}

    monkeypatch.setattr(methods.SeqIO, 'index_db', fake_index_db)
    prev_state = methods._FASTQ_INDEX_STATE.copy()
    try:
        methods._init_fastq_index(
            fwd_path=str(fwd),
            rev_path=str(rev),
            paired=True,
            tmpdir=str(tmpdir)
        )
        assert methods._FASTQ_INDEX_STATE['paired'] is True
        assert methods._get_fastq_record('read1/1') == 'forward'
        assert methods._get_fastq_record('read1/2') == 'reverse'
        assert methods._get_fastq_record('read1') == 'forward'
    finally:
        methods._FASTQ_INDEX_STATE.clear()
        methods._FASTQ_INDEX_STATE.update(prev_state)


def test_find_multibase_positions_without_fractional_cutoff_and_with_base_support():
    filtered_read_dict = {
        'congruent_SNV': {'C': 4},
        'forward_SNV_reverse_UM_QF': {'C': 2},
        'reverse_SNV_forward_UM_QF': {'C': 2},
        'forward_ref_reverse_UM_QF': {'A': 1},
        'reverse_ref_forward_UM_QF': {'A': 1},
    }
    base_support = {
        'C': {
            'quals': [30, 30],
            'forward': 1,
            'reverse': 1,
            'positions': [1, 2],
            'forward_positions': [1],
            'reverse_positions': [2],
            'mapqs': [20, 20],
        },
        'A': {
            'quals': [30],
            'forward': 1,
            'reverse': 0,
            'positions': [0],
            'forward_positions': [0],
            'reverse_positions': [],
            'mapqs': [20],
        }
    }
    snv_dict, passing_snv_dict, total_coverage, position_stats = methods.find_multibase_positions(
        ref_base='A',
        filtered_read_dict=filtered_read_dict,
        base_cutoff=2,
        base_fraction_cutoff=0.0,
        base_support=base_support
    )
    assert total_coverage == 10
    assert passing_snv_dict['congruent']['C'] == 4
    assert 'C' in position_stats


def test_read_contig_fasta_reference_fallback(monkeypatch, tmp_path):
    fasta = tmp_path / 'ref.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        calls = 0
        def __init__(self, path):
            FakeFastaFile.calls += 1
            self.references = [] if FakeFastaFile.calls == 1 else [b'abc_1']
        def close(self):
            pass

    def fake_faidx(path):
        Path(str(path) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods.pysam, 'faidx', fake_faidx)
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))

    result, report = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records={'abc_1': SeqRecord(Seq('ACGT'), id='abc_1')},
        fastq_records={},
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=1,
        fasta=False,
        nanopore=False
    )
    assert 'abc_1' in result
    assert result['abc_1']['_gene_stats']['num_sig_positions'] == 0
    assert report == ''


def test_read_contig_bytes_contig_name_decodes(monkeypatch, tmp_path):
    fasta = tmp_path / 'ref2.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            self.references = [b'abc_1']
        def close(self):
            pass

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))

    result, report = methods.read_contig(
        contig_name=b'abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records={'abc_1': SeqRecord(Seq('ACGT'), id='abc_1')},
        fastq_records={},
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=1,
        fasta=False,
        nanopore=False
    )
    assert report == ''
    assert 'abc_1' in result


def test_characterise_read_aggregates_existing_category_counts(tmp_path):
    def make_read(qname, seq, is_read1):
        return SimpleNamespace(
            query_position=0,
            alignment=SimpleNamespace(
                qname=qname,
                query_sequence=seq,
                query_alignment_end=len(seq),
                is_read1=is_read1,
                is_read2=not is_read1,
                is_paired=True,
                mate_is_unmapped=False,
                mapping_quality=30
            )
        )

    reads = [
        make_read('read1', 'C', True),
        make_read('read1', 'C', False),
        make_read('read2', 'C', True),
        make_read('read2', 'C', False),
    ]
    column = SimpleNamespace(pos=0, reference_name='gene1', pileups=reads)
    fastq_records = {
        'read1/1': SeqRecord(Seq('C'), id='read1/1', letter_annotations={'phred_quality': [30]}),
        'read1/2': SeqRecord(Seq('C'), id='read1/2', letter_annotations={'phred_quality': [30]}),
        'read2/1': SeqRecord(Seq('C'), id='read2/1', letter_annotations={'phred_quality': [30]}),
        'read2/2': SeqRecord(Seq('C'), id='read2/2', letter_annotations={'phred_quality': [30]}),
    }

    filtered, qualities, base_support = methods.characterise_read(
        column=column,
        reference_sequence='AC',
        fastq_records=fastq_records,
        quality_cutoff=20,
        min_quality=15,
        fasta=False,
        nanopore=False
    )
    assert filtered['congruent_SNV'] == {'C': 4}
    assert qualities == [30, 30, 30, 30]
    assert base_support['C']['forward'] == 2
    assert base_support['C']['reverse'] == 2


def test_find_multibase_positions_with_fractional_cutoff(monkeypatch):
    filtered_read_dict = {
        'congruent_SNV': {'C': 4},
        'forward_SNV_reverse_UM_QF': {'C': 4},
        'reverse_SNV_forward_UM_QF': {'G': 4},
        'congruent_ref': {'A': 4}
    }
    base_support = {
        'C': {
            'quals': [30, 30, 30, 30],
            'forward': 2,
            'reverse': 2,
            'positions': [1, 2, 3, 4],
            'forward_positions': [1, 2],
            'reverse_positions': [3, 4],
            'mapqs': [20, 20, 20, 20],
        },
        'G': {
            'quals': [30, 30, 30, 30],
            'forward': 2,
            'reverse': 2,
            'positions': [1, 2, 3, 4],
            'forward_positions': [1, 2],
            'reverse_positions': [3, 4],
            'mapqs': [20, 20, 20, 20],
        },
        'A': {
            'quals': [30, 30, 30, 30],
            'forward': 2,
            'reverse': 2,
            'positions': [1, 2, 3, 4],
            'forward_positions': [1, 2],
            'reverse_positions': [3, 4],
            'mapqs': [20, 20, 20, 20],
        }
    }
    snv_dict, passing_snv_dict, total_coverage, position_stats = methods.find_multibase_positions(
        ref_base='A',
        filtered_read_dict=filtered_read_dict,
        base_cutoff=2,
        base_fraction_cutoff=0.2,
        base_support=base_support
    )
    assert passing_snv_dict['congruent']['C'] == 4
    assert passing_snv_dict['forward']['C'] == 4
    assert passing_snv_dict['reverse']['G'] == 4
    assert total_coverage == 16
    assert 'C' in position_stats


def test_position_entry_passes_probabilistic_gating():
    position_stats = {
        'A': {'q_value': 0.01, 'strand_p': 0.5, 'pos_p': 0.5},
        'C': {'q_value': 0.1, 'strand_p': 0.5, 'pos_p': 0.5}
    }
    assert methods._position_entry_passes_probabilistic_gating(
        position_stats=position_stats,
        q_threshold=0.05,
        strand_p_threshold=0.1,
        pos_p_threshold=0.1
    ) is True


def test_determine_cutoff_empty_qualities_returns_minimum():
    assert methods.determine_cutoff(
        qualities=[],
        reference_sequence='ACGT',
        base_cutoff=1,
        error_cutoff=1.0
    ) == (1, 0.0, 0.0)


def test_poisson_binomial_tail_degenerate_cases():
    assert methods.poisson_binomial_tail(k=0, p_list=[]) == 1.0
    assert methods.poisson_binomial_tail(k=1, p_list=[0.0, 0.0]) == 0.0


def test_fisher_two_sided_p_for_equal_counts():
    assert methods.fisher_two_sided_p(a=1, b=1, c=1, d=1) == pytest.approx(1.0)


def test_get_fastq_record_strips_description_fields(monkeypatch):
    prev_state = methods._FASTQ_INDEX_STATE.copy()
    try:
        methods._FASTQ_INDEX_STATE.clear()
        methods._FASTQ_INDEX_STATE.update({
            'fwd': {'read1': 'forward'},
            'rev': {'read1': 'reverse'},
            'paired': True
        })
        assert methods._get_fastq_record('read1/1') == 'forward'
        assert methods._get_fastq_record('read1/2') == 'reverse'
        assert methods._get_fastq_record('read1  extra') == 'forward'
    finally:
        methods._FASTQ_INDEX_STATE.clear()
        methods._FASTQ_INDEX_STATE.update(prev_state)


def test_count_fastq_reads_uses_wc_and_zcat(monkeypatch, tmp_path):
    content = '@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nIIII\n'
    fastq = tmp_path / 'reads.fastq'
    fastq.write_text(content, encoding='utf-8')
    gz = tmp_path / 'reads.fastq.gz'
    with gzip.open(gz, 'wt', encoding='utf-8') as handle:
        handle.write(content)

    commands = []

    def fake_run_cmd(cmd):
        commands.append(cmd)
        return '8\n', ''

    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    assert count_fastq_reads(fastq_path=str(fastq)) == 2
    assert 'wc -l' in commands[0]
    assert count_fastq_reads(fastq_path=str(gz)) == 2
    assert 'zcat' in commands[1]


def test_estimate_mean_read_length_plain_and_gz(tmp_path):
    content = '@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nIIII\n'
    fastq = tmp_path / 'reads.fastq'
    fastq.write_text(content, encoding='utf-8')
    gz = tmp_path / 'reads.fastq.gz'
    with gzip.open(gz, 'wt', encoding='utf-8') as handle:
        handle.write(content)

    assert methods.estimate_mean_read_length(fastq_path=str(fastq)) == 50
    assert methods.estimate_mean_read_length(fastq_path=str(gz)) == 50


def test_estimate_mean_read_length_returns_default_on_missing_file(tmp_path):
    assert methods.estimate_mean_read_length(fastq_path=str(tmp_path / 'missing.fastq')) == 150


def test_find_paired_and_unpaired_reads(tmp_path):
    r1 = tmp_path / 'sample_R1.fastq'
    r2 = tmp_path / 'sample_R2.fastq'
    other = tmp_path / 'single.fastq'
    r1.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    r2.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    other.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    pairs = methods.find_paired_reads(fastq_directory=str(tmp_path))
    assert [str(r1), str(r2)] in pairs

    unpaired = methods.find_unpaired_reads(fastq_directory=str(tmp_path))
    assert [str(other)] in unpaired


def test_run_cmd_success_and_failure():
    out, err = methods.run_cmd(cmd='echo hello')
    assert 'hello' in out
    assert err == ''
    with pytest.raises(Exception):
        methods.run_cmd(cmd='false')


def test_downsample_reads_paired_and_unpaired(monkeypatch, tmp_path):
    log = tmp_path / 'downsample.log'
    commands = []

    def fake_run_cmd(cmd):
        commands.append(cmd)
        return 'stdout', 'stderr'

    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)

    paired_output = downsample_reads(
        pair=['/input/reads_R1.fastq.gz', '/input/reads_R2.fastq.gz'],
        sample_tmp_dir=str(tmp_path),
        sample_name='sample',
        target_reads=100,
        log=str(log),
        xmx='4g',
        seed=42
    )
    assert len(paired_output) == 2
    assert 'reformat.sh' in commands[0]
    assert 'samplereadstarget=100' in commands[0]
    assert 'sampleseed=42' in commands[0]
    assert str(paired_output[0]).endswith('_downsampled_R1.fastq.gz')
    assert 'Command used:' in log.read_text()

    commands.clear()
    single_log = tmp_path / 'downsample_single.log'
    single_output = downsample_reads(
        pair=['/input/reads.fastq.gz'],
        sample_tmp_dir=str(tmp_path),
        sample_name='single',
        target_reads=50,
        log=str(single_log),
        xmx='2g'
    )
    assert len(single_output) == 1
    assert 'in=/input/reads.fastq.gz' in commands[0]
    assert 'samplereadstarget=50' in commands[0]


def test_write_to_logfile(tmp_path):
    log = tmp_path / 'log.txt'
    write_to_logfile(
        logfile=str(log),
        out='out data',
        err='err data',
        cmd='echo hello'
    )
    content = log.read_text()
    assert 'echo hello' in content
    assert 'STDOUT: out data' in content
    assert 'STDERR: err data' in content


def test_kwargs_to_string_and_bbtools_command_generation(monkeypatch, tmp_path):
    assert bbtools.kwargs_to_string({'xmx': '4g', 'threads': 2}) == ' xmx=4g threads=2'

    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    commands = []

    def fake_run_subprocess(cmd):
        commands.append(cmd)
        return 'out', 'err'

    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)

    _, _, cmd = bbtools.bbmap(
        reference=str(tmp_path / 'ref.fasta'),
        forward_in=str(forward),
        out_bam=str(tmp_path / 'out.bam'),
        returncmd=True
    )
    assert 'in2=' in cmd
    assert 'out=' in cmd

    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda args: b'/usr/bin/bbduk.sh\n')
    _, _, cmd = bbtools.bbduk_trim(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'trimmed_R1.fastq'),
        returncmd=True
    )
    assert 'in1=' in cmd
    assert 'in2=' in cmd
    assert 'out1=' in cmd
    assert 'out2=' in cmd

    _, _, cmd = bbtools.bbduk_bait(
        reference=str(tmp_path / 'ref.fasta'),
        forward_in=str(forward),
        forward_out=str(tmp_path / 'baited_R1.fastq'),
        returncmd=True
    )
    assert 'ref=' in cmd
    assert 'in2=' in cmd
    assert 'outm=' in cmd


def test_bbtools_tadpole_and_bbnorm_and_bbmerge(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_run_subprocess(cmd):
        return 'out', 'err'

    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)

    _, _, tadpole_cmd = bbtools.tadpole(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'tadpole_R1.fastq'),
        returncmd=True
    )
    assert 'tadpole.sh' in tadpole_cmd

    _, _, bbnorm_cmd = bbtools.bbnorm(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'bbnorm_R1.fastq'),
        returncmd=True
    )
    assert 'bbnorm.sh' in bbnorm_cmd

    _, _, bbmerge_cmd = bbtools.bbmerge(
        forward_in=str(forward),
        merged_reads=str(tmp_path / 'merged.fastq'),
        returncmd=True
    )
    assert 'bbmerge.sh' in bbmerge_cmd


def test_bbtools_bbduk_filter_and_subsample_and_validate(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_run_subprocess(cmd):
        return 'out', 'err'

    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)

    _, _, filter_cmd = bbtools.bbduk_filter(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'filtered_R1.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh' in filter_cmd

    _, _, subsample_cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'subsampled_R1.fastq'),
        num_bases=100,
        returncmd=True
    )
    assert 'reformat.sh' in subsample_cmd

    _, _, validate_cmd = bbtools.validate_reads(
        forward_in=str(forward),
        returncmd=True
    )
    assert 'reformat.sh' in validate_cmd

    _, _, reformat_cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'reformat_R1.fastq'),
        returncmd=True
    )
    assert 'reformat.sh' in reformat_cmd


def test_bbtools_tadpole_reverse_out_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: (_ for _ in ()).throw(AssertionError('should not run')))

    with pytest.raises(ValueError):
        bbtools.tadpole(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq'),
            reverse_in=str(reverse),
            returncmd=True
        )


def test_bbtools_tadpole_skips_when_output_exists(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    output = tmp_path / 'tadpole_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output.write_text('existing', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: (_ for _ in ()).throw(AssertionError('should not run')))

    out, err, cmd = bbtools.tadpole(
        forward_in=str(forward),
        forward_out=str(output),
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'tadpole.sh' in cmd


def test_bbtools_bbnorm_reverse_out_error(tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.bbnorm(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq'),
            reverse_in=str(reverse),
            returncmd=True
        )


def test_bbtools_tadpole_auto_reverse_out_missing_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.tadpole(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq')
        )


def test_bbtools_bbnorm_auto_reverse_out_missing_error(tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.bbnorm(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq')
        )


def test_bbtools_bbduk_trim_auto_reverse_out_missing_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda args: b'/usr/bin/bbduk.sh\n')

    with pytest.raises(ValueError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq')
        )


def test_bbtools_bbduk_bait_auto_reverse_out_missing_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.bbduk_bait(
            reference='ref.fasta',
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq')
        )


def test_bbtools_bbduk_filter_auto_reverse_out_missing_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.bbduk_filter(
            reference='ref.fasta',
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq')
        )


def test_bbtools_subsample_reads_auto_reverse_out_missing_error(tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.subsample_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq'),
            num_bases=100
        )


def test_bbtools_reformat_reads_auto_reverse_out_missing_error(tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.reformat_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq')
        )


def test_bbtools_repair_reads_auto_reverse_out_missing_error(tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.repair_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq')
        )


def test_bbtools_tadpole_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.tadpole(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'tadpole_R1.fastq')
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_bbnorm_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.bbnorm(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'bbnorm_R1.fastq')
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_bbmerge_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.bbmerge(
        forward_in=str(forward),
        merged_reads=str(tmp_path / 'merged.fastq')
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_bbduk_bait_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.bbduk_bait(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'baited_R1.fastq')
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_bbduk_filter_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.bbduk_filter(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'filtered_R1.fastq')
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_seal_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.seal(
        reference='ref.fasta',
        forward_in=str(forward),
        output_file=str(tmp_path / 'stats.txt')
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_kmercountexact_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.kmercountexact(
        forward_in=str(forward)
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_subsample_reads_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'subsampled_R1.fastq'),
        num_bases=100
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_reformat_reads_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'reformat_R1.fastq')
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_repair_reads_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.repair_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'repaired_R1.fastq')
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_validate_reads_returns_out_err_with_run_subprocess(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err = bbtools.validate_reads(
        forward_in=str(forward)
    )
    assert out == 'out'
    assert err == 'err'


def test_bbtools_bbmerge_skips_when_output_exists(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    merged = tmp_path / 'merged.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    merged.write_text('existing', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: (_ for _ in ()).throw(AssertionError('should not run')))

    out, err, cmd = bbtools.bbmerge(
        forward_in=str(forward),
        merged_reads=str(merged),
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'bbmerge.sh' in cmd


def test_bbtools_bbduk_trim_raises_reverse_out_missing_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda args: b'/usr/bin/bbduk.sh\n')

    with pytest.raises(ValueError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq'),
            reverse_in=str(reverse)
        )


def test_bbtools_bbduk_bait_reverse_out_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: (_ for _ in ()).throw(AssertionError('should not run')))

    with pytest.raises(ValueError):
        bbtools.bbduk_bait(
            reference='ref.fasta',
            forward_in=str(forward),
            forward_out=str(tmp_path / 'baited.fastq'),
            reverse_in=str(reverse)
        )


def test_bbtools_bbduk_filter_reverse_out_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: (_ for _ in ()).throw(AssertionError('should not run')))

    with pytest.raises(ValueError):
        bbtools.bbduk_filter(
            reference='ref.fasta',
            forward_in=str(forward),
            forward_out=str(tmp_path / 'filtered.fastq'),
            reverse_in=str(reverse)
        )


def test_bbtools_seal_uses_reverse_in_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err, cmd = bbtools.seal(
        reference='ref.fasta',
        forward_in=str(forward),
        output_file=str(tmp_path / 'stats.txt'),
        reverse_in=str(reverse),
        returncmd=True
    )
    assert 'in2=' in cmd
    assert 'rpkm=' in cmd


def test_bbtools_kmercountexact_uses_reverse_in_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err, cmd = bbtools.kmercountexact(
        forward_in=str(forward),
        reverse_in=str(reverse),
        returncmd=True
    )
    assert 'in2=' in cmd
    assert 'kmercountexact.sh' in cmd


def test_bbtools_subsample_reads_reverse_in_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.subsample_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'subsampled.fastq'),
            num_bases=100,
            reverse_in=str(reverse)
        )


def test_bbtools_subsample_reads_skips_when_output_exists(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output = tmp_path / 'subsampled.fastq'
    output.write_text('existing', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: (_ for _ in ()).throw(AssertionError('should not run')))

    out, err, cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(output),
        num_bases=100,
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'reformat.sh' in cmd


def test_bbtools_validate_reads_with_reverse_in(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: ('out', 'err'))

    out, err, cmd = bbtools.validate_reads(
        forward_in=str(forward),
        reverse_in=str(reverse),
        returncmd=True
    )
    assert 'vpair' in cmd


def test_bbtools_reformat_reads_reverse_in_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.reformat_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'reformat.fastq'),
            reverse_in=str(reverse)
        )


def test_bbtools_reformat_reads_skips_when_output_exists(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output = tmp_path / 'reformat.fastq'
    output.write_text('existing', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: (_ for _ in ()).throw(AssertionError('should not run')))

    out, err, cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(output),
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'reformat.sh' in cmd


def test_bbtools_repair_reads_reverse_in_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.repair_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'repaired.fastq'),
            reverse_in=str(reverse)
        )


def test_bbtools_repair_reads_skips_when_output_exists(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output = tmp_path / 'repaired.fastq'
    output.write_text('existing', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda cmd: (_ for _ in ()).throw(AssertionError('should not run')))

    out, err, cmd = bbtools.repair_reads(
        forward_in=str(forward),
        forward_out=str(output),
        reverse_in=str(tmp_path / 'reads_R2.fastq'),
        reverse_out=str(output.parent / 'repaired_R2.fastq'),
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'repair.sh' in cmd


def test_bbtools_genome_size_haploid_false(tmp_path):
    peaks = tmp_path / 'peaks.txt'
    peaks.write_text('#genome_size 900000\n', encoding='utf-8')
    assert bbtools.genome_size(str(peaks), haploid=False) == 900000


def test_bbtools_backward_reverse_out_error(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda args: b'/usr/bin/bbduk.sh\n')
    with pytest.raises(ValueError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'out.fastq'),
            reverse_in='reads_R2.fastq'
        )


def test_dependency_check(monkeypatch):
    monkeypatch.setattr(methods.shutil, 'which', lambda dep: '/usr/bin/python' if dep == 'python' else None)
    assert dependency_check(dependency='python') is True
    assert dependency_check(dependency='nonexistent') is False


def test_find_genus_specific_allele_list_and_setup_database(tmp_path):
    profiles = tmp_path / 'profiles.txt'
    profiles.write_text(
        'Listeria:abc_1,def_2,\n'
        'Other:xyz_3,\n',
        encoding='utf-8'
    )
    assert find_genus_specific_allele_list(
        profiles_file=str(profiles),
        target_genus='Listeria'
    ) == ['abc_1', 'def_2']

    database_folder = tmp_path
    fasta_file = database_folder / 'rMLST_combined.fasta'
    fasta_file.write_text(
        '>abc_1\nACGT\n>def_2\nTGCA\n',
        encoding='utf-8'
    )
    output_database = tmp_path / 'subset.fasta'
    setup_allelespecific_database(
        fasta_file=str(output_database),
        database_folder=str(database_folder),
        allele_list=['abc_1', 'missing_1']
    )
    output_text = output_database.read_text(encoding='utf-8')
    assert '>abc_1' in output_text
    assert '>def_2' not in output_text


def test_setup_confindr_database_combines_fake_loci_and_profiles(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')
    output_folder = tmp_path / 'db'
    output_folder.mkdir()

    class FakeRmlstRest:
        def __init__(self, consumer_secret_file, output_folder, unverified=False):
            self.output_folder = output_folder

        def get_request_token(self):
            return None

        def get_access_token(self):
            return None

        def get_session_token(self):
            return None

        def get_loci_and_scheme_url(self):
            return None

        def download_loci(self):
            with open(
                os.path.join(self.output_folder, 'BACT000001.tfa'),
                'w', encoding='utf-8'
            ) as f:
                f.write('>abc_1\nACGT\n')

        def download_profile(self):
            header = 'genus'
            for i in range(1, 66):
                header += f'\tBACT{i:06d}'
            row = 'Escherichia' + '\t' + '\t'.join(['1'] + ['N'] * 64)
            with open(
                os.path.join(self.output_folder, 'profiles.txt'),
                'w', encoding='utf-8'
            ) as f:
                f.write(header + '\n' + row + '\n')

    monkeypatch.setattr(dbsetup, 'RmlstRest', FakeRmlstRest)
    monkeypatch.setattr(dbsetup, 'index', lambda **kwargs: None)

    dbsetup.setup_confindr_database(
        output_folder=str(output_folder),
        consumer_secret=str(secret_file),
        index_databases=False,
        unverified=False
    )

    assert (output_folder / 'rMLST_combined.fasta').exists()
    assert (output_folder / 'gene_allele.txt').exists()
    combined_text = (output_folder / 'rMLST_combined.fasta').read_text(encoding='utf-8')
    assert '>abc_1' in combined_text


def test_rmlst_rest_get_loci_and_scheme_url_json(monkeypatch, tmp_path):
    class FakeResponse:
        def __init__(self, status_code, headers, payload):
            self.status_code = status_code
            self.headers = headers
            self._payload = payload
        def json(self):
            return self._payload
        @property
        def text(self):
            return json.dumps(self._payload)

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(
                200,
                {'content-type': 'application/json'},
                {'loci': ['url1'], 'schemes': 'profile_url'}
            )

    secret_path = tmp_path / 'secret.txt'
    secret_path.write_text('key\nsecret\n', encoding='utf-8')
    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)
    rest = dbsetup.RmlstRest(consumer_secret_file=str(secret_path), output_folder=str(tmp_path))
    rest.session_token = 'session'
    rest.session_secret = 'secret'
    rest.get_loci_and_scheme_url()
    assert rest.loci == ['url1']
    assert rest.profile == 'profile_url'


def test_rmlst_rest_get_loci_and_scheme_url_text(monkeypatch, tmp_path):
    class FakeResponse:
        def __init__(self, status_code, headers, payload):
            self.status_code = status_code
            self.headers = headers
            self._payload = payload
        def json(self):
            return self._payload
        @property
        def text(self):
            return json.dumps(self._payload)

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(
                200,
                {'content-type': 'text/plain'},
                {'loci': ['url2'], 'schemes': 'profile_text'}
            )

    secret_path = tmp_path / 'secret.txt'
    secret_path.write_text('key\nsecret\n', encoding='utf-8')
    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)
    rest = dbsetup.RmlstRest(consumer_secret_file=str(secret_path), output_folder=str(tmp_path))
    rest.session_token = 'session'
    rest.session_secret = 'secret'
    rest.get_loci_and_scheme_url()
    assert rest.loci == ['url2']
    assert rest.profile == 'profile_text'


def test_read_contig_dynamic_base_cutoff_and_seqio_fallback(monkeypatch, tmp_path):
    fasta = tmp_path / 'reference.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            self.references = []
        def close(self):
            pass

    class FakeFastaFileIndexed(FakeFastaFile):
        def __init__(self, path):
            self.references = [b'abc_1']
        def close(self):
            pass

    def fake_pysam_fasta(path):
        if not hasattr(fake_pysam_fasta, 'called'):
            fake_pysam_fasta.called = True
            return FakeFastaFile(path)
        return FakeFastaFileIndexed(path)

    monkeypatch.setattr(methods.pysam, 'FastaFile', fake_pysam_fasta)
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: None)
    monkeypatch.setattr(methods.SeqIO, 'index', lambda path, fmt: (_ for _ in ()).throw(FileNotFoundError('missing')))
    monkeypatch.setattr(methods.SeqIO, 'to_dict', lambda recs: {'abc_1': SimpleNamespace(seq='ACGT')})
    monkeypatch.setattr(methods.SeqIO, 'parse', lambda path, fmt: [SimpleNamespace(id='abc_1', seq='ACGT')])

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'abc_1'
            self.pileups = []
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), [FakeColumn()]))
    monkeypatch.setattr(methods, 'characterise_read', lambda **kwargs: ({}, [30], {}))
    monkeypatch.setattr(methods, 'find_multibase_positions', lambda **kwargs: ({}, {}, 0, {}))

    multibase, report = methods.read_contig(
        contig_name='abc_1',
        bamfile_name=str(tmp_path / 'dummy.bam'),
        reference_fasta=str(fasta),
        allele_records=None,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=0,
        error_cutoff=1.0
    )
    assert 'abc_1' in multibase
    assert '_gene_stats' in multibase['abc_1']
    assert report == ''


def test_number_of_bases_above_threshold():
    counts = {'A': 3, 'C': 1, 'G': 2}
    assert methods.number_of_bases_above_threshold(
        high_quality_base_count=counts,
        base_count_cutoff=2
    ) == 2
    assert methods.number_of_bases_above_threshold(
        high_quality_base_count=counts,
        base_count_cutoff=2,
        base_fraction_cutoff=0.4
    ) == 1


def test_poisson_binomial_and_beta_binomial_helpers():
    pmf = methods._poisson_binomial_pmf_fft(p_list=[0.5, 0.5])
    assert pytest.approx(sum(pmf), rel=1e-9) == 1.0
    assert len(pmf) == 3
    assert methods.poisson_binomial_tail(k=1, p_list=[0.5, 0.5]) == pytest.approx(0.75)
    assert methods.poisson_binomial_tail(k=0, p_list=[]) == 1.0
    assert methods._beta_binomial_tail(k=0, n=1, a=1.0, b=1.0) == pytest.approx(1.0)
    a, b = methods._estimate_beta_params(p_list=[0.1, 0.5, 0.9])
    assert a is not None and b is not None
    assert methods._estimate_beta_params(p_list=[0.5, 0.5, 0.5]) == (None, None)


def test_statistical_and_multiple_testing_helpers():
    pvalue = methods.mann_whitney_u_p(x=[1, 2], y=[3, 4])
    assert isinstance(pvalue, float)
    assert 0.0 <= pvalue <= 1.0

    fisher_p = methods.fisher_two_sided_p(a=1, b=9, c=1, d=9)
    assert 0.0 <= fisher_p <= 1.0

    qvals = methods.benjamini_hochberg(pvals=[0.01, 0.04, 0.03])
    assert qvals == [0.01, 0.12, 0.045]

    combined = methods.combine_pvalues_fisher(pvals=[0.1, 0.2])
    assert 0.0 < combined < 1.0


def test_position_entry_passes_probabilistic_gating():
    stats = {
        'A': {'q_value': 0.01, 'strand_p': 0.5, 'pos_p': 0.5},
        'C': {'q_value': 0.2}
    }
    assert methods._position_entry_passes_probabilistic_gating(
        position_stats=stats,
        q_threshold=0.05
    ) is True


def test_find_multibase_positions_and_position_details():
    filtered_read_dict = {
        'congruent_ref': {'A': 2},
        'forward_SNV_reverse_ref': {'G': 3}
    }
    snv_dict, passing_snv_dict, total_coverage, position_stats = methods.find_multibase_positions(
        ref_base='A',
        filtered_read_dict=filtered_read_dict,
        base_cutoff=2,
        base_fraction_cutoff=0.1
    )
    assert total_coverage == 5
    assert snv_dict['total_forward_SNV'] == 3
    assert passing_snv_dict['forward']['G'] == 3
    assert passing_snv_dict['paired']['G'] == 3
    assert position_stats == {}

    row = methods.position_details(
        actual_position=1,
        passing_snv_dict={
            'congruent': {'A': 2},
            'paired': {'G': 3},
            'forward': {},
            'reverse': {}
        },
        contig_name='contig1',
        ref_base='A',
        total_coverage=5,
        base_cutoff=2,
        error_perc=0.1234,
        p_value=0.001,
        adj_p_value=0.002,
        strand_p=0.003,
        pos_p=0.004,
        mean_q=30.5,
        mean_mapq=60.5
    )
    assert row.startswith('contig1\t1\tA\t')
    assert '0.12' in row
    assert '1.000e-03' in row
    assert '60.50' in row


def test_find_multibase_positions_with_base_support_and_fraction():
    filtered_read_dict = {
        'congruent_ref': {'A': 2},
        'forward_SNV_reverse_ref': {'C': 3},
        'reverse_SNV_forward_ref': {'G': 1}
    }
    base_support = {
        'C': {'quals': [30, 30, 30], 'forward': 3, 'reverse': 0, 'positions': [1, 2, 3]},
        'G': {'quals': [30, 30], 'forward': 0, 'reverse': 2, 'positions': [1, 2]}
    }
    snv_dict, passing_snv_dict, total_coverage, position_stats = methods.find_multibase_positions(
        ref_base='A',
        filtered_read_dict=filtered_read_dict,
        base_cutoff=2,
        base_fraction_cutoff=0.1,
        base_support=base_support
    )
    assert total_coverage == 6
    assert passing_snv_dict['forward']['C'] == 3
    assert passing_snv_dict['reverse']['G'] == 0
    assert passing_snv_dict['paired']['C'] == 3
    assert 'C' in position_stats
    assert position_stats['C']['p_value'] is not None


def test_find_multibase_positions_without_fractional_cutoff_calls_count_cutoff():
    filtered_read_dict = {
        'congruent_ref': {'A': 2},
        'congruent_SNV': {'G': 3}
    }
    snv_dict, passing_snv_dict, total_coverage, position_stats = methods.find_multibase_positions(
        ref_base='A',
        filtered_read_dict=filtered_read_dict,
        base_cutoff=2,
        base_fraction_cutoff=0.0
    )

    assert total_coverage == 5
    assert snv_dict['total_congruent'] == 5
    assert snv_dict['total_congruent_SNV'] == 3
    assert snv_dict['total_forward'] == 2
    assert snv_dict['total_reverse'] == 2
    assert snv_dict['total_SNV'] == 3
    assert passing_snv_dict['congruent']['G'] == 3
    assert passing_snv_dict['paired']['G'] == 3
    assert position_stats == {}


def test_characterise_read_with_missing_fastq_record_uses_quality_zero():
    class FakeAlignment:
        def __init__(self, qname, seq, is_read1):
            self.qname = qname
            self.query_sequence = seq
            self.is_read1 = is_read1
            self.is_read2 = not is_read1
            self.mate_is_unmapped = False
            self.is_paired = True
            self.mapping_quality = 30
            self.query_alignment_end = len(seq)

    class FakePileupRead:
        def __init__(self, alignment, query_position=0):
            self.query_position = query_position
            self.alignment = alignment

    column = SimpleNamespace(
        pos=0,
        reference_name='gene1',
        pileups=[FakePileupRead(FakeAlignment('read1', 'C', True))]
    )

    original_state = methods._FASTQ_INDEX_STATE.copy()
    methods._FASTQ_INDEX_STATE.clear()
    try:
        filtered, qualities, support = methods.characterise_read(
            column=column,
            reference_sequence='A',
            fastq_records={},
            quality_cutoff=20,
            min_quality=15,
            fasta=False,
            nanopore=False
        )
    finally:
        methods._FASTQ_INDEX_STATE.clear()
        methods._FASTQ_INDEX_STATE.update(original_state)

    assert filtered['forward_quality_filtered'] == {'C': 1}
    assert qualities == []
    assert support == {}


def test_find_multibase_positions_ref_absent_returns_empty():
    filtered_read_dict = {
        'forward_SNV_reverse_ref': {'C': 3}
    }
    snv_dict, passing_snv_dict, total_coverage, position_stats = methods.find_multibase_positions(
        ref_base='A',
        filtered_read_dict=filtered_read_dict,
        base_cutoff=2,
        base_fraction_cutoff=0.0
    )

    assert total_coverage == 3
    assert passing_snv_dict == {}
    assert snv_dict['total'] == 3
    assert position_stats == {}


def test_find_multibase_positions_reverse_snv_fractional_cutoff():
    filtered_read_dict = {
        'congruent_ref': {'A': 1},
        'reverse_SNV_forward_ref': {'G': 3}
    }
    snv_dict, passing_snv_dict, total_coverage, position_stats = methods.find_multibase_positions(
        ref_base='A',
        filtered_read_dict=filtered_read_dict,
        base_cutoff=2,
        base_fraction_cutoff=0.5
    )

    assert total_coverage == 4
    assert passing_snv_dict['reverse']['G'] == 3
    assert passing_snv_dict['paired']['G'] == 3
    assert position_stats == {}


def test_find_contamination_skips_downsampling_for_ambiguous_genus(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    output_folder = tmp_path / 'out'
    output_folder.mkdir()
    reads = tmp_path / 'reads_R1.fastq.gz'
    with gzip.open(reads, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'ND')

    called = {'wrote': False}
    def fake_write_output(**kwargs):
        called['wrote'] = True
        assert kwargs['genus'] == 'ND'
    monkeypatch.setattr(methods, 'write_output', fake_write_output)

    methods.find_contamination(
        pair=[str(reads)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        downsample_depth=10
    )

    assert called['wrote'] is True


def test_find_contamination_uses_multiprocessing_pool_when_fastq_records_missing(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_database = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')
    Path(str(sample_database) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_downsample_reads(pair, sample_tmp_dir, sample_name, target_reads, log, xmx=None, seed=None):
        out1 = Path(sample_tmp_dir) / f'{sample_name}_downsampled_R1.fastq.gz'
        out2 = Path(sample_tmp_dir) / f'{sample_name}_downsampled_R2.fastq.gz'
        out1.write_bytes(b'@r1/1\nACGT\n+\nIIII\n')
        out2.write_bytes(b'@r1/2\nACGT\n+\nIIII\n')
        return [str(out1), str(out2)]

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return '', '', 'bbduk_bait'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return '', '', 'bbduk_trim'

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            outpath = cmd.split('out=')[1].split()[0]
            Path(outpath).write_bytes(b'')
        return '', ''

    def fake_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        Path(output).write_bytes(b'')

    def fake_index(path):
        Path(str(path) + '.bai').write_bytes(b'')

    class FakePool:
        def __init__(self, processes, initializer=None, initargs=None):
            pass
        def map(self, func, iterable, chunksize=1):
            return [({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n')]
        def close(self):
            pass
        def join(self):
            pass

    monkeypatch.setattr(methods, 'estimate_genome_size', lambda genus: 100)
    monkeypatch.setattr(methods, 'count_fastq_reads', lambda fastq_path: 5000)
    monkeypatch.setattr(methods, 'estimate_mean_read_length', lambda fastq_path: 100)
    monkeypatch.setattr(methods, 'downsample_reads', fake_downsample_reads)
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'sort', fake_sort)
    monkeypatch.setattr(methods.pysam, 'index', fake_index)
    monkeypatch.setattr(methods, '_ensure_fastq_index', lambda *args, **kwargs: None)
    def fake_load_fastq_records(*args, **kwargs):
        if kwargs.get('forward'):
            return {'read1/1': SimpleNamespace()}
        return {'read1/2': SimpleNamespace()}
    monkeypatch.setattr(methods, 'load_fastq_records', fake_load_fastq_records)
    monkeypatch.setattr(methods.multiprocessing, 'Pool', FakePool)
    monkeypatch.setattr(methods, '_read_contig_chunk_dispatch', lambda kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        threads=1,
        keep_files=False,
        min_matching_hashes=40,
        downsample_depth=1,
        subreplicates=2
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_characterise_read_same_snv_in_forward_and_reverse(tmp_path):
    def make_read(qname, seq, is_read1, query_position=0):
        alignment = SimpleNamespace(
            qname=qname,
            query_sequence=seq,
            query_alignment_end=len(seq),
            is_read1=is_read1,
            is_read2=not is_read1,
            is_paired=True,
            mate_is_unmapped=False,
            mapping_quality=30
        )
        return SimpleNamespace(
            query_position=query_position,
            alignment=alignment
        )

    forward = make_read('read1', 'C', True)
    reverse = make_read('read1', 'C', False)
    column = SimpleNamespace(
        pos=0,
        reference_name='gene1',
        pileups=[forward, reverse]
    )
    fastq_records = {
        'read1/1': SeqRecord(Seq('C'), id='read1/1', letter_annotations={'phred_quality': [30]}),
        'read1/2': SeqRecord(Seq('C'), id='read1/2', letter_annotations={'phred_quality': [30]})
    }
    filtered, qualities, base_support = methods.characterise_read(
        column=column,
        reference_sequence='AC',
        fastq_records=fastq_records,
        quality_cutoff=20,
        min_quality=15,
        fasta=False,
        nanopore=False
    )

    assert filtered['congruent_SNV'] == {'C': 2}
    assert filtered['congruent_ref'] == {}
    assert qualities == [30, 30]
    assert base_support['C']['forward'] == 1
    assert base_support['C']['reverse'] == 1


def test_characterise_read_forward_snv_reverse_ref(tmp_path):
    def make_read(qname, seq, is_read1, query_position=0):
        alignment = SimpleNamespace(
            qname=qname,
            query_sequence=seq,
            query_alignment_end=len(seq),
            is_read1=is_read1,
            is_read2=not is_read1,
            is_paired=True,
            mate_is_unmapped=False,
            mapping_quality=30
        )
        return SimpleNamespace(
            query_position=query_position,
            alignment=alignment
        )

    forward = make_read('read2', 'C', True)
    reverse = make_read('read2', 'A', False)
    column = SimpleNamespace(
        pos=0,
        reference_name='gene1',
        pileups=[forward, reverse]
    )
    fastq_records = {
        'read2/1': SeqRecord(Seq('C'), id='read2/1', letter_annotations={'phred_quality': [30]}),
        'read2/2': SeqRecord(Seq('A'), id='read2/2', letter_annotations={'phred_quality': [30]})
    }
    filtered, qualities, base_support = methods.characterise_read(
        column=column,
        reference_sequence='AC',
        fastq_records=fastq_records,
        quality_cutoff=20,
        min_quality=15,
        fasta=False,
        nanopore=False
    )

    assert filtered['forward_SNV_reverse_ref'] == {'C': 1, 'A': 1}
    assert qualities == [30, 30]
    assert base_support['C']['forward'] == 1
    assert base_support['A']['reverse'] == 1


def test_characterise_read_different_snvs_both_pass_quality(tmp_path):
    def make_read(qname, seq, is_read1):
        alignment = SimpleNamespace(
            qname=qname,
            query_sequence=seq,
            query_alignment_end=len(seq),
            is_read1=is_read1,
            is_read2=not is_read1,
            is_paired=True,
            mate_is_unmapped=False,
            mapping_quality=30
        )
        return SimpleNamespace(
            query_position=0,
            alignment=alignment
        )

    forward = make_read('read3', 'C', True)
    reverse = make_read('read3', 'G', False)
    column = SimpleNamespace(
        pos=0,
        reference_name='gene1',
        pileups=[forward, reverse]
    )
    fastq_records = {
        'read3/1': SeqRecord(Seq('C'), id='read3/1', letter_annotations={'phred_quality': [30]}),
        'read3/2': SeqRecord(Seq('G'), id='read3/2', letter_annotations={'phred_quality': [30]})
    }
    filtered, qualities, base_support = methods.characterise_read(
        column=column,
        reference_sequence='AC',
        fastq_records=fastq_records,
        quality_cutoff=20,
        min_quality=15,
        fasta=False,
        nanopore=False
    )

    assert filtered['forward_SNV_reverse_SNV1'] == {'C': 1}
    assert filtered['reverse_SNV_forward_SNV1'] == {'G': 1}
    assert qualities == [30, 30]


def test_characterise_read_different_snvs_only_forward_quality(tmp_path):
    def make_read(qname, seq, is_read1):
        alignment = SimpleNamespace(
            qname=qname,
            query_sequence=seq,
            query_alignment_end=len(seq),
            is_read1=is_read1,
            is_read2=not is_read1,
            is_paired=True,
            mate_is_unmapped=False,
            mapping_quality=30
        )
        return SimpleNamespace(
            query_position=0,
            alignment=alignment
        )

    forward = make_read('read4', 'C', True)
    reverse = make_read('read4', 'G', False)
    column = SimpleNamespace(
        pos=0,
        reference_name='gene1',
        pileups=[forward, reverse]
    )
    fastq_records = {
        'read4/1': SeqRecord(Seq('C'), id='read4/1', letter_annotations={'phred_quality': [30]}),
        'read4/2': SeqRecord(Seq('G'), id='read4/2', letter_annotations={'phred_quality': [10]})
    }
    filtered, qualities, base_support = methods.characterise_read(
        column=column,
        reference_sequence='AC',
        fastq_records=fastq_records,
        quality_cutoff=20,
        min_quality=15,
        fasta=False,
        nanopore=False
    )

    assert filtered['forward_SNV_reverse_UM_QF'] == {'C': 1}
    assert filtered['reverse_quality_filtered'] == {'G': 1}
    assert qualities == [30]


def test_characterise_read_quality_filtered_unpaired_reads():
    class FakeAlignment:
        def __init__(self, qname, seq, is_read1):
            self.qname = qname
            self.query_sequence = seq
            self.mapping_quality = 5
            self.is_read1 = is_read1
            self.is_read2 = not is_read1
            self.is_paired = False
            self.mate_is_unmapped = True
            self.query_alignment_end = len(seq)

    class FakePileupRead:
        def __init__(self, alignment, query_position=0):
            self.query_position = query_position
            self.alignment = alignment

    forward = FakeAlignment('read1', 'C', True)
    column = SimpleNamespace(
        pos=0,
        reference_name='gene1',
        pileups=[FakePileupRead(forward)]
    )
    fastq_records = {
        'read1': SeqRecord(Seq('C'), id='read1', letter_annotations={'phred_quality': [10]})
    }

    filtered, qualities, base_support = methods.characterise_read(
        column=column,
        reference_sequence='A',
        fastq_records=fastq_records,
        quality_cutoff=20,
        min_quality=15,
        fasta=False,
        nanopore=False
    )

    assert filtered['forward_quality_filtered'] == {'C': 1}
    assert qualities == []
    assert base_support == {}


def test_characterise_read_different_snvs_only_reverse_quality(tmp_path):
    def make_read(qname, seq, is_read1):
        alignment = SimpleNamespace(
            qname=qname,
            query_sequence=seq,
            query_alignment_end=len(seq),
            is_read1=is_read1,
            is_read2=not is_read1,
            is_paired=True,
            mate_is_unmapped=False,
            mapping_quality=30
        )
        return SimpleNamespace(
            query_position=0,
            alignment=alignment
        )

    forward = make_read('read5', 'C', True)
    reverse = make_read('read5', 'G', False)
    column = SimpleNamespace(
        pos=0,
        reference_name='gene1',
        pileups=[forward, reverse]
    )
    fastq_records = {
        'read5/1': SeqRecord(Seq('C'), id='read5/1', letter_annotations={'phred_quality': [10]}),
        'read5/2': SeqRecord(Seq('G'), id='read5/2', letter_annotations={'phred_quality': [30]})
    }
    filtered, qualities, base_support = methods.characterise_read(
        column=column,
        reference_sequence='AC',
        fastq_records=fastq_records,
        quality_cutoff=20,
        min_quality=15,
        fasta=False,
        nanopore=False
    )

    assert filtered['reverse_SNV_forward_UM_QF'] == {'G': 1}
    assert filtered['forward_quality_filtered'] == {'C': 1}
    assert qualities == [30]


def test_characterise_read_different_snvs_neither_quality(tmp_path):
    def make_read(qname, seq, is_read1):
        alignment = SimpleNamespace(
            qname=qname,
            query_sequence=seq,
            query_alignment_end=len(seq),
            is_read1=is_read1,
            is_read2=not is_read1,
            is_paired=True,
            mate_is_unmapped=False,
            mapping_quality=30
        )
        return SimpleNamespace(
            query_position=0,
            alignment=alignment
        )

    forward = make_read('read6', 'C', True)
    reverse = make_read('read6', 'G', False)
    column = SimpleNamespace(
        pos=0,
        reference_name='gene1',
        pileups=[forward, reverse]
    )
    fastq_records = {
        'read6/1': SeqRecord(Seq('C'), id='read6/1', letter_annotations={'phred_quality': [10]}),
        'read6/2': SeqRecord(Seq('G'), id='read6/2', letter_annotations={'phred_quality': [10]})
    }
    filtered, qualities, base_support = methods.characterise_read(
        column=column,
        reference_sequence='AC',
        fastq_records=fastq_records,
        quality_cutoff=20,
        min_quality=15,
        fasta=False,
        nanopore=False
    )

    assert filtered['forward_quality_filtered'] == {'C': 1}
    assert filtered['reverse_quality_filtered'] == {'G': 1}
    assert qualities == []


def test_find_contamination_returns_early_when_database_missing(tmp_path, monkeypatch):
    input_file = tmp_path / 'reads_R1.fastq.gz'
    input_file.write_bytes(b'@r1/1\nACGT\n+\nIIII\n')
    output_folder = tmp_path / 'out'
    output_folder.mkdir()
    db_folder = tmp_path / 'db'
    db_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda **kwargs: 'Escherichia')
    def fake_write_output(**kwargs):
        with open(kwargs['output_report'], 'w', encoding='utf-8') as f:
            f.write('Sample\tGenus\tNumContamSNVs\tContamStatus\tBasesExamined\tDatabaseDownloadDate\tScore\n')
            f.write('sample\tEscherichia\t0\tFalse\t0\tND\t0\n')
    monkeypatch.setattr(methods, 'write_output', fake_write_output)

    methods.find_contamination(
        pair=[str(input_file)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path / 'tmp'),
        data_type='Illumina'
    )

    assert (output_folder / 'confindr_report.tsv').exists()


def test_find_contamination_creates_rmlst_database_fallback_and_writes_report(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    (db_folder / 'gene_allele.txt').write_text(
        'Escherichia:abc_1,\n',
        encoding='utf-8'
    )
    (db_folder / 'rMLST_combined.fasta').write_text(
        '>abc_1\nACGT\n',
        encoding='utf-8'
    )

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    forward = tmp_path / 'reads_R1.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')

    def fake_find_cross_contamination(**kwargs):
        return 'Escherichia'

    def fake_find_genus_specific_allele_list(*args, **kwargs):
        return ['abc_1']

    def fake_setup_allelespecific_database(**kwargs):
        with open(kwargs['fasta_file'], 'w', encoding='utf-8') as f:
            f.write('>abc_1\nACGT\n')

    def fake_bbduk_bait(**kwargs):
        out_file = kwargs.get('forward_out')
        with gzip.open(out_file, 'wt', encoding='utf-8') as f:
            f.write('@r1/1\nACGT\n+\nIIII\n')
        return ('', '', '')

    def fake_bbduk_trim(**kwargs):
        out_file = kwargs.get('forward_out')
        with gzip.open(out_file, 'wt', encoding='utf-8') as f:
            f.write('@r1/1\nACGT\n+\nIIII\n')
        return ('', '', '')

    monkeypatch.setattr(methods, 'find_cross_contamination', fake_find_cross_contamination)
    monkeypatch.setattr(methods, 'find_genus_specific_allele_list', fake_find_genus_specific_allele_list)
    monkeypatch.setattr(methods, 'setup_allelespecific_database', fake_setup_allelespecific_database)
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: sample_database.replace('.fasta', '_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'find_total_sequence_length', lambda fasta_file: 100)
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}}}, ''))
    monkeypatch.setattr(methods, 'write_to_logfile', lambda **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', lambda cmd: ('', ''))
    monkeypatch.setattr(methods.pysam, 'faidx', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods.pysam, 'index', lambda *args, **kwargs: None)

    methods.find_contamination(
        pair=[str(forward)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path / 'tmp'),
        data_type='Illumina'
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_R1_contamination.tsv').exists()


def test_find_multibase_positions_returns_empty_when_ref_absent():
    filtered_read_dict = {
        'forward_SNV_reverse_ref': {'C': 3},
        'reverse_SNV_forward_ref': {'G': 2}
    }
    snv_dict, passing_snv_dict, total_coverage, position_stats = methods.find_multibase_positions(
        ref_base='A',
        filtered_read_dict=filtered_read_dict,
        base_cutoff=2,
        base_fraction_cutoff=0.2
    )
    assert total_coverage == 5
    assert snv_dict['total'] == 5
    assert passing_snv_dict == {}
    assert position_stats == {}


def test_position_entry_passes_probabilistic_gating_false():
    stats = {
        'A': {'q_value': 0.2, 'strand_p': 0.5, 'pos_p': 0.5},
        'C': {'q_value': 0.01, 'strand_p': 0.001, 'pos_p': 0.5}
    }
    assert methods._position_entry_passes_probabilistic_gating(
        position_stats=stats,
        q_threshold=0.05,
        strand_p_threshold=0.01,
        pos_p_threshold=0.01
    ) is False


def test_benjamini_hochberg_empty_and_unsorted():
    assert methods.benjamini_hochberg(pvals=[]) == []
    qvals = methods.benjamini_hochberg(pvals=[0.2, 0.1])
    assert qvals[0] >= qvals[1]


def test_combine_pvalues_fisher_empty_and_single():
    assert methods.combine_pvalues_fisher(pvals=[]) is None
    assert 0.0 < methods.combine_pvalues_fisher(pvals=[0.5]) < 1.0


def test_poisson_binomial_tail_large_n():
    p_list = [0.5] * 1001
    tail = methods.poisson_binomial_tail(k=500, p_list=p_list)
    assert 0.0 <= tail <= 1.0


def test_estimate_beta_params_edge_cases():
    assert methods._estimate_beta_params(p_list=[]) == (None, None)
    a, b = methods._estimate_beta_params(p_list=[0.1, 0.5, 0.9])
    assert a is not None and b is not None


def test_format_seconds_invalid_inputs_returns_na():
    assert methods._format_seconds(s='not a number') == 'N/A'
    assert methods._format_seconds(s=float('nan')) == 'N/A'


def test_format_seconds_unroundable_input_returns_na():
    class BadNumber:
        def __round__(self):
            raise TypeError('cannot round')

    assert methods._format_seconds(s=BadNumber()) == 'N/A'


def test_determine_cutoff_dynamic_base_cutoff():
    k, expected_positions, error_perc = methods.determine_cutoff(
        qualities=[30] * 10,
        reference_sequence='A' * 100,
        base_cutoff=0,
        error_cutoff=1.0,
        max_expected_positions=0.01
    )
    assert k >= methods.MIN_DYNAMIC_CUTOFF
    assert expected_positions >= 0.0
    assert 0.0 <= error_perc <= 100.0


def test_determine_cutoff_no_qualities_uses_minimum_cutoff():
    assert methods.determine_cutoff(
        qualities=[],
        reference_sequence='A' * 100,
        base_cutoff=1
    ) == (1, 0.0, 0.0)


def test_poisson_binomial_pmf_fft_edge_cases():
    result = methods._poisson_binomial_pmf_fft(p_list=[])
    assert len(result) == 1
    assert result[0] == 1.0

    result = methods._poisson_binomial_pmf_fft(p_list=[0.5])
    assert result.shape == (2,)
    assert pytest.approx(result.sum(), rel=1e-7) == 1.0
    assert pytest.approx(result[0] + result[1], rel=1e-7) == 1.0


def test_beta_binomial_tail_and_poisson_binomial_tail_degenerate():
    assert methods._beta_binomial_tail(k=0, n=0, a=1.0, b=1.0) == 1.0
    assert methods._beta_binomial_tail(k=1, n=0, a=1.0, b=1.0) == 0.0
    assert methods.poisson_binomial_tail(k=0, p_list=[]) == 1.0
    assert methods.poisson_binomial_tail(k=1, p_list=[]) == 0.0
    assert 0.0 <= methods.poisson_binomial_tail(k=1, p_list=[0.2]) <= 1.0


def test_mann_whitney_u_p_empty_inputs_returns_none():
    assert methods.mann_whitney_u_p(x=[], y=[1, 2, 3]) is None
    assert methods.mann_whitney_u_p(x=[1, 2, 3], y=[]) is None


def test_fisher_two_sided_p_invalid_inputs_returns_one():
    assert methods.fisher_two_sided_p(a=-1, b=1, c=1, d=1) == 1.0
    assert methods.fisher_two_sided_p(a=1, b=-1, c=1, d=1) == 1.0


def test_find_contamination_unpaired_nanopore_path(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_db = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_db.write_text('>abc_1\nACGT\n', encoding='utf-8')
    Path(str(sample_db) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    forward = tmp_path / 'reads.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    def fake_run_cmd(cmd):
        if '-o ' in cmd:
            out_path = cmd.split('-o ')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return 'out', 'err'

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(path + '.bai').write_bytes(b''))
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', lambda *args, **kwargs: SimpleNamespace(references=[b'abc_1'], close=lambda: None))
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        threads=1,
        keep_files=False,
        data_type='Nanopore',
        min_matching_hashes=40
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_find_contamination_skips_when_no_present_alleles(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_db = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_db.write_text('>abc_1\nACGT\n', encoding='utf-8')
    Path(str(sample_db) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    forward = tmp_path / 'reads.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return 'out', 'err'

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(path + '.bai').write_bytes(b''))
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({}, ''))
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', lambda *args, **kwargs: SimpleNamespace(references=[b'other'], close=lambda: None))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        keep_files=False,
        min_matching_hashes=40
    )

    assert (output_folder / 'confindr_report.tsv').exists()


def test_find_contamination_aggregates_subreplicate_reports(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_db = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_db.write_text('>abc_1\nACGT\n', encoding='utf-8')
    Path(str(sample_db) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    forward = tmp_path / 'reads.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_bbduk_bait(*args, **kwargs):
        out_file = kwargs.get('forward_out')
        if out_file:
            Path(out_file).write_bytes(b'')
        reverse_out = kwargs.get('reverse_out')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return '', '', ''

    def fake_bbduk_trim(*args, **kwargs):
        out_file = kwargs.get('forward_out')
        if out_file:
            Path(out_file).write_bytes(b'')
        reverse_out = kwargs.get('reverse_out')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return '', '', ''

    def fake_run_cmd(cmd):
        if 'out=' in cmd and cmd.endswith('.bam'):
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'find_total_sequence_length', lambda fasta_file: 100)
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))
    monkeypatch.setattr(methods, 'write_to_logfile', lambda **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: None)
    def fake_sort(*args, **kwargs):
        if len(args) > 1:
            Path(args[1]).write_bytes(b'')
    monkeypatch.setattr(methods.pysam, 'sort', fake_sort)
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(path + '.bai').write_bytes(b''))
    monkeypatch.setattr(methods.SeqIO, 'parse', lambda path, fmt: [SeqRecord(Seq('ACGT'), id='abc_1')])

    methods.find_contamination(
        pair=[str(forward)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        threads=1,
        keep_files=False,
        min_matching_hashes=40,
        subreplicates=2
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_build_contig_chunks_missing_contig_length():
    chunks = methods._build_contig_chunks(
        contig_lengths={},
        contigs=['x'],
        threads=1,
        multiplier=1,
        max_chunk_bases=100
    )
    assert chunks == [[{'contig': 'x', 'start': None, 'end': None, 'length': 1}]]


def test_read_contig_chunk_dispatch_error(monkeypatch):
    def fake_read_contig(**kwargs):
        raise ValueError('boom')

    monkeypatch.setattr(methods, 'read_contig', fake_read_contig)
    with pytest.raises(ValueError):
        methods._read_contig_chunk_dispatch({
            'contig_chunk': [{'contig': 'a', 'start': 0, 'end': 1}],
            'bamfile_name': 'dummy.bam',
            'reference_fasta': 'dummy.fasta',
            'fastq_records': {},
            'quality_cutoff': 20,
            'base_cutoff': 1,
            'base_fraction_cutoff': 0.1,
            'fasta': False,
            'error_cutoff': 1.0,
            'nanopore': False,
            'max_expected_positions': 0.001
        })


def test_read_contig_dynamic_cutoff_and_fasta_mode(monkeypatch):
    allele_records = {'abc_1': SimpleNamespace(seq='ACGT')}

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            return [b'abc_1']

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'abc_1'
            self.pileups = []

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), [FakeColumn()]))
    monkeypatch.setattr(methods, 'characterise_read', lambda **kwargs: ({'congruent_ref': {'A': 2}}, [30, 30], {}))
    monkeypatch.setattr(methods, 'determine_cutoff', lambda **kwargs: (3, 0.0, 0.0))
    monkeypatch.setattr(methods, 'find_multibase_positions', lambda **kwargs: ({'total': 3}, {'forward': {'C': 2}, 'congruent': {}, 'reverse': {}, 'paired': {}}, 3, {'C': {'p_value': 0.01, 'strand_p': 0.5, 'pos_p': 0.5}}))
    monkeypatch.setattr(methods, 'benjamini_hochberg', lambda pvals: [0.01])
    monkeypatch.setattr(methods, '_position_entry_passes_probabilistic_gating', lambda position_stats: True)
    monkeypatch.setattr(methods, 'combine_pvalues_fisher', lambda pvals: 0.05)
    monkeypatch.setattr(methods, 'position_details', lambda **kwargs: 'abc_1\t1\tA\n')

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta='dummy.fasta',
        allele_records=allele_records,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=0,
        fasta=True
    )

    assert 'abc_1' in result
    assert report_text == 'abc_1\t1\tA\n'


def test_bbtools_bbduk_trim_raises_file_not_found_when_binary_missing(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    forward_out = tmp_path / 'trimmed_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_check_output(cmd):
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(bbtools.subprocess, 'check_output', fake_check_output)
    with pytest.raises(FileNotFoundError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(forward_out),
            reverse_in=str(reverse),
            reverse_out=str(tmp_path / 'trimmed_R2.fastq')
        )


def test_bbtools_reformat_reads_requires_reverse_out_if_reverse_in_provided(monkeypatch, tmp_path):
    assert methods.determine_cutoff(
        qualities=[],
        reference_sequence='A' * 10,
        base_cutoff=1
    ) == (1, 0.0, 0.0)
    k, expected_positions, error_perc = methods.determine_cutoff(
        qualities=[30] * 10,
        reference_sequence='A' * 100,
        base_cutoff=1,
        error_cutoff=1.0,
        max_expected_positions=0.01
    )
    assert k >= 1
    assert expected_positions >= 0.0
    assert 0.0 <= error_perc <= 100.0


def test_write_output_writes_summary_line(tmp_path):
    report = tmp_path / 'report.tsv'
    methods.write_output(
        output_report=str(report),
        sample_name='sample1',
        multi_positions=5,
        genus='Escherichia',
        total_gene_length=1000,
        database_download_date='2026-07-08',
        snp_cutoff=3,
        pysam_pass=True,
        sample_score=2.5,
        use_probabilistic=False,
        score_threshold=None
    )
    content = report.read_text(encoding='utf-8').splitlines()
    assert content[0].startswith('Sample\tGenus\tNumContamSNVs\tContamStatus')
    assert content[1].startswith('sample1\tEscherichia\t5\tTrue\t1000\t2026-07-08\t2.500')


def test_write_output_subreplicates_summary(tmp_path):
    report = tmp_path / 'report.tsv'
    methods.write_output(
        output_report=str(report),
        sample_name='sample1',
        multi_positions=0,
        genus='Escherichia',
        total_gene_length=1000,
        database_download_date='2026-07-08',
        subreplicate_counts=[1, 2, 3],
        pysam_pass=True,
        sample_score=2.0,
        use_probabilistic=True,
        score_threshold=2.0
    )
    content = report.read_text(encoding='utf-8').splitlines()
    assert 'MeanContamSNVs' in content[0]
    fields = content[1].split('\t')
    assert fields[0] == 'sample1'
    assert fields[5] == 'True'
    assert fields[-1] == '2.000'


def test_write_output_pysam_fails(tmp_path):
    report = tmp_path / 'report.tsv'
    methods.write_output(
        output_report=str(report),
        sample_name='sample2',
        multi_positions=0,
        genus='Escherichia',
        total_gene_length=500,
        database_download_date='2026-07-08',
        pysam_pass=False,
        sample_score=None,
        use_probabilistic=False,
        score_threshold=None
    )
    content = report.read_text(encoding='utf-8').splitlines()
    assert content[1].split('\t')[3] == 'Pysam SamtoolsError'
    assert content[1].split('\t')[2] == 'ND'


def test_get_version():
    assert methods.get_version().startswith('ConFindr ')


def test_check_for_databases_and_download_optional_warns(monkeypatch, tmp_path):
    for name in [
        'Escherichia_db_cgderived.fasta',
        'Listeria_db_cgderived.fasta',
        'Salmonella_db_cgderived.fasta',
        'refseq.msh'
    ]:
        (tmp_path / name).write_text('x', encoding='utf-8')

    monkeypatch.setattr(methods, 'download_mash_sketch', lambda output_folder: (_ for _ in ()).throw(AssertionError('should not be called')))
    monkeypatch.setattr(methods, 'download_cgmlst_derived_data', lambda output_folder: (_ for _ in ()).throw(AssertionError('should not be called')))

    methods.check_for_databases_and_download(database_location=str(tmp_path))


def test_find_contamination_exits_on_missing_cgmlst_db(monkeypatch, tmp_path):
    output_folder = tmp_path / 'out'
    output_folder.mkdir()
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    reads = tmp_path / 'reads_R1.fastq.gz'
    reads.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    with pytest.raises(SystemExit) as excinfo:
        methods.find_contamination(
            pair=[str(reads)],
            output_folder=str(output_folder),
            databases_folder=str(db_folder),
            base_cutoff=3,
            xmx='4g',
            cgmlst_db=str(tmp_path / 'missing.fasta'),
            threads=1
        )
    assert excinfo.value.code == 1


def test_check_for_databases_and_download_invokes_downloads(monkeypatch, tmp_path):
    called = {'sketch': False, 'cgmlst': False}

    def fake_download_mash_sketch(output_folder):
        called['sketch'] = True
        (Path(output_folder) / 'refseq.msh').write_text('x', encoding='utf-8')

    def fake_download_cgmlst_derived_data(output_folder):
        called['cgmlst'] = True
        (Path(output_folder) / 'gene_allele.txt').write_text('x', encoding='utf-8')
        (Path(output_folder) / 'rMLST_combined.fasta').write_text('>a\nA\n', encoding='utf-8')

    monkeypatch.setattr(methods, 'download_mash_sketch', fake_download_mash_sketch)
    monkeypatch.setattr(methods, 'download_cgmlst_derived_data', fake_download_cgmlst_derived_data)

    methods.check_for_databases_and_download(database_location=str(tmp_path))
    assert called['sketch']
    assert called['cgmlst']


def test_check_for_databases_and_download_skips_when_present(monkeypatch, tmp_path):
    for name in [
        'Escherichia_db_cgderived.fasta',
        'Listeria_db_cgderived.fasta',
        'Salmonella_db_cgderived.fasta',
        'refseq.msh'
    ]:
        (tmp_path / name).write_text('x', encoding='utf-8')

    monkeypatch.setattr(methods, 'download_mash_sketch', lambda output_folder: (_ for _ in ()).throw(AssertionError('should not be called')))
    monkeypatch.setattr(methods, 'download_cgmlst_derived_data', lambda output_folder: (_ for _ in ()).throw(AssertionError('should not be called')))

    methods.check_for_databases_and_download(database_location=str(tmp_path))


def test_find_contamination_writes_output_when_no_database(monkeypatch, tmp_path):
    pair = [str(tmp_path / 'reads_R1.fastq.gz'), str(tmp_path / 'reads_R2.fastq.gz')]
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'ND')
    called = {'wrote': False}

    def fake_write_output(**kwargs):
        called['wrote'] = True
        assert kwargs['genus'] == 'ND'
        assert kwargs['sample_name'] == 'reads_' or kwargs['sample_name'].startswith('reads')

    monkeypatch.setattr(methods, 'write_output', fake_write_output)
    result = methods.find_contamination(
        pair=pair,
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False
    )
    assert result is None
    assert called['wrote']


def test_find_contamination_generates_reports_single_run(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_database = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    def fake_pysam_faidx(path):
        Path(path + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    def fake_pysam_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        Path(output).write_bytes(b'')

    def fake_pysam_index(path):
        Path(path + '.bai').write_bytes(b'')

    class FakeAlignmentFile:
        def __init__(self, *args, **kwargs):
            self.references = [b'abc_1']
        def close(self):
            pass

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    class FakePool:
        def __init__(self, processes):
            pass
        def map(self, func, iterable, chunksize=1):
            return [({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n')]
        def close(self):
            pass
        def join(self):
            pass

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', fake_pysam_faidx)
    monkeypatch.setattr(methods.pysam, 'sort', fake_pysam_sort)
    monkeypatch.setattr(methods.pysam, 'index', fake_pysam_index)
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))
    monkeypatch.setattr(methods, 'ThreadPool', FakePool)

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_find_contamination_downsampled_subreplicates_aggregates_results(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_database = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    def fake_find_cross_contamination(*args, **kwargs):
        return 'Escherichia'

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return '', '', 'bbduk_bait'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return '', '', 'bbduk_trim'

    def fake_downsample_reads(pair, sample_tmp_dir, sample_name, target_reads, log, xmx=None, seed=None):
        out1 = Path(sample_tmp_dir) / f'{sample_name}_downsampled_R1.fastq.gz'
        out2 = Path(sample_tmp_dir) / f'{sample_name}_downsampled_R2.fastq.gz'
        out1.write_bytes(b'')
        out2.write_bytes(b'')
        return [str(out1), str(out2)]

    def fake_run_cmd(cmd):
        if '-o ' in cmd:
            outpath = cmd.split('-o ')[1].split()[0]
            Path(outpath + '.res').write_text('', encoding='utf-8')
            return '', ''
        if 'bbmap.sh' in cmd or 'minimap2' in cmd:
            if 'out=' in cmd:
                out_path = cmd.split('out=')[1].split()[0]
                Path(out_path).write_bytes(b'')
            return '', ''
        return '', ''

    def fake_faidx(path):
        Path(str(path) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    def fake_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        Path(output).write_bytes(b'')

    def fake_index(path):
        Path(str(path) + '.bai').write_bytes(b'')

    class FakeAlignmentFile:
        def __init__(self, *args, **kwargs):
            self.references = [b'abc_1']
        def close(self):
            pass

    monkeypatch.setattr(methods, 'find_cross_contamination', fake_find_cross_contamination)
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'downsample_reads', fake_downsample_reads)
    monkeypatch.setattr(methods, 'count_fastq_reads', lambda fastq_path: 500000)
    monkeypatch.setattr(methods, 'estimate_mean_read_length', lambda fastq_path: 100)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'find_total_sequence_length', lambda fasta_file: 4)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', fake_faidx)
    monkeypatch.setattr(methods.pysam, 'sort', fake_sort)
    monkeypatch.setattr(methods.pysam, 'index', fake_index)
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    monkeypatch.setattr(methods, '_read_contig_chunk_dispatch', lambda kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False,
        use_prob_scoring=True,
        score_threshold=1.0,
        downsample_depth=10,
        subreplicates=2,
        subreplicate_seed=1,
        subreplicate_consensus=0.5
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_find_contamination_nonreplicate_with_fake_pipeline(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_database = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    def fake_find_cross_contamination(*args, **kwargs):
        return 'Escherichia'

    def fake_bbduk_bait(*args, **kwargs):
        if kwargs.get('forward_out'):
            Path(kwargs['forward_out']).write_bytes(b'')
        if kwargs.get('reverse_out'):
            Path(kwargs['reverse_out']).write_bytes(b'')
        return '', '', 'bbduk_bait'

    def fake_bbduk_trim(*args, **kwargs):
        if kwargs.get('forward_out'):
            Path(kwargs['forward_out']).write_bytes(b'')
        if kwargs.get('reverse_out'):
            Path(kwargs['reverse_out']).write_bytes(b'')
        return '', '', 'bbduk_trim'

    def fake_run_cmd(cmd):
        if 'kma -i' in cmd or 'kma -ipe' in cmd:
            out = ''
            # create the .res output file if not present
            if '-o ' in cmd:
                outpath = cmd.split('-o ')[1].split()[0]
                Path(outpath + '.res').write_text('', encoding='utf-8')
            return out, ''
        if 'bbmap.sh' in cmd or 'minimap2' in cmd:
            if 'out=' in cmd:
                out_path = cmd.split('out=')[1].split()[0]
                Path(out_path).write_bytes(b'')
            return '', ''
        return '', ''

    def fake_faidx(path):
        Path(str(path) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    def fake_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        Path(output).write_bytes(b'')

    def fake_index(path):
        Path(str(path) + '.bai').write_bytes(b'')

    class FakeAlignmentFile:
        def __init__(self, *args, **kwargs):
            self.references = [b'abc_1']
        def close(self):
            pass

    monkeypatch.setattr(methods, 'find_cross_contamination', fake_find_cross_contamination)
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'find_total_sequence_length', lambda fasta_file: 4)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', fake_faidx)
    monkeypatch.setattr(methods.pysam, 'sort', fake_sort)
    monkeypatch.setattr(methods.pysam, 'index', fake_index)
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}, '_gene_stats': {'combined_p': 0.05, 'gene_score': 1.0, 'num_sig_positions': 1}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False,
        use_prob_scoring=True,
        score_threshold=1.0
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_read_contig_dispatch_missing_contig_name():
    with pytest.raises(TypeError):
        methods._read_contig_dispatch(kwargs={})


def test_find_cross_contamination(monkeypatch, tmp_path):
    def fake_screen(*args, **kwargs):
        output_file = kwargs.get('output_file', 'screen.tab')
        (tmp_path / output_file).write_text(
            '99.9 40/100 1 1e-10 a/b/x/Shigella/c/d\n',
            encoding='utf-8'
        )
        return '', '', 'mash screen'

    monkeypatch.setattr(methods.mash, 'screen', fake_screen)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda logfile, out, err, cmd: None)

    result = methods.find_cross_contamination(
        databases=str(tmp_path),
        reads=str(tmp_path / 'reads.fastq'),
        sample_name='sample',
        tmpdir=str(tmp_path),
        log=str(tmp_path / 'log.txt'),
        threads=1,
        min_matching_hashes=40
    )
    assert result == 'Escherichia'


def test_find_cross_contamination_paired_reads_and_shigella_conversion(monkeypatch, tmp_path):
    def fake_screen(*args, **kwargs):
        output_file = kwargs.get('output_file', 'sample_screen.tab')
        (tmp_path / output_file).write_text(
            '99.9 40/100 1 1e-10 a/b/x/Shigella/c/d\n',
            encoding='utf-8'
        )
        return '', '', 'mash screen'

    monkeypatch.setattr(methods.mash, 'screen', fake_screen)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda logfile, out, err, cmd: None)

    result = methods.find_cross_contamination(
        databases=str(tmp_path),
        reads=[str(tmp_path / 'reads_R1.fastq'), str(tmp_path / 'reads_R2.fastq')],
        sample_name='sample',
        tmpdir=str(tmp_path),
        log=str(tmp_path / 'log.txt'),
        threads=1,
        min_matching_hashes=40
    )
    assert result == 'Escherichia'


def test_find_cross_contamination_uses_existing_screen_file(monkeypatch, tmp_path):
    screen_file = tmp_path / 'sample_screen.tab'
    screen_file.write_text(
        '99.9 40/100 1 1e-10 a/b/x/Escherichia/c/d\n',
        encoding='utf-8'
    )

    monkeypatch.setattr(methods.mash, 'screen', lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError('should not be called')))
    monkeypatch.setattr(methods.mash, 'read_mash_screen', lambda path: [SimpleNamespace(query_id='a/b/x/Escherichia/c/d', shared_hashes='40/100')])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda logfile, out, err, cmd: None)

    result = methods.find_cross_contamination(
        databases=str(tmp_path),
        reads=str(tmp_path / 'reads.fastq'),
        sample_name='sample',
        tmpdir=str(tmp_path),
        log=str(tmp_path / 'log.txt'),
        threads=1,
        min_matching_hashes=40
    )
    assert result == 'Escherichia'


def test_parse_bam_uses_bytes_contig_name(monkeypatch):
    called = {}

    class FakeAlignmentFile:
        def __init__(self, path, mode):
            called['path'] = path
            called['mode'] = mode
        def pileup(self, contig_name, *args, **kwargs):
            called['contig_name'] = contig_name
            called['kwargs'] = kwargs
            return []

    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    fake_fasta = object()

    bamfile, pileup = methods.parse_bam(
        bamfile_name='sample.bam',
        contig_name=b'abc_1',
        pysam_fasta=fake_fasta
    )

    assert isinstance(bamfile, FakeAlignmentFile)
    assert called['contig_name'] == 'abc_1'
    assert called['kwargs']['fastafile'] is fake_fasta
    assert pileup == []


def test_parse_bam_with_start_end(monkeypatch):
    called = {}

    class FakeAlignmentFile:
        def __init__(self, path, mode):
            pass
        def pileup(self, contig_name, start, end, **kwargs):
            called['contig_name'] = contig_name
            called['start'] = start
            called['end'] = end
            called['kwargs'] = kwargs
            return []

    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    fake_fasta = object()

    bamfile, pileup = methods.parse_bam(
        bamfile_name='sample.bam',
        contig_name='abc_1',
        pysam_fasta=fake_fasta,
        start=10,
        end=20
    )

    assert called['start'] == 10
    assert called['end'] == 20
    assert pileup == []


def test_characterise_read_in_fasta_mode_records_forward_ref(monkeypatch):
    class FakeAlignment:
        def __init__(self):
            self.qname = 'read1'
            self.query_sequence = 'ACGT'
            self.is_read1 = True
            self.is_read2 = False
            self.mate_is_unmapped = False
            self.is_paired = True
            self.mapping_quality = 30
            self.query_alignment_end = 4

    class FakePileupRead:
        def __init__(self):
            self.query_position = 0
            self.alignment = FakeAlignment()

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'contig1'
            self.pileups = [FakePileupRead()]

    rec = SimpleNamespace(letter_annotations={'phred_quality': [30, 20, 20, 20]})
    filtered, qualities, support = methods.characterise_read(
        column=FakeColumn(),
        reference_sequence='ACGT',
        fastq_records={'read1': rec},
        quality_cutoff=20,
        min_quality=15,
        fasta=True
    )

    assert filtered['forward_ref_reverse_UM_QF']['A'] == 1
    assert qualities == [30]
    assert 'A' in support


def test_read_contig_with_allele_records_and_fake_parse_bam(monkeypatch, tmp_path):
    allele_records = {'abc_1': SimpleNamespace(seq='ACGT')}

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'abc_1'
            self.pileups = []

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            return [b'abc_1']

    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), [FakeColumn()]))
    monkeypatch.setattr(methods, 'find_multibase_positions', lambda **kwargs: ({}, {}, 0, {}))
    monkeypatch.setattr(methods, 'benjamini_hochberg', lambda pvals: [0.05])
    monkeypatch.setattr(methods, '_position_entry_passes_probabilistic_gating', lambda position_stats: True)
    monkeypatch.setattr(methods, 'combine_pvalues_fisher', lambda pvals: 0.1)
    monkeypatch.setattr(methods, 'position_details', lambda **kwargs: 'abc_1\t1\tA\n')
    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta='dummy.fasta',
        allele_records=allele_records,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15
    )

    assert 'abc_1' in result
    assert '_gene_stats' in result['abc_1']
    assert report_text == ''


def test_read_contig_missing_reference_contig_returns_empty(monkeypatch, tmp_path):
    fasta = tmp_path / 'reference.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            return [b'wrong_contig']

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods.SeqIO, 'index', lambda path, fmt: {'wrong_contig': SimpleNamespace(seq='ACGT')})
    monkeypatch.setattr(methods.SeqIO, 'parse', lambda path, fmt: [SimpleNamespace(id='abc_1', seq='ACGT')])

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'abc_1'
            self.pileups = []

    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), [FakeColumn()]))
    monkeypatch.setattr(methods, 'find_multibase_positions', lambda **kwargs: ({}, {}, 0, {}))
    monkeypatch.setattr(methods, 'benjamini_hochberg', lambda pvals: [0.05])
    monkeypatch.setattr(methods, '_position_entry_passes_probabilistic_gating', lambda position_stats: True)
    monkeypatch.setattr(methods, 'combine_pvalues_fisher', lambda pvals: 0.1)
    monkeypatch.setattr(methods, 'position_details', lambda **kwargs: 'abc_1\t1\tA\n')

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records=None,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15
    )

    assert result == {}
    assert report_text == ''


def test_read_contig_falls_back_to_seqio_index_when_pysam_missing(monkeypatch, tmp_path):
    fasta = tmp_path / 'reference.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            raise AttributeError('pysam error')

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods.SeqIO, 'index', lambda path, fmt: {'abc_1': SimpleNamespace(seq='ACGT')})

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'abc_1'
            self.pileups = []

    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), [FakeColumn()]))
    monkeypatch.setattr(methods, 'find_multibase_positions', lambda **kwargs: ({}, {}, 0, {}))
    monkeypatch.setattr(methods, 'benjamini_hochberg', lambda pvals: [0.05])
    monkeypatch.setattr(methods, '_position_entry_passes_probabilistic_gating', lambda position_stats: True)
    monkeypatch.setattr(methods, 'combine_pvalues_fisher', lambda pvals: 0.1)
    monkeypatch.setattr(methods, 'position_details', lambda **kwargs: 'abc_1\t1\tA\n')

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records=None,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15
    )

    assert 'abc_1' in result
    assert report_text == ''


def test_read_contig_computes_dynamic_base_cutoff_when_zero(monkeypatch, tmp_path):
    fasta = tmp_path / 'reference.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            return [b'abc_1']

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    called = {'determine_cutoff': False}

    def fake_determine_cutoff(**kwargs):
        called['determine_cutoff'] = True
        return (2, 0.5, 0.1)

    monkeypatch.setattr(methods, 'determine_cutoff', fake_determine_cutoff)

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'abc_1'
            self.pileups = []

    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), [FakeColumn()]))
    monkeypatch.setattr(methods, 'find_multibase_positions', lambda **kwargs: ({}, {}, 0, {}))
    monkeypatch.setattr(methods, 'benjamini_hochberg', lambda pvals: [0.05])
    monkeypatch.setattr(methods, '_position_entry_passes_probabilistic_gating', lambda position_stats: True)
    monkeypatch.setattr(methods, 'combine_pvalues_fisher', lambda pvals: 0.1)
    monkeypatch.setattr(methods, 'position_details', lambda **kwargs: 'abc_1\t1\tA\n')

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records={'abc_1': SimpleNamespace(seq='ACGT')},
        fastq_records={},
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=0,
        base_fraction_cutoff=0.05
    )

    assert called['determine_cutoff'] is True
    assert 'abc_1' in result
    assert report_text == ''


def test_read_contig_fasta_index_creation_failure(monkeypatch, tmp_path):
    fasta = tmp_path / 'reference.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            return [b'wrong_contig']

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: (_ for _ in ()).throw(methods.SamtoolsError('index failed')))
    monkeypatch.setattr(methods.SeqIO, 'parse', lambda path, fmt: [SimpleNamespace(id='abc_1', seq='ACGT')])
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records=None,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15
    )

    assert result == {}
    assert report_text == ''


def test_read_contig_seqio_parse_failure(monkeypatch, tmp_path):
    fasta = tmp_path / 'reference.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            return [b'wrong_contig']

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods.SeqIO, 'parse', lambda path, fmt: (_ for _ in ()).throw(ValueError('parse failed')))
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records=None,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15
    )

    assert result == {}
    assert report_text == ''


def test_find_contamination_unpaired_illumina_bbmap_path(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_db = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_db.write_text('>abc_1\nACGT\n', encoding='utf-8')
    Path(str(sample_db) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    forward = tmp_path / 'reads.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            Path(str(tmp_path / 'sorted.bam')).write_bytes(b'')
        return 'out', 'err'

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: Path(str(tmp_path / 'sorted.bam')).write_bytes(b''))
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(path + '.bai').write_bytes(b''))
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', lambda *args, **kwargs: SimpleNamespace(references=[b'abc_1'], close=lambda: None))
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        threads=1,
        keep_files=False,
        data_type='Illumina',
        min_matching_hashes=40
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_find_contamination_direct_missing_allele_warning(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_db = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_db.write_text('>abc_1\nACGT\n', encoding='utf-8')
    Path(str(sample_db) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return ('out', 'err', 'cmd')

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return ('out2', 'err2', 'cmd2')

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1', 'missing_2'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', lambda cmd: ('out', 'err'))
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: Path(path + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8'))
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(path + '.bai').write_bytes(b''))
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', lambda *args, **kwargs: SimpleNamespace(references=[b'abc_1'], close=lambda: None))
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        threads=1,
        keep_files=False,
        min_matching_hashes=40
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_write_output_probabilistic_threshold():
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nTGCA\n+\nIIII\n', encoding='utf-8')

    def fake_index_db(path, fasta, fmt):
        return {'path': path}

    monkeypatch.setattr(methods.SeqIO, 'index_db', fake_index_db)
    methods._FASTQ_INDEX_STATE.clear()
    methods._init_fastq_index(
        fwd_path=str(forward),
        rev_path=str(reverse),
        paired=True,
        tmpdir=str(tmp_path)
    )

    assert methods._FASTQ_INDEX_STATE['paired'] is True
    assert methods._FASTQ_INDEX_STATE['fwd'] == {'path': str(tmp_path / 'reads_R1.fastq.fwd.idx')}
    assert methods._FASTQ_INDEX_STATE['rev'] == {'path': str(tmp_path / 'reads_R2.fastq.rev.idx')}


def test_get_fastq_record_suffix_fallback(monkeypatch):
    methods._FASTQ_INDEX_STATE.clear()
    methods._FASTQ_INDEX_STATE['paired'] = True
    methods._FASTQ_INDEX_STATE['fwd'] = {'read1': 'forward_record'}
    methods._FASTQ_INDEX_STATE['rev'] = {'read2/2': 'reverse_record'}

    assert methods._get_fastq_record('read1/1') == 'forward_record'
    assert methods._get_fastq_record('read2/2') == 'reverse_record'
    assert methods._get_fastq_record('read1') == 'forward_record'


def test_get_fastq_record_swallows_index_exceptions():
    class BadIndex:
        def __contains__(self, key):
            return True
        def __getitem__(self, key):
            raise KeyError('missing')

    original_state = methods._FASTQ_INDEX_STATE.copy()
    methods._FASTQ_INDEX_STATE.clear()
    methods._FASTQ_INDEX_STATE.update({'fwd': BadIndex(), 'rev': None, 'paired': False})
    try:
        assert methods._get_fastq_record('read1') is None
    finally:
        methods._FASTQ_INDEX_STATE.clear()
        methods._FASTQ_INDEX_STATE.update(original_state)


def test_init_fastq_index_handles_index_failure(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nTGCA\n+\nIIII\n', encoding='utf-8')

    def fake_index_db(path, fasta, fmt):
        raise ValueError('index failure')

    monkeypatch.setattr(methods.SeqIO, 'index_db', fake_index_db)
    methods._FASTQ_INDEX_STATE.clear()
    methods._init_fastq_index(
        fwd_path=str(forward),
        rev_path=str(reverse),
        paired=True,
        tmpdir=str(tmp_path)
    )

    assert methods._FASTQ_INDEX_STATE['fwd'] is None
    assert methods._FASTQ_INDEX_STATE['rev'] is None


def test_ensure_fastq_index_skips_existing_index(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nTGCA\n+\nIIII\n', encoding='utf-8')

    Path(str(forward) + '.fwd.idx').write_text('', encoding='utf-8')
    Path(str(reverse) + '.rev.idx').write_text('', encoding='utf-8')

    def fake_index_db(path, fasta, fmt):
        raise AssertionError('Should not build index when it already exists')

    monkeypatch.setattr(methods.SeqIO, 'index_db', fake_index_db)
    methods._ensure_fastq_index(
        fwd_path=str(forward),
        rev_path=str(reverse),
        paired=True,
        tmpdir=str(tmp_path)
    )


def test_index_databases_writes_log_on_kma_failure(monkeypatch, tmp_path):
    sample_database = tmp_path / 'Escherichia_db.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')

    def fake_faidx(path):
        Path(path + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    captured = {}

    def fake_write_to_logfile(logfile, out, err, cmd):
        captured['err'] = err
        captured['cmd'] = cmd

    def fake_run_cmd(cmd):
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(methods.pysam, 'faidx', fake_faidx)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods, 'write_to_logfile', fake_write_to_logfile)

    kma_database = methods.index_databases(sample_database=str(sample_database))

    assert kma_database.endswith('_kma')
    assert isinstance(captured['err'], subprocess.CalledProcessError)
    assert 'kma index' in captured['cmd']


def test_characterise_read_same_snv_and_one_sided_branches():
    class FakeAlignment:
        def __init__(self, qname, query_sequence, is_read1, is_read2):
            self.qname = qname
            self.query_sequence = query_sequence
            self.mapping_quality = 30
            self.is_read1 = is_read1
            self.is_read2 = is_read2
            self.mate_is_unmapped = False
            self.is_paired = True
            self.query_alignment_end = len(query_sequence)

    class FakePileupRead:
        def __init__(self, alignment, query_position=0):
            self.alignment = alignment
            self.query_position = query_position

    class FakeColumn:
        def __init__(self, pileups):
            self.pos = 0
            self.reference_name = 'contig1'
            self.pileups = pileups

    rec1 = SimpleNamespace(letter_annotations={'phred_quality': [30, 30, 30, 30, 30, 30]})
    rec2 = SimpleNamespace(letter_annotations={'phred_quality': [30, 30, 30, 30, 30, 30]})

    # Both reads support same SNV
    forward = FakeAlignment('read1', 'CAAAAA', True, False)
    reverse = FakeAlignment('read1', 'CAAAAA', False, True)
    column = FakeColumn([FakePileupRead(forward), FakePileupRead(reverse)])
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='AAAAAA',
        fastq_records={'read1/1': rec1, 'read1/2': rec2},
        quality_cutoff=20,
        min_quality=15,
        fasta=False
    )

    assert filtered['congruent_SNV']['C'] == 2
    assert qualities == [30, 30]
    assert 'C' in support

    # Forward SNV only, reverse matches reference
    forward = FakeAlignment('read1', 'CAAAAA', True, False)
    reverse = FakeAlignment('read1', 'AAAAAA', False, True)
    column = FakeColumn([FakePileupRead(forward), FakePileupRead(reverse)])
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='AAAAAA',
        fastq_records={'read1/1': rec1, 'read1/2': rec2},
        quality_cutoff=20,
        min_quality=15,
        fasta=False
    )

    assert filtered['forward_SNV_reverse_ref']['C'] == 1
    assert filtered['forward_SNV_reverse_ref']['A'] == 1

    # Reverse SNV only, forward matches reference
    forward = FakeAlignment('read1', 'AAAAAA', True, False)
    reverse = FakeAlignment('read1', 'CAAAAA', False, True)
    column = FakeColumn([FakePileupRead(forward), FakePileupRead(reverse)])
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='AAAAAA',
        fastq_records={'read1/1': rec1, 'read1/2': rec2},
        quality_cutoff=20,
        min_quality=15,
        fasta=False
    )

    assert filtered['reverse_SNV_forward_ref']['C'] == 1
    assert filtered['reverse_SNV_forward_ref']['A'] == 1


def test_characterise_read_quality_filtered_both_reads(monkeypatch):
    class FakeAlignment:
        def __init__(self, qname, seq, is_read1, is_read2):
            self.qname = qname
            self.query_sequence = seq
            self.mapping_quality = 10
            self.is_read1 = is_read1
            self.is_read2 = is_read2
            self.mate_is_unmapped = False
            self.is_paired = True
            self.query_alignment_end = len(seq)

    class FakePileupRead:
        def __init__(self, alignment, query_position=0):
            self.query_position = query_position
            self.alignment = alignment

    class FakeColumn:
        def __init__(self, pileups):
            self.pos = 0
            self.reference_name = 'gene1'
            self.pileups = pileups

    forward = FakeAlignment('read1', 'C', True, False)
    reverse = FakeAlignment('read1', 'G', False, True)
    column = FakeColumn([FakePileupRead(forward), FakePileupRead(reverse)])
    rec1 = SimpleNamespace(letter_annotations={'phred_quality': [10]})
    rec2 = SimpleNamespace(letter_annotations={'phred_quality': [10]})

    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='A',
        fastq_records={'read1/1': rec1, 'read1/2': rec2},
        quality_cutoff=20,
        min_quality=15,
        fasta=False
    )

    assert filtered['forward_quality_filtered']['C'] == 1
    assert filtered['reverse_quality_filtered']['G'] == 1
    assert qualities == []
    assert support == {}


def test_read_contig_fasta_index_creation_failure(monkeypatch, tmp_path):
    fasta = tmp_path / 'reference.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            return [b'wrong_contig']

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: (_ for _ in ()).throw(methods.SamtoolsError('index failed')))
    monkeypatch.setattr(methods, 'SeqIO', methods.SeqIO)
    monkeypatch.setattr(methods.SeqIO, 'parse', lambda path, fmt: [SimpleNamespace(id='abc_1', seq='ACGT')])
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records=None,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15
    )

    assert result == {}
    assert report_text == ''


def test_read_contig_seqio_parse_failure(monkeypatch, tmp_path):
    fasta = tmp_path / 'reference.fasta'
    fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeFastaFile:
        def __init__(self, path):
            pass
        def close(self):
            pass
        @property
        def references(self):
            return [b'wrong_contig']

    monkeypatch.setattr(methods.pysam, 'FastaFile', FakeFastaFile)
    monkeypatch.setattr(methods.SeqIO, 'parse', lambda path, fmt: (_ for _ in ()).throw(ValueError('parse failed')))
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(fasta),
        allele_records=None,
        fastq_records={},
        quality_cutoff=20,
        min_quality=15
    )

    assert result == {}
    assert report_text == ''


def test_find_contamination_parallel_pool_branch_uses_ensure_fastq_index(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_db = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_db.write_text('>abc_1\nACGT\n', encoding='utf-8')
    Path(str(sample_db) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    def fake_run_cmd(cmd):
        outbam = tmp_path / 'outbam.bam'
        outbam.write_bytes(b'')
        return 'out', 'err'

    class FakePool:
        def __init__(self, processes, initializer=None, initargs=None):
            self.processes = processes
            self.initializer = initializer
            self.initargs = initargs
            self.closed = False
        def map(self, func, iterable, chunksize=1):
            return [({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n')]
        def close(self):
            self.closed = True
        def join(self):
            self.closed = True

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    def fake_load_fastq_records(*args, **kwargs):
        if kwargs.get('forward'):
            return {'read1/1': SimpleNamespace()}
        return {'read1/2': SimpleNamespace()}
    monkeypatch.setattr(methods, 'load_fastq_records', fake_load_fastq_records)
    monkeypatch.setattr(methods, '_ensure_fastq_index', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: Path(args[1]).write_bytes(b''))
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(path + '.bai').write_bytes(b''))
    class FakeAlignmentFile:
        def __init__(self, *args, **kwargs):
            self.references = [b'abc_1']
        def close(self):
            pass
        def __enter__(self):
            return self
        def __exit__(self, exc_type, exc, tb):
            return False
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    monkeypatch.setattr(methods.multiprocessing, 'Pool', FakePool)
    monkeypatch.setattr(methods, '_read_contig_chunk_dispatch', lambda kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        threads=2,
        keep_files=False,
        min_matching_hashes=40
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_write_output_probabilistic_threshold():
    report = Path('prob_report.tsv')
    try:
        methods.write_output(
            output_report=str(report),
            sample_name='sample1',
            multi_positions=0,
            genus='Escherichia',
            total_gene_length=1000,
            database_download_date='2026-07-08',
            pysam_pass=True,
            sample_score=5.0,
            use_probabilistic=True,
            score_threshold=2.0
        )

        content = report.read_text(encoding='utf-8').splitlines()
        assert content[1].split('\t')[3] == 'True'
    finally:
        if report.exists():
            report.unlink()


def test_downsample_reads_includes_seed(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nTGCA\n+\nIIII\n', encoding='utf-8')

    commands = []

    def fake_run_cmd(cmd):
        commands.append(cmd)
        return 'out', 'err'

    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)

    methods.downsample_reads(
        pair=[str(forward), str(reverse)],
        sample_tmp_dir=str(tmp_path),
        sample_name='sample',
        target_reads=10,
        log=str(tmp_path / 'downsample.log'),
        xmx='4g',
        seed=42
    )

    assert 'sampleseed=42' in commands[0]
    assert '-Xmx=4g' in commands[0] or '-Xmx4g' in commands[0]


def test_find_unpaired_reads_fasta(tmp_path):
    fasta = tmp_path / 'sample_R1.fasta'
    fasta.write_text('>seq1\nACGT\n', encoding='utf-8')
    results = find_unpaired_reads(
        fastq_directory=str(tmp_path),
        find_fasta=True
    )
    assert results == [[str(fasta)]]


def test_download_mash_sketch(monkeypatch, tmp_path):
    def fake_urlretrieve(url, dest):
        Path(dest).write_text('sketch', encoding='utf-8')
        return url, dest

    monkeypatch.setattr(methods.urllib.request, 'urlretrieve', fake_urlretrieve)
    methods.download_mash_sketch(output_folder=str(tmp_path))
    assert (tmp_path / 'refseq.msh').read_text(encoding='utf-8') == 'sketch'


def test_download_cgmlst_derived_data_extracts_and_indexes(monkeypatch, tmp_path):
    tar_path = tmp_path / 'confindr_db.tar.gz'

    def fake_urlretrieve(url, dest):
        with tarfile.open(dest, 'w:gz') as tar:
            gene_allele = tmp_path / 'gene_allele.txt'
            rmlst = tmp_path / 'rMLST_combined.fasta'
            gene_allele.write_text('Escherichia:abc_1,\n', encoding='utf-8')
            rmlst.write_text('>abc_1\nACGT\n', encoding='utf-8')
            tar.add(str(gene_allele), arcname='gene_allele.txt')
            tar.add(str(rmlst), arcname='rMLST_combined.fasta')
        return url, dest

    monkeypatch.setattr(methods.urllib.request, 'urlretrieve', fake_urlretrieve)
    called = {}

    def fake_index(output_folder, genera, cgderived):
        called['args'] = (output_folder, tuple(genera), cgderived)

    monkeypatch.setattr(methods, 'index', fake_index)
    methods.download_cgmlst_derived_data(output_folder=str(tmp_path))
    assert called['args'] == (str(tmp_path), ('Escherichia', 'Listeria', 'Salmonella'), True)
    assert not tar_path.exists()
    assert (tmp_path / 'gene_allele.txt').exists()
    assert (tmp_path / 'rMLST_combined.fasta').exists()


def test_download_cgmlst_derived_data_ignores_remove_failure(monkeypatch, tmp_path):
    tarball = tmp_path / 'confindr_db.tar.gz'

    def fake_urlretrieve(url, dest):
        with tarfile.open(dest, 'w:gz') as tar:
            tmp_file = tmp_path / 'file.txt'
            tmp_file.write_text('content', encoding='utf-8')
            tar.add(str(tmp_file), arcname='file.txt')
        return url, dest

    def fake_index(output_folder, genera, cgderived):
        pass

    monkeypatch.setattr(methods.urllib.request, 'urlretrieve', fake_urlretrieve)
    monkeypatch.setattr(methods, 'index', fake_index)
    monkeypatch.setattr(methods.os, 'remove', lambda dest: (_ for _ in ()).throw(OSError('remove failed')))

    methods.download_cgmlst_derived_data(output_folder=str(tmp_path))
    assert (tmp_path / 'file.txt').exists()


def test_index_creates_and_indexes_database(monkeypatch, tmp_path):
    profiles = tmp_path / 'gene_allele.txt'
    profiles.write_text('Escherichia:abc_1,\n', encoding='utf-8')
    rmlst = tmp_path / 'rMLST_combined.fasta'
    rmlst.write_text('>abc_1\nACGT\n', encoding='utf-8')

    called = {'find': False, 'setup': False, 'index': False}

    def fake_find_genus_specific_allele_list(profiles_file, target_genus):
        called['find'] = True
        assert target_genus == 'Escherichia'
        return ['abc_1']

    def fake_setup_allelespecific_database(fasta_file, database_folder, allele_list):
        called['setup'] = True
        assert allele_list == ['abc_1']
        assert database_folder == str(tmp_path)

    def fake_index_databases(sample_database):
        called['index'] = True
        assert sample_database == str(tmp_path / 'Escherichia_db_cgderived.fasta')

    monkeypatch.setattr(methods, 'find_genus_specific_allele_list', fake_find_genus_specific_allele_list)
    monkeypatch.setattr(methods, 'setup_allelespecific_database', fake_setup_allelespecific_database)
    monkeypatch.setattr(methods, 'index_databases', fake_index_databases)

    methods.index(output_folder=str(tmp_path), genera=['Escherichia'], cgderived=True)

    assert called['find']
    assert called['setup']
    assert called['index']


def test_index_skips_database_creation_if_database_exists(monkeypatch, tmp_path):
    db_file = tmp_path / 'Escherichia_db_cgderived.fasta'
    db_file.write_text('>abc_1\nACGT\n', encoding='utf-8')

    def fake_index_databases(sample_database):
        assert sample_database == str(db_file)

    monkeypatch.setattr(methods, 'find_genus_specific_allele_list', lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError('should not be called')))
    monkeypatch.setattr(methods, 'setup_allelespecific_database', lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError('should not be called')))
    monkeypatch.setattr(methods, 'index_databases', fake_index_databases)

    methods.index(output_folder=str(tmp_path), genera=['Escherichia'], cgderived=True)


def test_get_fastq_record_with_paired_indices():
    original_state = methods._FASTQ_INDEX_STATE.copy()
    methods._FASTQ_INDEX_STATE.update({
        'fwd': {'read1': 'FWD', 'read1/1': 'FWD1'},
        'rev': {'read1/2': 'REV2'},
        'paired': True
    })
    try:
        assert methods._get_fastq_record('read1/1') == 'FWD1'
        assert methods._get_fastq_record('read1/2') == 'REV2'
        assert methods._get_fastq_record('read1') == 'FWD'
        assert methods._get_fastq_record('missing/1') is None
    finally:
        methods._FASTQ_INDEX_STATE.clear()
        methods._FASTQ_INDEX_STATE.update(original_state)


def test_get_fastq_record_returns_none_without_index():
    original_state = methods._FASTQ_INDEX_STATE.copy()
    methods._FASTQ_INDEX_STATE.update({'fwd': None, 'rev': None, 'paired': False})
    try:
        assert methods._get_fastq_record('read1/1') is None
    finally:
        methods._FASTQ_INDEX_STATE.clear()
        methods._FASTQ_INDEX_STATE.update(original_state)


def test_index_databases_creates_kma_and_logs(monkeypatch, tmp_path):
    sample_database = tmp_path / 'sample.fasta'
    sample_database.write_text('>seq1\nACGT\n', encoding='utf-8')
    fai_path = tmp_path / 'sample.fasta.fai'
    name_path = tmp_path / 'sample_kma.name'

    def fake_faidx(path):
        assert path == str(sample_database)
        fai_path.write_text('index', encoding='utf-8')

    called = {'run': False, 'log': False}

    def fake_run_cmd(cmd):
        called['run'] = True
        assert 'kma index' in cmd
        name_path.write_text('name', encoding='utf-8')
        return 'out', 'err'

    def fake_write_to_logfile(logfile, out, err, cmd):
        called['log'] = True
        assert logfile == str(sample_database) + '_log.txt'
        assert out == 'out'
        assert err == 'err'
        assert 'kma index' in cmd

    monkeypatch.setattr(methods.pysam, 'faidx', fake_faidx)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods, 'write_to_logfile', fake_write_to_logfile)

    result = methods.index_databases(sample_database=str(sample_database))
    assert result.endswith('_kma')
    assert called['run'] is True
    assert called['log'] is True
    assert name_path.exists()


def test_index_databases_ignores_faidx_failure(monkeypatch, tmp_path):
    sample_database = tmp_path / 'sample.fasta'
    sample_database.write_text('>seq1\nACGT\n', encoding='utf-8')

    def fake_faidx(path):
        raise methods.SamtoolsError('faidx failure')

    called = {'run': False}

    def fake_run_cmd(cmd):
        called['run'] = True
        return 'out', 'err'

    monkeypatch.setattr(methods.pysam, 'faidx', fake_faidx)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)

    result = methods.index_databases(sample_database=str(sample_database))
    assert result.endswith('_kma')
    assert called['run'] is True

    assert (str(sample_database).replace('.fasta', '') + '_kma') == result


def test_index_databases_skips_kma_when_already_indexed(monkeypatch, tmp_path):
    sample_database = tmp_path / 'sample.fasta'
    sample_database.write_text('>seq1\nACGT\n', encoding='utf-8')
    fai_path = tmp_path / 'sample.fasta.fai'
    name_path = tmp_path / 'sample_kma.name'
    fai_path.write_text('index', encoding='utf-8')
    name_path.write_text('name', encoding='utf-8')

    monkeypatch.setattr(methods, 'run_cmd', lambda cmd: (_ for _ in ()).throw(AssertionError('run_cmd should not be called')))
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError('write_to_logfile should not be called')))

    result = methods.index_databases(sample_database=str(sample_database))
    assert result == str(tmp_path / 'sample_kma')


def test_ensure_fastq_index_builds_indices(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nTGCA\n+\nIIII\n', encoding='utf-8')

    def fake_index_db(path, fasta, fmt):
        Path(path).write_text('', encoding='utf-8')

    monkeypatch.setattr(methods.SeqIO, 'index_db', fake_index_db)

    methods._ensure_fastq_index(
        fwd_path=str(forward),
        rev_path=str(reverse),
        paired=True,
        tmpdir=str(tmp_path)
    )

    assert (tmp_path / 'reads_R1.fastq.fwd.idx').exists()
    assert (tmp_path / 'reads_R2.fastq.rev.idx').exists()


def test_characterise_read_for_single_forward_read():
    class FakeAlignment:
        def __init__(self):
            self.qname = 'read1'
            self.query_sequence = 'ACGT'
            self.is_read1 = True
            self.is_read2 = False
            self.mate_is_unmapped = False
            self.is_paired = True
            self.mapping_quality = 30
            self.query_alignment_end = 4

    class FakePileupRead:
        def __init__(self):
            self.query_position = 0
            self.alignment = FakeAlignment()

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'contig1'
            self.pileups = [FakePileupRead()]

    rec = SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]})
    filtered, qualities, support = methods.characterise_read(
        column=FakeColumn(),
        reference_sequence='ACGT',
        fastq_records={'read1/1': rec},
        quality_cutoff=20,
        min_quality=15
    )

    assert filtered['forward_ref_reverse_UM_QF'] == {'A': 1}
    assert qualities == [20]
    assert 'A' in support


def make_pileup_column(reads):
    class FakeAlignment:
        def __init__(self, qname, sequence, is_read1, mapping_quality=30):
            self.qname = qname
            self.query_sequence = sequence
            self.is_read1 = is_read1
            self.is_read2 = not is_read1
            self.mate_is_unmapped = False
            self.is_paired = True
            self.mapping_quality = mapping_quality
            self.query_alignment_end = len(sequence)

    class FakePileupRead:
        def __init__(self, alignment):
            self.query_position = 0
            self.alignment = alignment

    class FakeColumn:
        def __init__(self, reads):
            self.pos = 0
            self.reference_name = 'contig1'
            self.pileups = [FakePileupRead(FakeAlignment(*r)) for r in reads]

    return FakeColumn(reads)


def test_characterise_read_pair_same_snv_counts_as_congruent_snv():
    reads = [
        ('read1/1', 'CCGT', True),
        ('read1/2', 'CCGT', False)
    ]
    column = make_pileup_column(reads)
    recs = {
        'read1/1': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]}),
        'read1/2': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]})
    }
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='ACGT',
        fastq_records=recs,
        quality_cutoff=20,
        min_quality=15
    )

    assert filtered['forward_SNV_reverse_UM_QF'] == {'C': 1}
    assert filtered['reverse_SNV_forward_UM_QF'] == {'C': 1}
    assert qualities == [20, 20]
    assert 'C' in support


def test_characterise_read_pair_different_snvs_creates_snv1_groups():
    reads = [
        ('read1/1', 'CCGT', True),
        ('read1/2', 'GCGT', False)
    ]
    column = make_pileup_column(reads)
    recs = {
        'read1/1': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]}),
        'read1/2': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]})
    }
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='ACGT',
        fastq_records=recs,
        quality_cutoff=20,
        min_quality=15
    )

    assert filtered['forward_SNV_reverse_UM_QF'] == {'C': 1}
    assert filtered['reverse_SNV_forward_UM_QF'] == {'G': 1}
    assert qualities == [20, 20]


def test_characterise_read_pair_one_side_snv_and_one_side_ref():
    reads = [
        ('read1/1', 'CCGT', True),
        ('read1/2', 'ACGT', False)
    ]
    column = make_pileup_column(reads)
    recs = {
        'read1/1': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]}),
        'read1/2': SimpleNamespace(letter_annotations={'phred_quality': [10, 20, 20, 20]})
    }
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='ACGT',
        fastq_records=recs,
        quality_cutoff=20,
        min_quality=15
    )

    assert filtered['forward_SNV_reverse_UM_QF'] == {'C': 1}
    assert filtered['reverse_quality_filtered'] == {'A': 1}


def test_characterise_read_pair_both_match_reference_counts_congruent_ref():
    reads = [
        ('read1/1', 'ACGT', True),
        ('read1/2', 'ACGT', False)
    ]
    column = make_pileup_column(reads)
    recs = {
        'read1/1': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]}),
        'read1/2': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]})
    }
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='ACGT',
        fastq_records=recs,
        quality_cutoff=20,
        min_quality=15
    )

    assert filtered['forward_ref_reverse_UM_QF'] == {'A': 1}
    assert filtered['reverse_ref_forward_UM_QF'] == {'A': 1}
    assert qualities == [20, 20]


def test_characterise_read_pair_same_snv_counts_as_congruent_snv():
    reads = [
        ('read1', 'CCGT', True),
        ('read1', 'CCGT', False)
    ]
    column = make_pileup_column(reads)
    recs = {
        'read1/1': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]}),
        'read1/2': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]})
    }
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='ACGT',
        fastq_records=recs,
        quality_cutoff=20,
        min_quality=15
    )

    assert filtered['congruent_SNV'] == {'C': 2}
    assert qualities == [20, 20]
    assert support['C']['forward'] == 1
    assert support['C']['reverse'] == 1


def test_characterise_read_pair_different_snvs_both_high_quality():
    reads = [
        ('read1', 'CCGT', True),
        ('read1', 'GCGT', False)
    ]
    column = make_pileup_column(reads)
    recs = {
        'read1/1': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]}),
        'read1/2': SimpleNamespace(letter_annotations={'phred_quality': [20, 20, 20, 20]})
    }
    filtered, qualities, support = methods.characterise_read(
        column=column,
        reference_sequence='ACGT',
        fastq_records=recs,
        quality_cutoff=20,
        min_quality=15
    )

    assert filtered['forward_SNV_reverse_SNV1'] == {'C': 1}
    assert filtered['reverse_SNV_forward_SNV1'] == {'G': 1}
    assert qualities == [20, 20]


def test_read_contig_returns_empty_for_missing_contig(tmp_path):
    reference_fasta = tmp_path / 'ref.fasta'
    reference_fasta.write_text('>chr1\nACGT\n', encoding='utf-8')

    result, text = methods.read_contig(
        contig_name='chr2',
        bamfile_name=str(tmp_path / 'dummy.bam'),
        reference_fasta=str(reference_fasta)
    )
    assert result == {}
    assert text == ''


def test_bbtools_kwargs_to_string():
    assert bbtools.kwargs_to_string({'Xmx': '4g', 'threads': 8}) == ' Xmx=4g threads=8'


def test_mash_kwargs_to_string():
    assert mash.kwargs_to_string({'k': 21, 's': ''}) == ' -k 21 -s '


def test_bbtools_run_subprocess_success_and_failure(monkeypatch):
    class FakeProcess:
        def __init__(self, command, shell, stdout, stderr):
            self.returncode = 0
            self.command = command
        def communicate(self):
            return b'out', b'err'
    monkeypatch.setattr(bbtools, 'Popen', FakeProcess)
    out, err = bbtools.run_subprocess('echo hi')
    assert out == 'out'
    assert err == 'err'

    class FakeProcessFail(FakeProcess):
        def __init__(self, command, shell, stdout, stderr):
            super().__init__(command, shell, stdout, stderr)
            self.returncode = 1
    monkeypatch.setattr(bbtools, 'Popen', FakeProcessFail)
    with pytest.raises(Exception):
        bbtools.run_subprocess('false')


def test_bbtools_bbmap_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_run_subprocess(command):
        return 'stdout', 'stderr'
    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)
    out, err, cmd = bbtools.bbmap(
        reference='ref.fasta',
        forward_in=str(forward),
        out_bam=str(tmp_path / 'out.bam'),
        returncmd=True
    )
    assert 'in2=' in cmd
    assert 'out=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_bbduck_trim_autodetects_pair(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_check_output(cmd):
        return b'/usr/bin/bbduk.sh\n'
    def fake_run_subprocess(command):
        return 'stdout', 'stderr'

    monkeypatch.setattr(bbtools.subprocess, 'check_output', fake_check_output)
    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)

    forward_out = tmp_path / 'reads_R1_trimmed.fastq'
    reverse_out = tmp_path / 'reads_R2_trimmed.fastq'
    out, err = bbtools.bbduk_trim(
        forward_in=str(forward),
        forward_out=str(forward_out),
        reverse_in='NA',
        reverse_out=str(reverse_out),
        returncmd=False
    )
    assert out == 'stdout'
    assert err == 'stderr'


@pytest.mark.parametrize('func,args,expected', [
    (mash.sketch, ('file.fasta',), 'mash sketch file.fasta -o sketch.msh -p 1 '),
    (mash.dist, ('a.msh', 'b.msh'), 'mash dist a.msh b.msh  -p 1  > distances.tab'),
    (mash.screen, ('a.msh', 'b.msh'), 'mash screen a.msh b.msh  -p 1  | sort -gr > screen.tab'),
])
def test_mash_command_strings(monkeypatch, func, args, expected):
    def fake_run_subprocess(command):
        return '', ''
    monkeypatch.setattr(mash, 'run_subprocess', fake_run_subprocess)
    out, err, cmd = func(*args, returncmd=True)
    assert expected in cmd


def test_read_mash_output_and_screen(tmp_path):
    mash_file = tmp_path / 'dist.tab'
    mash_file.write_text('ref query 0.01 1e-5 50\n', encoding='utf-8')
    results = mash.read_mash_output(str(mash_file))
    assert len(results) == 1
    assert results[0].reference == 'ref'

    screen_file = tmp_path / 'screen.tab'
    screen_file.write_text('99.9 100 1 1e-10 query1\n', encoding='utf-8')
    screen_results = mash.read_mash_screen(str(screen_file))
    assert len(screen_results) == 1
    assert screen_results[0].query_id == 'query1'


def test_mash_run_subprocess(monkeypatch):
    class FakeProcess:
        def __init__(self, command, shell, stdout, stderr):
            self.returncode = 0
        def communicate(self):
            return b'out', b'err'
    monkeypatch.setattr(mash, 'Popen', FakeProcess)
    out, err = mash.run_subprocess('echo hi')
    assert out == 'out'
    assert err == 'err'


def test_mash_sketch_and_dist_require_args():
    with pytest.raises(ValueError):
        mash.sketch(returncmd=True)
    with pytest.raises(ValueError):
        mash.dist(returncmd=True)


def test_bbtools_multiple_wrappers_and_auto_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_run_subprocess(command):
        return 'stdout', 'stderr'
    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)
    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda *args, **kwargs: b'/usr/bin/bbduk.sh')

    out, err, cmd = bbtools.tadpole(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'tadpole_R1.fastq'),
        returncmd=True
    )
    assert 'tadpole.sh' in cmd
    assert out == 'stdout'

    out, err, cmd = bbtools.bbnorm(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'bbnorm_R1.fastq'),
        returncmd=True
    )
    assert 'bbnorm.sh' in cmd

    out, err, cmd = bbtools.bbmerge(
        forward_in=str(forward),
        merged_reads=str(tmp_path / 'merged.fastq'),
        returncmd=True
    )
    assert 'bbmerge.sh' in cmd

    out, err, cmd = bbtools.bbduk_bait(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'baited_R1.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh' in cmd

    out, err, cmd = bbtools.bbduk_trim(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'trimmed_R1.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh' in cmd

    out, err, cmd = bbtools.dedupe(
        input_file=str(forward),
        output_file=str(tmp_path / 'deduped.fastq'),
        returncmd=True
    )
    assert 'dedupe.sh' in cmd

    out, err, cmd = bbtools.seal(
        reference='ref.fasta',
        forward_in=str(forward),
        output_file=str(tmp_path / 'seal.rpkm'),
        returncmd=True
    )
    assert 'seal.sh' in cmd

    out, err, cmd = bbtools.kmercountexact(
        forward_in=str(forward),
        returncmd=True
    )
    assert 'kmercountexact.sh' in cmd

    peaks_file = tmp_path / 'peaks.txt'
    peaks_file.write_text('#haploid_genome_size 12345\n', encoding='utf-8')
    assert bbtools.genome_size(str(peaks_file), haploid=True) == 12345

    out, err, cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'subsampled_R1.fastq'),
        num_bases=100,
        returncmd=True
    )
    assert 'reformat.sh' in cmd

    out, err, cmd = bbtools.validate_reads(
        forward_in=str(forward),
        returncmd=True
    )
    assert 'reformat.sh' in cmd

    out, err, cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'reformatted_R1.fastq'),
        returncmd=True
    )
    assert 'reformat.sh' in cmd

    out, err, cmd = bbtools.repair_reads(
        forward_in=str(forward),
        reverse_in=str(reverse),
        returncmd=True,
        forward_out=str(tmp_path / 'repair_R1.fastq'),
        reverse_out=str(tmp_path / 'repair_R2.fastq')
    )
    assert 'repair.sh' in cmd


def test_bbtools_bbmap_no_returncmd(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err = bbtools.bbmap(
        reference='ref.fasta',
        forward_in=str(forward),
        out_bam=str(tmp_path / 'out.bam')
    )
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_bbduk_trim_single_end(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda *args, **kwargs: b'/usr/bin/bbduk.sh\n')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))

    out, err, cmd = bbtools.bbduk_trim(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'trimmed.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh in=' in cmd
    assert 'out=' in cmd
    assert out == 'out'
    assert err == 'err'


def test_bbtools_bbduk_trim_requires_reverse_out(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda *args, **kwargs: b'/usr/bin/bbduk.sh\n')
    with pytest.raises(ValueError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'trimmed_R1.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_tadpole_output_exists_returns_empty(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output_file = tmp_path / 'tadpole.fastq'
    output_file.write_text('', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.tadpole(
        forward_in=str(forward),
        forward_out=str(output_file),
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'tadpole.sh' in cmd


def test_bbtools_bbduk_bait_single_end(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.bbduk_bait(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'baited.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh qin=33 in=' in cmd
    assert 'outm=' in cmd


def test_bbtools_bbduk_filter_single_end(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.bbduk_filter(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'filtered.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh qin=33 in=' in cmd
    assert 'out=' in cmd


def test_bbtools_bbnorm_output_exists_returns_empty(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output_file = tmp_path / 'bbnorm.fastq'
    output_file.write_text('', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.bbnorm(
        forward_in=str(forward),
        forward_out=str(output_file),
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'bbnorm.sh' in cmd


def test_bbtools_repair_reads_requires_reverse_out(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.repair_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'repair_R1.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_reformat_reads_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'reformatted_R1.fastq'),
        reverse_in=str(reverse),
        reverse_out=str(tmp_path / 'reformatted_R2.fastq'),
        returncmd=True
    )
    assert 'reformat.sh in1=' in cmd
    assert 'out2=' in cmd


def test_bbtools_bbduk_trim_raises_file_not_found_when_binary_missing(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    forward_out = tmp_path / 'trimmed_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_check_output(cmd):
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(bbtools.subprocess, 'check_output', fake_check_output)
    with pytest.raises(FileNotFoundError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(forward_out),
            reverse_in=str(reverse),
            reverse_out=str(tmp_path / 'trimmed_R2.fastq')
        )


def test_bbtools_bbmap_explicit_reverse_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbmap(
        reference='ref.fasta',
        forward_in=str(forward),
        reverse_in=str(reverse),
        out_bam=str(tmp_path / 'out.bam'),
        returncmd=True
    )
    assert 'bbmap.sh ref=' in cmd
    assert 'in2=' in cmd


def test_bbtools_bbnorm_raises_reverse_out_missing(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.bbnorm(
            forward_in=str(forward),
            reverse_in=str(reverse),
            forward_out=str(tmp_path / 'bbnorm_R1.fastq'),
            reverse_out='NA'
        )


def test_bbtools_tadpole_explicit_reverse_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.tadpole(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'tadpole_R1.fastq'),
        reverse_in=str(reverse),
        reverse_out=str(tmp_path / 'tadpole_R2.fastq'),
        returncmd=True
    )
    assert 'tadpole.sh in1=' in cmd
    assert 'out2=' in cmd


def test_bbtools_bbduk_bait_explicit_reverse_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.bbduk_bait(
        reference='ref.fasta',
        forward_in=str(forward),
        reverse_in=str(reverse),
        forward_out=str(tmp_path / 'baited_R1.fastq'),
        reverse_out=str(tmp_path / 'baited_R2.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh qin=33 in=' in cmd
    assert 'outm2=' in cmd


def test_bbtools_bbduk_filter_explicit_reverse_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.bbduk_filter(
        reference='ref.fasta',
        forward_in=str(forward),
        reverse_in=str(reverse),
        forward_out=str(tmp_path / 'filter_R1.fastq'),
        reverse_out=str(tmp_path / 'filter_R2.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh qin=33 in=' in cmd
    assert 'in2=' in cmd
    assert 'out2=' in cmd


def test_read_contig_chunk_dispatch_merges_results(monkeypatch):
    def fake_read_contig(**kwargs):
        contig = kwargs['contig_name']
        report = f'{contig}\t1\tA\n'
        return ({contig: {1: {}}}, report)

    monkeypatch.setattr(methods, 'read_contig', fake_read_contig)

    combined, text = methods._read_contig_chunk_dispatch({
        'contig_chunk': [
            {'contig': 'a', 'start': 0, 'end': 10},
            {'contig': 'b', 'start': 10, 'end': 20}
        ],
        'bamfile_name': 'dummy.bam',
        'reference_fasta': 'dummy.fasta',
        'fastq_records': {},
        'quality_cutoff': 20,
        'base_cutoff': 1,
        'base_fraction_cutoff': 0.1,
        'fasta': False,
        'error_cutoff': 1.0,
        'nanopore': False,
        'max_expected_positions': 0.001
    })

    assert 'a' in combined and 'b' in combined
    assert 'a\t1\tA' in text
    assert 'b\t1\tA' in text


def test_bbtools_reformat_reads_requires_reverse_out_if_reverse_in_provided(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.reformat_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'reformatted.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_bbduk_trim_autodetect_reverse_requires_forward_out_r1(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda *args, **kwargs: b'/usr/bin/bbduk.sh\n')
    with pytest.raises(ValueError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'trimmed.fastq'),
            returncmd=True
        )


def test_bbtools_tadpole_requires_reverse_out_when_forward_out_not_r1(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.tadpole(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'tadpole.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_bbnorm_requires_reverse_out_when_forward_out_not_r1(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.bbnorm(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'bbnorm.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_bbmerge_output_exists_returns_empty(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    merged = tmp_path / 'merged.fastq'
    merged.write_text('', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.bbmerge(
        forward_in=str(forward),
        merged_reads=str(merged),
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'bbmerge.sh' in cmd


def test_bbtools_subsample_reads_requires_reverse_out_if_reverse_in_provided(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    with pytest.raises(ValueError):
        bbtools.subsample_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'subsampled.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA',
            num_bases=100
        )


def test_bbtools_reformat_reads_requires_reverse_out_if_reverse_in_provided(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.repair_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'repair_R1.fastq'),
        reverse_in='NA',
        reverse_out='NA',
        returncmd=True
    )
    assert 'repair.sh' in cmd
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.validate_reads(
        forward_in=str(forward),
        reverse_in=str(reverse),
        returncmd=True
    )
    assert 'reformat.sh in1=' in cmd
    assert 'in2=' in cmd


def test_bbtools_kmercountexact_explicit_reverse_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.kmercountexact(
        forward_in=str(forward),
        reverse_in=str(reverse),
        returncmd=True
    )
    assert 'kmercountexact.sh in=' in cmd
    assert 'in2=' in cmd


def test_bbtools_subsample_reads_single_end(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'subsampled.fastq'),
        num_bases=10,
        returncmd=True
    )
    assert 'reformat.sh in=' in cmd


def test_bbtools_reformat_reads_output_exists_empty(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output_file = tmp_path / 'reformatted.fastq'
    output_file.write_text('', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(output_file),
        returncmd=True
    )
    assert out == ''
    assert err == ''


def test_bbtools_repair_reads_output_exists_empty(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output_file = tmp_path / 'repair_R1.fastq'
    output_file.write_text('', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.repair_reads(
        forward_in=str(forward),
        forward_out=str(output_file),
        reverse_in=str(reverse),
        reverse_out=str(tmp_path / 'repair_R2.fastq'),
        returncmd=True
    )
    assert out == ''
    assert err == ''


def test_bbtools_validate_reads_auto_detects_reverse_use_cmd(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.validate_reads(
        forward_in=str(forward),
        returncmd=True
    )
    assert 'reformat.sh in1=' in cmd
    assert 'in2=' in cmd


def test_bbtools_reformat_reads_auto_detects_reverse_use_cmd(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'reformatted_R1.fastq'),
        returncmd=True
    )
    assert 'reformat.sh in1=' in cmd
    assert 'in2=' in cmd


def test_bbtools_subsample_reads_output_exists_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output_file = tmp_path / 'subsampled_R1.fastq'
    output_file.write_text('', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('out', 'err'))
    out, err, cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(output_file),
        reverse_out=str(tmp_path / 'subsampled_R2.fastq'),
        num_bases=100,
        returncmd=True
    )
    assert out == ''
    assert err == ''
    assert 'reformat.sh' in cmd


def test_bbtools_run_subprocess_raises_calledprocesserror_on_failure(monkeypatch):
    class FakeProcess:
        def __init__(self, command, shell, stdout, stderr):
            self.returncode = 1
        def communicate(self):
            return b'', b'err'

    monkeypatch.setattr(bbtools, 'Popen', FakeProcess)
    with pytest.raises(Exception):
        bbtools.run_subprocess('false')


def test_bbtools_genome_size_haploid_and_diploid(tmp_path):
    peaks_file = tmp_path / 'peaks.txt'
    peaks_file.write_text('#haploid_genome_size 12345\n#genome_size 24690\n', encoding='utf-8')
    assert bbtools.genome_size(str(peaks_file), haploid=True) == 12345
    assert bbtools.genome_size(str(peaks_file), haploid=False) == 24690


def test_mash_sketch_dist_require_args():
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_run_subprocess(command):
        return 'stdout', 'stderr'
    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)

    out, err, cmd = bbtools.bbnorm(
        forward_in=str(forward),
        reverse_in=str(reverse),
        forward_out=str(tmp_path / 'bbnorm_R1.fastq'),
        reverse_out=str(tmp_path / 'bbnorm_R2.fastq'),
        returncmd=True
    )
    assert 'bbnorm.sh in1=' in cmd
    assert 'out2=' in cmd

    out, err, cmd = bbtools.bbduk_filter(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'filter_R1.fastq'),
        reverse_out=str(tmp_path / 'filter_R2.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh qin=33' in cmd
    assert 'out2=' in cmd

    out, err, cmd = bbtools.seal(
        reference='ref.fasta',
        forward_in=str(forward),
        reverse_in=str(reverse),
        output_file=str(tmp_path / 'seal.rpkm'),
        returncmd=True
    )
    assert 'seal.sh' in cmd
    assert 'in2=' in cmd

    out, err, cmd = bbtools.kmercountexact(
        forward_in=str(forward),
        reverse_in=str(reverse),
        returncmd=True
    )
    assert 'kmercountexact.sh in=' in cmd
    assert 'in2=' in cmd

    out, err, cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'reformatted.fastq'),
        reverse_in=str(reverse),
        reverse_out=str(tmp_path / 'reformatted_R2.fastq'),
        returncmd=True
    )
    assert 'reformat.sh in1=' in cmd

    out, err, cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'subsampled_R1.fastq'),
        reverse_out=str(tmp_path / 'subsampled_R2.fastq'),
        num_bases=100,
        returncmd=True
    )
    assert 'reformat.sh in1=' in cmd


def test_build_contig_chunks_splits_and_balances():
    contig_lengths = {'a': 5000, 'b': 15000}
    chunks = methods._build_contig_chunks(
        contig_lengths=contig_lengths,
        contigs=['a', 'b', 'c'],
        threads=2,
        multiplier=1,
        max_chunk_bases=10000
    )

    total_length = sum(
        item['length'] for chunk in chunks for item in chunk
    )
    assert total_length == 5000 + 15000 + 1
    assert all(item['length'] <= 10000 for chunk in chunks for item in chunk)
    assert len(chunks) == 2
    assert any(item['contig'] == 'b' and item['start'] is not None for chunk in chunks for item in chunk)


def test_read_contig_dispatch_invalid_kwargs():
    with pytest.raises(AttributeError):
        methods._read_contig_dispatch(kwargs='not a dict')
    with pytest.raises(TypeError):
        methods._read_contig_dispatch(kwargs={})


def test_read_contig_dispatch_raises_type_error_for_non_dict_like_kwargs():
    class DictLike:
        def get(self, key):
            return None

    with pytest.raises(TypeError):
        methods._read_contig_dispatch(kwargs=DictLike())


def test_read_contig_dispatch_calls_read_contig(monkeypatch):
    def fake_read_contig(**kwargs):
        assert kwargs['contig_name'] == 'contig1'
        return ({'contig1': {1: {}}}, 'contig1\t1\tA\n')

    monkeypatch.setattr(methods, 'read_contig', fake_read_contig)

    combined, text = methods._read_contig_dispatch(
        kwargs={
            'contig_name': 'contig1',
            'bamfile_name': 'dummy.bam',
            'reference_fasta': 'dummy.fasta',
            'fastq_records': {},
            'quality_cutoff': 20,
            'base_cutoff': 1,
            'base_fraction_cutoff': 0.1,
            'fasta': False,
            'error_cutoff': 1.0,
            'nanopore': False,
            'max_expected_positions': 0.001
        }
    )

    assert combined == {'contig1': {1: {}}}
    assert text.strip() == 'contig1\t1\tA'


def test_find_contamination_subreplicate_aggregate_report(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_database = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    def fake_pysam_faidx(path):
        Path(path + '.fai').write_text(f'{sample_database.stem}\t4\t0\t0\t0\n', encoding='utf-8')

    def fake_pysam_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        Path(output).write_bytes(b'')

    def fake_pysam_index(path):
        Path(path + '.bai').write_bytes(b'')

    class FakeAlignmentFile:
        def __init__(self, *args, **kwargs):
            self.references = [b'abc_1']
        def close(self):
            pass

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', fake_pysam_faidx)
    monkeypatch.setattr(methods.pysam, 'sort', fake_pysam_sort)
    monkeypatch.setattr(methods.pysam, 'index', fake_pysam_index)
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)

    class FakePool:
        def __init__(self, processes):
            pass
        def map(self, func, iterable, chunksize=1):
            return [({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n')]
        def close(self):
            pass
        def join(self):
            pass

    monkeypatch.setattr(methods, 'ThreadPool', FakePool)
    monkeypatch.setattr(methods, 'downsample_reads', lambda *args, **kwargs: [str(forward), str(reverse)])

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False,
        downsample_depth=20,
        subreplicates=2,
        subreplicate_seed=1
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_find_contamination_downsampled_subreplicates_full_pipeline(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_db = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_db.write_text('>abc_1\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    def fake_pysam_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        if output:
            Path(output).write_bytes(b'')

    def fake_pysam_index(path):
        Path(path + '.bai').write_bytes(b'')

    def fake_faidx(path):
        Path(path + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods, 'estimate_mean_read_length', lambda fastq_path: 1)
    monkeypatch.setattr(methods, 'count_fastq_reads', lambda fastq_path: 100_000_000)
    monkeypatch.setattr(methods, 'estimate_genome_size', lambda genus: 4_000_000)
    monkeypatch.setattr(methods.pysam, 'faidx', fake_faidx)
    monkeypatch.setattr(methods.pysam, 'sort', fake_pysam_sort)
    monkeypatch.setattr(methods.pysam, 'index', fake_pysam_index)
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', lambda *args, **kwargs: SimpleNamespace(references=[b'abc_1'], close=lambda: None))
    monkeypatch.setattr(methods, '_read_contig_chunk_dispatch', lambda kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))
    monkeypatch.setattr(methods, 'downsample_reads', lambda *args, **kwargs: [str(forward), str(reverse)])

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False,
        downsample_depth=20,
        subreplicates=2,
        subreplicate_seed=1
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_find_contamination_reports_bait_failure_when_files_missing(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    (db_folder / 'Escherichia_db_cgderived.fasta').write_text(
        '>e\nACGT\n',
        encoding='utf-8'
    )

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', lambda *args, **kwargs: ('out', 'err', 'cmd'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False
    )

    report = output_folder / 'confindr_report.tsv'
    assert report.exists()
    assert 'Error processing sample' in report.read_text()


def test_find_contamination_single_run_writes_reports(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_db = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_db.write_text('>abc_1\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    def fake_pysam_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        if output:
            Path(output).write_bytes(b'')

    def fake_pysam_index(path):
        Path(path + '.bai').write_bytes(b'')

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: Path(path + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8'))
    monkeypatch.setattr(methods.pysam, 'sort', fake_pysam_sort)
    monkeypatch.setattr(methods.pysam, 'index', fake_pysam_index)
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', lambda *args, **kwargs: SimpleNamespace(references=[b'abc_1'], close=lambda: None))
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_find_contamination_builds_rmlst_database_when_use_rmlst_true(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    (db_folder / 'gene_allele.txt').write_text(
        'Escherichia:abc_1,\n',
        encoding='utf-8'
    )
    (db_folder / 'rMLST_combined.fasta').write_text(
        '>abc_1\nACGT\n',
        encoding='utf-8'
    )

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()
    temp_db_dir = tmp_path / 'tempdb'

    sample_db_path = temp_db_dir / 'Escherichia_db.fasta'
    called = {'built': False}

    def fake_setup_allelespecific_database(fasta_file, database_folder, allele_list):
        called['built'] = True
        assert fasta_file == str(sample_db_path)
        Path(fasta_file).write_text('>abc_1\nACGT\n', encoding='utf-8')

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    monkeypatch.setattr(methods, 'find_genus_specific_allele_list', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'setup_allelespecific_database', fake_setup_allelespecific_database)
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', lambda *args, **kwargs: ('out', 'err', 'cmd'))
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(temp_db_dir),
        threads=1,
        min_matching_hashes=40,
        keep_files=False,
        use_rmlst=True
    )

    assert called['built']
    assert sample_db_path.exists()
    assert 'Error processing sample' in (output_folder / 'confindr_report.tsv').read_text()


def test_find_contamination_creates_rmlst_db_when_cgderived_missing(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    (db_folder / 'gene_allele.txt').write_text(
        'Escherichia:abc_1,\n', encoding='utf-8'
    )
    (db_folder / 'rMLST_combined.fasta').write_text(
        '>abc_1\nACGT\n', encoding='utf-8'
    )

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()
    temp_db_dir = tmp_path / 'tempdb'

    sample_db_path = temp_db_dir / 'Escherichia_db.fasta'
    called = {'setup': False}

    def fake_setup_allelespecific_database(fasta_file, database_folder, allele_list):
        called['setup'] = True
        assert fasta_file == str(sample_db_path)
        assert database_folder == str(db_folder)
        assert allele_list == ['abc_1']
        Path(fasta_file).write_text('>abc_1\nACGT\n', encoding='utf-8')

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out', 'err', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return 'out2', 'err2', 'cmd2'

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    def fake_pysam_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        if output:
            Path(output).write_bytes(b'')

    def fake_pysam_index(path):
        Path(path + '.bai').write_bytes(b'')

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    monkeypatch.setattr(methods, 'find_genus_specific_allele_list', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'setup_allelespecific_database', fake_setup_allelespecific_database)
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: Path(path + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8'))
    monkeypatch.setattr(methods.pysam, 'sort', fake_pysam_sort)
    monkeypatch.setattr(methods.pysam, 'index', fake_pysam_index)
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', lambda *args, **kwargs: SimpleNamespace(references=[b'abc_1'], close=lambda: None))
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(temp_db_dir),
        threads=1,
        min_matching_hashes=40,
        keep_files=False
    )

    assert called['setup']
    assert sample_db_path.exists()
    assert (output_folder / 'confindr_report.tsv').exists()


def test_read_contig_chunk_dispatch_merges_results(monkeypatch):
    chunk = [{'contig': 'contig1', 'start': None, 'end': None, 'length': 4}]

    def fake_read_contig(**kwargs):
        return ({'contig1': {1: {}, '_gene_stats': {}}}, 'contig1\t1\tA\n')

    monkeypatch.setattr(methods, 'read_contig', fake_read_contig)

    combined, text = methods._read_contig_chunk_dispatch({
        'contig_chunk': chunk,
        'bamfile_name': 'dummy.bam',
        'reference_fasta': 'dummy.fasta',
        'fastq_records': {},
        'quality_cutoff': 20,
        'base_cutoff': 1,
        'base_fraction_cutoff': 0.1,
        'fasta': False,
        'error_cutoff': 1.0,
        'nanopore': False,
        'max_expected_positions': 0.001
    })

    assert 'contig1' in combined
    assert text.strip() == 'contig1\t1\tA'


def test_get_fastq_record_with_unpaired_index():
    original_state = methods._FASTQ_INDEX_STATE.copy()
    methods._FASTQ_INDEX_STATE.update({'fwd': {'read1': 'FWD'}, 'rev': None, 'paired': False})
    try:
        assert methods._get_fastq_record('read1') == 'FWD'
        assert methods._get_fastq_record('read1/1') is None
        assert methods._get_fastq_record('read2') is None
    finally:
        methods._FASTQ_INDEX_STATE.clear()
        methods._FASTQ_INDEX_STATE.update(original_state)


def test_find_contamination_downsample_skips_when_reads_low(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    (db_folder / 'Escherichia_db_cgderived.fasta').write_text('>e\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    content = '@r1/1\nACGT\n+\nIIII\n'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write(content)
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write(content)

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', lambda *args, **kwargs: ('out', 'err', 'cmd'))
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'count_fastq_reads', lambda fastq_path: 1)
    monkeypatch.setattr(methods, 'estimate_mean_read_length', lambda fastq_path: 100)

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False,
        downsample_depth=20
    )

    assert (output_folder / 'confindr_report.tsv').exists()


def test_bbtools_bbmap_explicit_reverse_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    out_bam = tmp_path / 'out.bam'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_run_subprocess(command):
        return 'stdout', 'stderr'

    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)
    out, err, cmd = bbtools.bbmap(
        reference='ref.fasta',
        forward_in=str(forward),
        reverse_in=str(reverse),
        out_bam=str(out_bam),
        returncmd=True
    )

    assert 'in2=' in cmd
    assert str(reverse) in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_bbduk_trim_requires_reverse_out_with_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools.subprocess, 'check_output', lambda *args, **kwargs: b'/usr/bin/bbduk.sh\n')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))

    with pytest.raises(ValueError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'trimmed_R1.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_subsample_reads_explicit_reverse_command(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_run_subprocess(command):
        return 'stdout', 'stderr'

    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)
    out, err, cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'subsampled_R1.fastq'),
        reverse_in=str(reverse),
        reverse_out=str(tmp_path / 'subsampled_R2.fastq'),
        num_bases=100,
        returncmd=True
    )

    assert 'reformat.sh in1=' in cmd
    assert 'in2=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_validate_reads_autodetects_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.validate_reads(
        forward_in=str(forward),
        returncmd=True
    )

    assert 'reformat.sh in1=' in cmd
    assert 'stdout' == out
    assert 'stderr' == err


def test_bbtools_bbmap_single_end(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbmap(
        reference='ref.fasta',
        forward_in=str(forward),
        out_bam=str(tmp_path / 'out.bam'),
        returncmd=True
    )
    assert 'bbmap.sh ref=ref.fasta in=' in cmd
    assert 'out=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_tadpole_requires_reverse_out_with_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    with pytest.raises(ValueError):
        bbtools.tadpole(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'tadpole_R1.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_bbnorm_single_end(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbnorm(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'bbnorm.fastq'),
        returncmd=True
    )
    assert 'bbnorm.sh in=' in cmd
    assert 'out=' in cmd


def test_bbtools_bbnorm_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbnorm(
        forward_in=str(forward),
        reverse_in=str(reverse),
        forward_out=str(tmp_path / 'bbnorm_R1.fastq'),
        reverse_out=str(tmp_path / 'bbnorm_R2.fastq'),
        returncmd=True
    )
    assert 'bbnorm.sh in1=' in cmd
    assert 'out2=' in cmd


def test_bbtools_bbmerge_single_end(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbmerge(
        forward_in=str(forward),
        merged_reads=str(tmp_path / 'merged.fastq'),
        returncmd=True
    )
    assert 'bbmerge.sh in=' in cmd


def test_bbtools_bbmerge_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbmerge(
        forward_in=str(forward),
        merged_reads=str(tmp_path / 'merged.fastq'),
        reverse_in=str(reverse),
        returncmd=True
    )
    assert 'bbmerge.sh in=' in cmd
    assert 'in2=' in cmd


def test_bbtools_seal_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.seal(
        reference='ref.fasta',
        forward_in=str(forward),
        reverse_in=str(reverse),
        output_file=str(tmp_path / 'seal.rpkm'),
        returncmd=True
    )
    assert 'seal.sh ref=ref.fasta in=' in cmd
    assert 'in2=' in cmd


def test_bbtools_kmercountexact_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.kmercountexact(
        forward_in=str(forward),
        reverse_in=str(reverse),
        returncmd=True
    )
    assert 'kmercountexact.sh in=' in cmd
    assert 'in2=' in cmd


def test_bbtools_subsample_reads_output_already_exists(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output_file = tmp_path / 'subsampled_R1.fastq'
    output_file.write_text('', encoding='utf-8')

    def fake_run_subprocess(command):
        raise AssertionError('should not run subprocess')

    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)
    out, err, cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(output_file),
        num_bases=100,
        returncmd=True
    )
    assert out == ''
    assert err == ''


def test_bbtools_reformat_reads_output_already_exists(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output_file = tmp_path / 'reformatted_R1.fastq'
    output_file.write_text('', encoding='utf-8')

    def fake_run_subprocess(command):
        raise AssertionError('should not run subprocess')

    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)
    out, err, cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(output_file),
        returncmd=True
    )
    assert out == ''
    assert err == ''


def test_bbtools_repair_reads_output_already_exists(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    output_file = tmp_path / 'repair_R1.fastq'
    output_file.write_text('', encoding='utf-8')

    def fake_run_subprocess(command):
        raise AssertionError('should not run subprocess')

    monkeypatch.setattr(bbtools, 'run_subprocess', fake_run_subprocess)
    out, err, cmd = bbtools.repair_reads(
        forward_in=str(forward),
        forward_out=str(output_file),
        reverse_in=str(reverse),
        reverse_out=str(tmp_path / 'repair_R2.fastq'),
        returncmd=True
    )
    assert out == ''
    assert err == ''


def test_bbtools_validate_reads_single_end(monkeypatch, tmp_path):
    forward = tmp_path / 'reads.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.validate_reads(
        forward_in=str(forward),
        returncmd=True
    )
    assert 'reformat.sh in=' in cmd


def test_bbtools_run_subprocess_failure(monkeypatch):
    class FakeProcess:
        def __init__(self, command, shell, stdout, stderr):
            self.returncode = 1
        def communicate(self):
            return b'', b'err'

    monkeypatch.setattr(bbtools, 'Popen', FakeProcess)
    with pytest.raises(Exception):
        bbtools.run_subprocess('false')


def test_bbtools_genome_size_diploid(tmp_path):
    peaks_file = tmp_path / 'peaks.txt'
    peaks_file.write_text('#genome_size 12345\n', encoding='utf-8')
    assert bbtools.genome_size(str(peaks_file), haploid=False) == 12345


def test_mash_sketch_dist_require_args():
    with pytest.raises(ValueError):
        mash.sketch()
    with pytest.raises(ValueError):
        mash.dist()
    out, err, cmd = mash.screen('a.msh', 'b.msh', returncmd=True)
    assert 'mash screen a.msh b.msh' in cmd

    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbduk_filter(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'filter_R1.fastq'),
        returncmd=True
    )

    assert 'in2=' in cmd
    assert 'out2=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_bbduk_filter_requires_reverse_out_with_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    with pytest.raises(ValueError):
        bbtools.bbduk_filter(
            reference='ref.fasta',
            forward_in=str(forward),
            forward_out=str(tmp_path / 'filter_R1.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_bbduk_trim_raises_file_not_found_when_bbduk_missing(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    def fake_check_output(*args, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=args)

    monkeypatch.setattr(bbtools.subprocess, 'check_output', fake_check_output)
    with pytest.raises(FileNotFoundError):
        bbtools.bbduk_trim(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'trimmed_R1.fastq')
        )


def test_bbtools_tadpole_explicit_reverse_out(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.tadpole(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'tadpole_R1.fastq'),
        returncmd=True
    )
    assert 'tadpole.sh' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_bbnorm_explicit_reverse_out(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbnorm(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'bbnorm_R1.fastq'),
        reverse_in=str(reverse),
        reverse_out=str(tmp_path / 'bbnorm_R2.fastq'),
        returncmd=True
    )
    assert 'bbnorm.sh in1=' in cmd
    assert 'out2=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_bbmerge_explicit_reverse_out(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbmerge(
        forward_in=str(forward),
        merged_reads=str(tmp_path / 'merged.fastq'),
        returncmd=True
    )
    assert 'bbmerge.sh in=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_bbduk_bait_autodetects_pair_and_reverse_out(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.bbduk_bait(
        reference='ref.fasta',
        forward_in=str(forward),
        forward_out=str(tmp_path / 'baited_R1.fastq'),
        returncmd=True
    )
    assert 'bbduk.sh qin=33' in cmd
    assert 'outm2=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_subsample_reads_no_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.subsample_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'subsampled_R1.fastq'),
        num_bases=100,
        returncmd=True
    )
    assert 'reformat.sh in=' in cmd
    assert 'out=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_reformat_reads_no_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    out, err, cmd = bbtools.reformat_reads(
        forward_in=str(forward),
        forward_out=str(tmp_path / 'reformatted_R1.fastq'),
        returncmd=True
    )
    assert 'reformat.sh in=' in cmd
    assert 'out=' in cmd
    assert out == 'stdout'
    assert err == 'stderr'


def test_bbtools_repair_reads_no_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    with pytest.raises(ValueError):
        bbtools.repair_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'repair_R1.fastq'),
            reverse_in='NA',
            reverse_out='NA'
        )


def test_bbtools_genome_size_haploid_missing_field(tmp_path):
    peaks_file = tmp_path / 'peaks.txt'
    peaks_file.write_text('#foo 12345\n', encoding='utf-8')
    assert bbtools.genome_size(str(peaks_file), haploid=True) == 0


def test_mash_sketch_dist_require_args():
    with pytest.raises(ValueError):
        mash.sketch()
    with pytest.raises(ValueError):
        mash.dist()
    out, err, cmd = mash.screen('a.msh', 'b.msh', returncmd=True)
    assert 'mash screen a.msh b.msh' in cmd


def test_find_genus_specific_allele_list_returns_only_target_genus(tmp_path):
    profiles = tmp_path / 'profiles.txt'
    profiles.write_text('Escherichia:BACT000001_1,BACT000002_2,\nShigella:BACT000001_3,\n', encoding='utf-8')
    alleles = methods.find_genus_specific_allele_list(
        profiles_file=str(profiles),
        target_genus='Escherichia'
    )
    assert alleles == ['BACT000001_1', 'BACT000002_2']


def test_setup_allelespecific_database_writes_selected_alleles(monkeypatch, tmp_path):
    db_dir = tmp_path / 'db'
    db_dir.mkdir()
    combined = db_dir / 'rMLST_combined.fasta'
    combined.write_text('>BACT000001_1\nACGT\n>BACT000002_2\nTGCA\n', encoding='utf-8')
    methods.setup_allelespecific_database(
        fasta_file=str(tmp_path / 'subset.fasta'),
        database_folder=str(db_dir),
        allele_list=['BACT000001_1']
    )
    written = Path(tmp_path / 'subset.fasta').read_text()
    assert '>BACT000001_1' in written
    assert '>BACT000002_2' not in written


def test_setup_allelespecific_database_ignores_write_file_not_found(monkeypatch, tmp_path):
    def fake_index(path, fmt):
        return {'abc_1': 'record'}

    def fake_write(seqs, fasta_file, fmt):
        raise FileNotFoundError('no such path')

    monkeypatch.setattr(methods.SeqIO, 'index', lambda path, fmt: {'abc_1': 'record'})
    monkeypatch.setattr(methods.SeqIO, 'write', fake_write)

    methods.setup_allelespecific_database(
        fasta_file=str(tmp_path / 'nonexistent_dir/subset.fasta'),
        database_folder=str(tmp_path),
        allele_list=['abc_1']
    )


def test_find_cross_contamination_converts_shigella_and_formats_multiple(monkeypatch, tmp_path):
    tmpdir = tmp_path / 'tmp'
    tmpdir.mkdir()
    screen_file = tmpdir / 'sample_screen.tab'
    screen_file.write_text('99.9 100/100 1 1e-5 species/Escherichia/strain\n99.8 50/100 1 1e-5 species/Shigella/strain\n', encoding='utf-8')

    class FakeScreenResult:
        def __init__(self, row):
            parts = row.split()
            self.query_id = parts[4]
            self.shared_hashes = parts[1]

    monkeypatch.setattr(methods.mash, 'screen', lambda *args, **kwargs: ('', '', ''))
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods.mash, 'read_mash_screen', lambda path: [
        FakeScreenResult('99.9 100/100 1 1e-5 foo/Escherichia/strain/x'),
        FakeScreenResult('99.8 50/100 1 1e-5 foo/Shigella/strain/x')
    ])
    genus = methods.find_cross_contamination(
        databases=str(tmp_path),
        reads=[str(tmp_path / 'reads_R1.fastq'), str(tmp_path / 'reads_R2.fastq')],
        sample_name='sample',
        tmpdir=str(tmpdir),
        log=str(tmpdir / 'log.txt'),
        threads=1,
        min_matching_hashes=20
    )
    assert genus == 'Escherichia'


def test_find_cross_contamination_returns_ND_when_no_matching_hashes(monkeypatch, tmp_path):
    tmpdir = tmp_path / 'tmp'
    tmpdir.mkdir()
    monkeypatch.setattr(methods.mash, 'screen', lambda *args, **kwargs: ('', '', ''))
    class FakeScreenResult:
        def __init__(self, row):
            parts = row.split()
            self.query_id = parts[4]
            self.shared_hashes = parts[1]
    monkeypatch.setattr(methods.mash, 'read_mash_screen', lambda path: [
        FakeScreenResult('99.9 10/100 1 1e-5 species/Escherichia/strain')
    ])
    genus = methods.find_cross_contamination(
        databases=str(tmp_path),
        reads=str(tmp_path / 'reads.fastq'),
        sample_name='sample',
        tmpdir=str(tmpdir),
        log=str(tmpdir / 'log.txt'),
        threads=1,
        min_matching_hashes=20
    )
    assert genus == 'ND'


def test_number_of_bases_above_threshold_with_fraction():
    result = methods.number_of_bases_above_threshold(
        high_quality_base_count={'A': 10, 'C': 2, 'G': 1},
        base_count_cutoff=2,
        base_fraction_cutoff=0.2
    )
    assert result == 1


def test_parse_bam_converts_bytes_contig_and_start_end(monkeypatch):
    class FakeAlignmentFile:
        def __init__(self, name, mode):
            self.name = name
        def pileup(self, contig, start=None, end=None, **kwargs):
            assert contig == 'chr1'
            assert start == 10
            assert end == 20
            return []
    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    class FakeFasta:
        pass
    bamfile, pileup = methods.parse_bam(
        bamfile_name='dummy.bam',
        contig_name=b'chr1',
        pysam_fasta=FakeFasta(),
        start=10,
        end=20
    )
    assert pileup == []


def test_position_details_formats_statistics_correctly():
    formatted = methods.position_details(
        actual_position=5,
        passing_snv_dict={
            'congruent': {'A': 1, 'C': 0, 'G': 0, 'T': 0},
            'paired': {'A': 0, 'C': 1, 'G': 0, 'T': 0},
            'forward': {'A': 0, 'C': 0, 'G': 1, 'T': 0},
            'reverse': {'A': 0, 'C': 0, 'G': 0, 'T': 1}
        },
        contig_name='chr1',
        ref_base='A',
        total_coverage=4,
        base_cutoff=2,
        error_perc=0.1234,
        p_value=1e-4,
        adj_p_value=2e-4,
        strand_p=3e-4,
        pos_p=4e-4,
        mean_q=25.5,
        mean_mapq=40.75
    )
    assert 'chr1\t5\tA' in formatted
    assert '0.12' in formatted
    assert '1.000e-04' in formatted
    assert '25.50' in formatted


def test_dependency_check_uses_shutil_which(monkeypatch):
    monkeypatch.setattr(methods.shutil, 'which', lambda dep: '/usr/bin/' + dep)
    assert methods.dependency_check(dependency='bash') is True


def test_find_paired_and_unpaired_reads_detect_missing_reverse(tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    assert methods.find_paired_reads(fastq_directory=str(tmp_path)) == []
    assert methods.find_unpaired_reads(fastq_directory=str(tmp_path)) == [[str(forward)]]


def test_write_output_probabilistic_decision_true(tmp_path):
    out = tmp_path / 'confindr_report.tsv'
    methods.write_output(
        output_report=str(out),
        sample_name='TestSample',
        multi_positions=0,
        genus='Fakella',
        total_gene_length=1000,
        database_download_date='ND',
        use_probabilistic=True,
        score_threshold=1.0,
        sample_score=1.0
    )
    assert 'TestSample' in out.read_text()
    assert 'ND' in out.read_text()


def test_bbtools_subsample_reads_explicit_reverse_requires_reverse_out(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    with pytest.raises(ValueError):
        bbtools.subsample_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'subsampled_R1.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA',
            num_bases=100
        )


def test_bbtools_reformat_reads_requires_reverse_out_with_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    with pytest.raises(ValueError):
        bbtools.reformat_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'reformatted_R1.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_repair_reads_requires_reverse_out_with_explicit_reverse(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    monkeypatch.setattr(bbtools, 'run_subprocess', lambda command: ('stdout', 'stderr'))
    with pytest.raises(ValueError):
        bbtools.repair_reads(
            forward_in=str(forward),
            forward_out=str(tmp_path / 'repair_R1.fastq'),
            reverse_in=str(reverse),
            reverse_out='NA'
        )


def test_bbtools_genome_size_diploid(tmp_path):
    peaks_file = tmp_path / 'peaks.txt'
    peaks_file.write_text('#genome_size 12345\n', encoding='utf-8')
    assert bbtools.genome_size(str(peaks_file), haploid=False) == 12345


def test_write_output_subreplicates_summary_outputs_mean_median_stddev(tmp_path):
    out = tmp_path / 'confindr_report.tsv'
    methods.write_output(
        output_report=str(out),
        sample_name='TestSample',
        multi_positions=3,
        genus='Fakella',
        total_gene_length=1000,
        database_download_date='ND',
        use_probabilistic=True,
        score_threshold=2.0,
        sample_score=4.0,
        subreplicate_counts=[1, 2, 3]
    )
    text = out.read_text()
    assert 'MeanContamSNVs' in text
    assert 'Subreplicate_1' in text
    assert 'TestSample' in text


def test_read_contig_chunk_dispatch_raises_on_read_contig_failure(monkeypatch):
    def fake_read_contig(**kwargs):
        raise ValueError('boom')
    monkeypatch.setattr(methods, 'read_contig', fake_read_contig)
    with pytest.raises(ValueError):
        methods._read_contig_chunk_dispatch(
            {'contig_chunk': [{'contig': 'contig1', 'start': None, 'end': None}]}
        )


def test_ensure_fastq_index_handles_index_failure(monkeypatch, tmp_path):
    forward = tmp_path / 'reads_R1.fastq'
    reverse = tmp_path / 'reads_R2.fastq'
    forward.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')
    reverse.write_text('@r1\nTGCA\n+\nIIII\n', encoding='utf-8')

    def fake_index_db(path, fasta, fmt):
        raise OSError('index failure')

    monkeypatch.setattr(methods.SeqIO, 'index_db', fake_index_db)
    methods._ensure_fastq_index(
        fwd_path=str(forward),
        rev_path=str(reverse),
        paired=True,
        tmpdir=str(tmp_path)
    )
    assert not (tmp_path / 'reads_R1.fastq.fwd.idx').exists()
    assert not (tmp_path / 'reads_R2.fastq.rev.idx').exists()


def test_find_contamination_nonreplicate_nanopore_fasta_runs_through_with_fake_pool(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    (db_folder / 'Escherichia_db_cgderived.fasta').write_text('>abc_1\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return '', '', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_bytes(b'')
        if reverse_out:
            Path(reverse_out).write_bytes(b'')
        return '', '', 'cmd'

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'read_contig', lambda **kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    def fake_run_cmd(cmd):
        return '', ''

    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)

    def fake_faidx(path):
        Path(path + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    monkeypatch.setattr(methods.pysam, 'faidx', fake_faidx)
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(path + '.bai').write_bytes(b''))

    monkeypatch.setattr(
        methods,
        '_read_contig_chunk_dispatch',
        lambda kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n')
    )

    class FakePool:
        def __init__(self, processes):
            pass
        def map(self, func, iterable, chunksize=1):
            return [({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n')]
        def close(self):
            pass
        def join(self):
            pass

    monkeypatch.setattr(methods, 'ThreadPool', FakePool)

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=1,
        min_matching_hashes=40,
        keep_files=False,
        data_type='Nanopore',
        fasta=True
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test__build_contig_chunks_splits_and_balances():
    chunks = methods._build_contig_chunks(
        contig_lengths={'a': 600, 'b': 300, 'c': 100},
        contigs=['a', 'b', 'c'],
        threads=1,
        multiplier=2,
        max_chunk_bases=500
    )
    assert len(chunks) == 2
    assert sum(item['length'] for item in chunks[0] + chunks[1]) == 1000


def test__read_contig_chunk_dispatch_merges_reports(monkeypatch):
    def fake_read_contig(**kwargs):
        return ({kwargs['contig_name']: {1: {}}}, f"{kwargs['contig_name']}\t1\tA\n")
    monkeypatch.setattr(methods, 'read_contig', fake_read_contig)

    combined, report = methods._read_contig_chunk_dispatch({
        'contig_chunk': [
            {'contig': 'abc_1', 'start': None, 'end': None, 'length': 4}
        ],
        'bamfile_name': 'dummy.bam',
        'reference_fasta': 'dummy.fasta',
        'quality_cutoff': 20,
        'min_quality': 15,
        'base_cutoff': 1,
        'base_fraction_cutoff': 0.1,
        'fasta': False,
        'error_cutoff': 1.0,
        'nanopore': False,
        'max_expected_positions': 0.001
    })
    assert combined == {'abc_1': {1: {}}}
    assert report.strip() == 'abc_1\t1\tA'


def test__read_contig_dispatch_invalid_kwargs_raises_type_error():
    with pytest.raises(TypeError):
        methods._read_contig_dispatch(kwargs={})


def test_count_multibase_positions_excludes_meta_keys():
    data = [
        {'g1': {1: {}, '_gene_stats': {}}},
        {'g2': {2: {}, 3: {}, '_gene_stats': {}}}
    ]
    assert methods.count_multibase_positions(multibase_dict_list=data) == 3


def test_find_rmlst_type_writes_sorted_report(tmp_path):
    kma_report = tmp_path / 'kma.res'
    rmlst_report = tmp_path / 'alleles.tsv'
    kma_report.write_text(
        '#Template\tScore\n'
        'abc_1\t10\n'
        'abc_2\t20\n'
        'def_1\t5\n',
        encoding='utf-8'
    )
    alleles = methods.find_rmlst_type(
        kma_report=str(kma_report),
        rmlst_report=str(rmlst_report)
    )
    assert alleles == ['abc_2', 'def_1']
    text = rmlst_report.read_text(encoding='utf-8').splitlines()
    assert text[0] == 'Gene\tAllele'
    assert text[1] == 'abc\t2'


def test__valid_downsample_depth_accepts_and_rejects():
    assert methods._valid_downsample_depth('10') == 10
    assert methods._valid_downsample_depth('100') == 100
    with pytest.raises(argparse.ArgumentTypeError):
        methods._valid_downsample_depth('9')
    with pytest.raises(argparse.ArgumentTypeError):
        methods._valid_downsample_depth('bad')


def test_recommend_xmx_uses_available_memory(monkeypatch):
    monkeypatch.setattr(methods.psutil, 'virtual_memory', lambda: SimpleNamespace(available=5 * 1024**3))
    assert methods.recommend_xmx(fraction=0.5) == '2g'


def test_blast_result_parses_fields():
    line = 'query subject 80.0 50 100 1 50 1 50 1e-5'
    br = cgdb.BlastResult(line)
    assert br.query_name == 'query'
    assert br.subject_name == 'subject'
    assert br.percent_identity == 80.0
    assert br.query_coverage == 50.0


def test_get_potential_genes_selects_top_genes(tmp_path):
    report = tmp_path / 'gene_hit_report.tsv'
    report.write_text(
        'Gene\tOneHitPerGenome\n'
        'g1\t1.0\n'
        'g2\t0.2\n'
        'g3\t0.1\n',
        encoding='utf-8'
    )
    genes = cgdb.get_potential_genes(
        gene_report=str(report),
        desired_genes=2
    )
    assert genes == ['g1', 'g2']


def test_database_setup_get_loci_and_scheme_url_handles_text_and_json(monkeypatch, tmp_path):
    class FakeResponse:
        def __init__(self, status_code, headers, text, json_data=None):
            self.status_code = status_code
            self.headers = headers
            self.text = text
            self._json_data = json_data
        def json(self):
            return self._json_data

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            if url.endswith('/oauth/get_session_token'):
                return FakeResponse(200, {'content-type': 'application/json'}, '', {'oauth_token': 'tok', 'oauth_token_secret': 'sec'})
            return FakeResponse(200, {'content-type': 'text/plain'}, '{"loci": ["http://example.com/locus"], "schemes": "http://example.com/schemes"}')

    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)
    rest = dbsetup.RmlstRest(str(secret_file), str(tmp_path), unverified=False)
    rest.test_rest_url = 'http://example.com'
    rest.session_token = 'tok'
    rest.session_secret = 'sec'
    rest.get_loci_and_scheme_url()
    assert rest.loci == ['http://example.com/locus']
    assert rest.profile == 'http://example.com/schemes'


def test_find_contamination_parallel_branch(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_database = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')

    forward = tmp_path / 'reads_R1.fastq.gz'
    reverse = tmp_path / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')
    def fake_bbduk_bait(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            with gzip.open(forward_out, 'wt', encoding='utf-8') as f:
                f.write('@r1/1\nACGT\n+\nIIII\n')
        if reverse_out:
            with gzip.open(reverse_out, 'wt', encoding='utf-8') as f:
                f.write('@r1/2\nACGT\n+\nIIII\n')
        return '', '', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            with gzip.open(forward_out, 'wt', encoding='utf-8') as f:
                f.write('@r1/1\nACGT\n+\nIIII\n')
        if reverse_out:
            with gzip.open(reverse_out, 'wt', encoding='utf-8') as f:
                f.write('@r1/2\nACGT\n+\nIIII\n')
        return '', '', 'cmd'

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: sample_database.replace('.fasta', '_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'find_total_sequence_length', lambda fasta_file: 4)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)

    def fake_run_cmd(cmd):
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', lambda *args, **kwargs: None)
    def fake_pysam_sort(*args, **kwargs):
        output = args[1] if len(args) > 1 else kwargs.get('out')
        if output:
            Path(output).write_bytes(b'')
    monkeypatch.setattr(methods.pysam, 'sort', fake_pysam_sort)
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(str(path) + '.bai').write_bytes(b''))

    class FakeAlignmentFile:
        def __init__(self, path, mode):
            self.path = path
        def __enter__(self):
            return self
        def __exit__(self, exc_type, exc_val, exc_tb):
            return False
        @property
        def references(self):
            return [b'abc_1']
        def close(self):
            pass

    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), []))
    monkeypatch.setattr(methods, '_ensure_fastq_index', lambda *args, **kwargs: None)

    class FakePool:
        def __init__(self, processes, initializer=None, initargs=None):
            pass
        def map(self, func, iterable, chunksize=1):
            return [({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n')]
        def close(self):
            pass
        def join(self):
            pass

    monkeypatch.setattr(methods.multiprocessing, 'Pool', FakePool)

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path),
        threads=2,
        min_matching_hashes=40,
        keep_files=False,
        data_type='Illumina'
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()
