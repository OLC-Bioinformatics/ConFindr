#!/usr/bin/env python3

import gzip
import os
import shutil
import tarfile
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import confindr_src.create_genus_specific_db as cgdb
import confindr_src.database_setup as dbsetup
import confindr_src.methods as methods
from confindr_src.wrappers import bbtools


def test_find_unpaired_reads_reverse_only(tmp_path):
    reverse = tmp_path / 'sample_R2.fastq'
    reverse.write_text('@r1\nACGT\n+\nIIII\n', encoding='utf-8')

    result = methods.find_unpaired_reads(
        fastq_directory=str(tmp_path),
        forward_id='_R1',
        reverse_id='_R2'
    )

    assert [str(reverse)] in result


def test_find_cross_contamination_multiple_genera(monkeypatch, tmp_path):
    def fake_screen(*args, **kwargs):
        output_file = kwargs.get('output_file')
        Path(output_file).write_text('', encoding='utf-8')
        return '', '', 'mash screen'

    monkeypatch.setattr(methods.mash, 'screen', fake_screen)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(
        methods.mash,
        'read_mash_screen',
        lambda path: [
            SimpleNamespace(query_id='a/b/x/Escherichia/c/d', shared_hashes='40/100'),
            SimpleNamespace(query_id='a/b/x/Bacillus/c/d', shared_hashes='40/100'),
            SimpleNamespace(query_id='a/b/x/Shigella/c/d', shared_hashes='40/100')
        ]
    )

    result = methods.find_cross_contamination(
        databases=str(tmp_path),
        reads=str(tmp_path / 'reads.fastq'),
        sample_name='sample',
        tmpdir=str(tmp_path),
        log=str(tmp_path / 'log.txt'),
        threads=1,
        min_matching_hashes=40
    )

    assert result == 'Escherichia:Bacillus'


def test_download_mash_sketch_uses_urlretrieve(monkeypatch, tmp_path):
    called = []

    def fake_urlretrieve(url, dest):
        called.append((url, dest))
        Path(dest).write_text('dummy', encoding='utf-8')
        return dest, None

    monkeypatch.setattr(methods.urllib.request, 'urlretrieve', fake_urlretrieve)
    methods.download_mash_sketch(output_folder=str(tmp_path))

    assert called
    assert called[0][0].endswith('refseq.msh')
    assert (tmp_path / 'refseq.msh').exists()


def test_download_cgmlst_derived_data_extracts_tarball_and_indexes(monkeypatch, tmp_path):
    tarball = tmp_path / 'source.tar.gz'
    payload = tmp_path / 'file.txt'
    payload.write_text('content', encoding='utf-8')
    with tarfile.open(tarball, 'w:gz') as tar:
        tar.add(str(payload), arcname='file.txt')

    def fake_urlretrieve(url, dest):
        shutil.copy(str(tarball), dest)
        return dest, None

    called = {'indexed': False}

    def fake_index(**kwargs):
        called['indexed'] = True
        assert kwargs['output_folder'] == str(tmp_path)
        assert kwargs['cgderived'] is True

    monkeypatch.setattr(methods.urllib.request, 'urlretrieve', fake_urlretrieve)
    monkeypatch.setattr(methods, 'index', fake_index)

    methods.download_cgmlst_derived_data(output_folder=str(tmp_path))

    assert called['indexed'] is True
    assert (tmp_path / 'file.txt').exists()
    assert not (tmp_path / 'confindr_db.tar.gz').exists()


def test_index_uses_existing_cgderived_database(tmp_path, monkeypatch):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_database = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')

    recorded = {}

    def fake_index_databases(*args, **kwargs):
        recorded['sample_database'] = kwargs.get('sample_database') if kwargs else args[0]

    monkeypatch.setattr(methods, 'index_databases', fake_index_databases)

    methods.index(output_folder=str(db_folder), genera=['Escherichia'], cgderived=True)

    assert recorded['sample_database'] == str(sample_database)


def test_index_builds_rmlst_genus_database_when_cgderived_missing(tmp_path, monkeypatch):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    (db_folder / 'gene_allele.txt').write_text('Escherichia:abc_1,\n', encoding='utf-8')
    (db_folder / 'rMLST_combined.fasta').write_text('>abc_1\nACGT\n', encoding='utf-8')

    calls = {'alleles': False, 'indexed': False}

    monkeypatch.setattr(
        methods,
        'find_genus_specific_allele_list',
        lambda *args, **kwargs: ['abc_1']
    )

    def fake_setup_allelespecific_database(**kwargs):
        calls['alleles'] = True
        assert kwargs['fasta_file'].endswith('Escherichia_db_cgderived.fasta')
        assert kwargs['database_folder'] == str(db_folder)
        assert kwargs['allele_list'] == ['abc_1']

    def fake_index_databases(*args, **kwargs):
        calls['indexed'] = True

    monkeypatch.setattr(methods, 'setup_allelespecific_database', fake_setup_allelespecific_database)
    monkeypatch.setattr(methods, 'index_databases', fake_index_databases)

    methods.index(output_folder=str(db_folder), genera=['Escherichia'], cgderived=True)

    assert calls['alleles'] is True
    assert calls['indexed'] is True


def test_find_contamination_uses_cgmlst_db_and_threadpool(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    cgmlst_db = tmp_path / 'cgmlst.fasta'
    cgmlst_db.write_text('>abc_1\nACGT\n', encoding='utf-8')

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
            Path(forward_out).write_text('@r1/1\nACGT\n+\nIIII\n', encoding='utf-8')
        if reverse_out:
            Path(reverse_out).write_text('@r1/2\nACGT\n+\nIIII\n', encoding='utf-8')
        return '', '', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        forward_out = kwargs.get('forward_out')
        reverse_out = kwargs.get('reverse_out')
        if forward_out:
            Path(forward_out).write_text('@r1/1\nACGT\n+\nIIII\n', encoding='utf-8')
        if reverse_out:
            Path(reverse_out).write_text('@r1/2\nACGT\n+\nIIII\n', encoding='utf-8')
        return '', '', 'cmd'

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'find_total_sequence_length', lambda fasta_file: 4)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)

    def fake_run_cmd(cmd):
        if ' -o ' in cmd:
            out_path = cmd.split(' -o ')[1].split()[0]
            if out_path.endswith('_kma'):
                Path(out_path + '.res').write_text('', encoding='utf-8')
            else:
                Path(out_path).write_bytes(b'')
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: Path(str(path) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8'))
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: Path(args[1] if len(args) > 1 else kwargs.get('out')).write_bytes(b''))
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(str(path) + '.bai').write_bytes(b''))

    class FakeAlignmentFile:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        @property
        def references(self):
            return [b'abc_1']

        def pileup(self, *args, **kwargs):
            return []

        def close(self):
            pass

    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
    monkeypatch.setattr(methods, 'load_fastq_records', lambda *args, **kwargs: {'read1/1': SimpleNamespace()} if kwargs.get('forward') else {'read1/2': SimpleNamespace()})
    monkeypatch.setattr(methods, 'ThreadPool', type('FakePool', (), {
        '__init__': lambda self, processes: None,
        'map': lambda self, func, iterable, chunksize=1: [({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n')],
        'close': lambda self: None,
        'join': lambda self: None
    }))
    monkeypatch.setattr(methods, '_read_contig_chunk_dispatch', lambda kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward), str(reverse)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        tmpdir=str(tmp_path / 'tmp'),
        threads=1,
        min_matching_hashes=40,
        keep_files=False,
        data_type='Illumina',
        cgmlst_db=str(cgmlst_db)
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_find_contamination_unpaired_illumina_uses_multiprocessing_pool(monkeypatch, tmp_path):
    db_folder = tmp_path / 'db'
    db_folder.mkdir()
    sample_database = db_folder / 'Escherichia_db_cgderived.fasta'
    sample_database.write_text('>abc_1\nACGT\n', encoding='utf-8')
    Path(str(sample_database) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8')

    forward = tmp_path / 'reads.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1\nACGT\n+\nIIII\n')

    output_folder = tmp_path / 'out'
    output_folder.mkdir()

    monkeypatch.setattr(methods, 'find_cross_contamination', lambda *args, **kwargs: 'Escherichia')

    def fake_bbduk_bait(*args, **kwargs):
        out = kwargs.get('forward_out')
        if out:
            Path(out).write_bytes(b'')
        return '', '', 'cmd'

    def fake_bbduk_trim(*args, **kwargs):
        out = kwargs.get('forward_out')
        if out:
            Path(out).write_bytes(b'')
        return '', '', 'cmd'

    def fake_run_cmd(cmd):
        if ' -o ' in cmd:
            out_path = cmd.split(' -o ')[1].split()[0]
            if out_path.endswith('_kma'):
                Path(out_path + '.res').write_text('', encoding='utf-8')
            else:
                Path(out_path).write_bytes(b'')
        if 'out=' in cmd:
            out_path = cmd.split('out=')[1].split()[0]
            Path(out_path).write_bytes(b'')
        return '', ''

    monkeypatch.setattr(methods.bbtools, 'bbduk_bait', fake_bbduk_bait)
    monkeypatch.setattr(methods.bbtools, 'bbduk_trim', fake_bbduk_trim)
    monkeypatch.setattr(methods, 'index_databases', lambda sample_database: str(tmp_path / 'sample_kma'))
    monkeypatch.setattr(methods, 'find_rmlst_type', lambda *args, **kwargs: ['abc_1'])
    monkeypatch.setattr(methods, 'find_total_sequence_length', lambda fasta_file: 4)
    monkeypatch.setattr(methods, 'write_to_logfile', lambda *args, **kwargs: None)
    monkeypatch.setattr(methods, 'run_cmd', fake_run_cmd)
    monkeypatch.setattr(methods.pysam, 'faidx', lambda path: Path(str(path) + '.fai').write_text('abc_1\t4\t0\t0\t0\n', encoding='utf-8'))
    monkeypatch.setattr(methods.pysam, 'sort', lambda *args, **kwargs: Path(args[1] if len(args) > 1 else kwargs.get('out')).write_bytes(b''))
    monkeypatch.setattr(methods.pysam, 'index', lambda path: Path(str(path) + '.bai').write_bytes(b''))

    class FakeAlignmentFile:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        @property
        def references(self):
            return [b'abc_1']

        def pileup(self, *args, **kwargs):
            return []

        def close(self):
            pass

    monkeypatch.setattr(methods.pysam, 'AlignmentFile', FakeAlignmentFile)
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
    monkeypatch.setattr(methods, '_read_contig_chunk_dispatch', lambda kwargs: ({'abc_1': {1: {}, '_gene_stats': {}}}, 'abc_1\t1\tA\n'))

    methods.find_contamination(
        pair=[str(forward)],
        output_folder=str(output_folder),
        databases_folder=str(db_folder),
        base_cutoff=3,
        xmx='4g',
        threads=2,
        min_matching_hashes=40,
        keep_files=False,
        data_type='Illumina',
        fasta=False
    )

    assert (output_folder / 'confindr_report.tsv').exists()
    assert (output_folder / 'reads_contamination.tsv').exists()
    assert (output_folder / 'reads_gene_summary.tsv').exists()


def test_read_contig_fasta_index_fallback_to_full_parse(monkeypatch, tmp_path):
    reference_fasta = tmp_path / 'ref.fasta'
    reference_fasta.write_text('>abc_1\nACGT\n', encoding='utf-8')

    class FakeColumn:
        def __init__(self):
            self.pos = 0
            self.reference_name = 'abc_1'
            self.pileups = []

    monkeypatch.setattr(methods.SeqIO, 'index', lambda path, fmt: (_ for _ in ()).throw(FileNotFoundError('fail')))
    monkeypatch.setattr(methods, 'parse_bam', lambda **kwargs: (SimpleNamespace(close=lambda: None), [FakeColumn()]))
    monkeypatch.setattr(methods, 'characterise_read', lambda **kwargs: ({'congruent_ref': {'A': 2}}, [30], {}))
    monkeypatch.setattr(methods, 'determine_cutoff', lambda **kwargs: (3, 0.0, 0.0))
    monkeypatch.setattr(methods, 'find_multibase_positions', lambda **kwargs: ({'total': 2}, {'congruent': {'A': 2}, 'forward': {}, 'reverse': {}, 'paired': {}}, 2, {}))
    monkeypatch.setattr(methods, 'benjamini_hochberg', lambda pvals: [0.01])
    monkeypatch.setattr(methods, '_position_entry_passes_probabilistic_gating', lambda **kwargs: True)
    monkeypatch.setattr(methods, 'combine_pvalues_fisher', lambda pvals: 0.05)
    monkeypatch.setattr(methods, 'position_details', lambda **kwargs: 'abc_1\t1\tA\n')

    result, report_text = methods.read_contig(
        contig_name='abc_1',
        bamfile_name='dummy.bam',
        reference_fasta=str(reference_fasta),
        allele_records={'abc_1': SimpleNamespace(seq='ACGT')},
        fastq_records={},
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=1,
        fasta=True
    )

    assert 'abc_1' in result
    assert report_text == 'abc_1\t1\tA\n'


def test_characterise_read_forward_ref_reverse_um_qf_quality_fail():
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
        return SimpleNamespace(query_position=0, alignment=alignment)

    forward = make_read('read1', 'C', True)
    reverse = make_read('read1', 'A', False)
    column = SimpleNamespace(pos=0, reference_name='gene1', pileups=[forward, reverse])
    fastq_records = {
        'read1/1': SeqRecord(Seq('C'), id='read1/1', letter_annotations={'phred_quality': [10]}),
        'read1/2': SeqRecord(Seq('A'), id='read1/2', letter_annotations={'phred_quality': [30]})
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

    assert filtered['reverse_ref_forward_UM_QF'] == {'A': 1}
    assert qualities == [30]


def test_characterise_read_reverse_ref_forward_um_qf_quality_fail():
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
        return SimpleNamespace(query_position=0, alignment=alignment)

    forward = make_read('read1', 'A', True)
    reverse = make_read('read1', 'C', False)
    column = SimpleNamespace(pos=0, reference_name='gene1', pileups=[forward, reverse])
    fastq_records = {
        'read1/1': SeqRecord(Seq('A'), id='read1/1', letter_annotations={'phred_quality': [30]}),
        'read1/2': SeqRecord(Seq('C'), id='read1/2', letter_annotations={'phred_quality': [10]})
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

    assert filtered['forward_ref_reverse_UM_QF'] == {'A': 1}
    assert qualities == [30]


def test_determine_cutoff_tightens_to_max_expected_positions():
    qualities = [30] * 50
    k, expected_positions, error_perc = methods.determine_cutoff(
        qualities=qualities,
        reference_sequence='A',
        base_cutoff=1,
        error_cutoff=1.0,
        max_expected_positions=1e-6
    )
    assert isinstance(k, int)
    assert expected_positions <= 1e-6 * len('A') or k == 1
    assert 0.0 <= error_perc <= 100.0


def test_determine_cutoff_uses_normal_approximation_for_large_depth():
    qualities = [30] * 1001
    k, expected_positions, error_perc = methods.determine_cutoff(
        qualities=qualities,
        reference_sequence='A',
        base_cutoff=1,
        error_cutoff=1.0,
        max_expected_positions=0.01
    )

    assert isinstance(k, int)
    assert expected_positions >= 0.0
    assert 0.0 <= error_perc <= 100.0


def test_blast_result_parses_tab_delimited_line():
    line = 'query1 subject1 99.5 50 50 1 50 1 50 1e-10'
    result = cgdb.BlastResult(line)

    assert result.query_name == 'query1'
    assert result.subject_name == 'subject1'
    assert result.percent_identity == 99.5
    assert result.query_coverage == 100.0


def test_create_gene_allele_file_handles_slash_genus(tmp_path):
    profiles = tmp_path / 'profiles.txt'
    profiles.write_text(
        'genus\tBACT000001\tBACT000002\n'
        'Escherichia/Shigella\t1\t2\n',
        encoding='utf-8'
    )
    output = tmp_path / 'gene_allele.txt'

    genera = dbsetup.create_gene_allele_file(
        profiles_file=str(profiles),
        gene_allele_file=str(output)
    )

    assert 'Escherichia' in genera
    content = output.read_text(encoding='utf-8')
    assert content.startswith('Escherichia:BACT000001_1,BACT000002_2,')
