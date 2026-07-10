#!/usr/bin/env python3

"""
Tests for the ConFindr CLI entrypoint and database setup helpers.
"""

import os
import csv
import gzip
import json
import ssl
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

import confindr_src.confindr as conf
import confindr_src.create_genus_specific_db as cgdb
import confindr_src.database_setup as dbsetup
import confindr_src.methods as methods


def test_confindr_main_help_exits_success(monkeypatch):
    monkeypatch.setattr('sys.argv', ['confindr.py', '-h'])
    with pytest.raises(SystemExit) as excinfo:
        conf.main()
    assert excinfo.value.code == 0


def test_confindr_main_propagates_chunking_options(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()

    called = {}
    def fake_confindr(args):
        called['args'] = args
    monkeypatch.setattr(conf, 'confindr', fake_confindr)
    monkeypatch.setattr(sys, 'argv', [
        'confindr.py',
        '-i', str(input_dir),
        '-o', str(output_dir),
        '--contig-chunk-multiplier', '5',
        '--contig-chunk-bases', '1000'
    ])
    monkeypatch.setattr(conf, 'coloredlogs', None)

    conf.main()

    assert called['args'].forward_id == '_R1'
    assert methods.CONTIG_CHUNK_MULTIPLIER == 5
    assert methods.CONTIG_CHUNK_MAX_BASES == 1000


def test_confindr_main_invalid_chunking_errors(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    monkeypatch.setattr(sys, 'argv', [
        'confindr.py',
        '-i', str(input_dir),
        '-o', str(output_dir),
        '--contig-chunk-multiplier', 'foo',
        '--contig-chunk-bases', 'bar'
    ])
    monkeypatch.setattr(conf, 'coloredlogs', None)

    with pytest.raises(SystemExit) as excinfo:
        conf.main()
    assert excinfo.value.code == 2


def test_confindr_handles_pipeline_failure_and_writes_error_report(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    forward = input_dir / 'reads_R1.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='4g',
        cgmlst=None,
        data_type='Illumina',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)
    monkeypatch.setattr(conf, 'check_for_databases_and_download', lambda database_location: None)
    monkeypatch.setattr(conf, 'find_paired_reads', lambda fastq_directory, forward_id, reverse_id: [[str(forward)]])
    monkeypatch.setattr(conf, 'find_unpaired_reads', lambda *args, **kwargs: [])

    def fake_find_contamination(*args, **kwargs):
        raise subprocess.CalledProcessError(1, 'cmd')

    monkeypatch.setattr(conf, 'find_contamination', fake_find_contamination)
    monkeypatch.setattr(conf.shutil, 'rmtree', lambda path: None)

    conf.confindr(args=args)

    assert (output_dir / 'confindr_report.tsv').exists()


def test_confindr_confindr_calls_find_contamination(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    forward = input_dir / 'reads_R1.fastq.gz'
    reverse = input_dir / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='4g',
        cgmlst=None,
        data_type='Illumina',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    called = {'find_contamination': False}

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)
    monkeypatch.setattr(conf, 'check_for_databases_and_download', lambda database_location: None)
    monkeypatch.setattr(conf, 'find_paired_reads', lambda fastq_directory, forward_id, reverse_id: [[str(forward), str(reverse)]])
    monkeypatch.setattr(conf, 'find_unpaired_reads', lambda *args, **kwargs: [])

    def fake_find_contamination(*args, **kwargs):
        called['find_contamination'] = True
        return None

    monkeypatch.setattr(conf, 'find_contamination', fake_find_contamination)

    conf.confindr(args=args)

    assert called['find_contamination'] is True
    assert (output_dir).exists()


def test_confindr_cgmlst_nanopore_exits(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='4g',
        cgmlst=True,
        data_type='Nanopore',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)
    monkeypatch.setattr(conf, 'check_for_databases_and_download', lambda database_location: None)

    with pytest.raises(SystemExit) as excinfo:
        conf.confindr(args=args)
    assert excinfo.value.code == 1


def test_confindr_exits_when_dependency_missing(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='4g',
        cgmlst=None,
        data_type='Illumina',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: False)

    with pytest.raises(SystemExit) as excinfo:
        conf.confindr(args=args)
    assert excinfo.value.code == 1


def test_confindr_exits_on_invalid_xmx(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='5x',
        cgmlst=None,
        data_type='Illumina',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)

    with pytest.raises(SystemExit) as excinfo:
        conf.confindr(args=args)
    assert excinfo.value.code == 1


def test_confindr_debug_scans_position_files(monkeypatch, caplog, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()
    sample_dir = output_dir / 'reads'
    sample_dir.mkdir()
    pos_file = sample_dir / 'reads_positions.tsv'
    pos_file.write_text('contig\tposition\nchr1\t1\n', encoding='utf-8')

    forward = input_dir / 'reads_R1.fastq.gz'
    reverse = input_dir / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='4g',
        cgmlst=None,
        data_type='Illumina',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='debug',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)
    monkeypatch.setattr(conf, 'check_for_databases_and_download', lambda database_location: None)
    monkeypatch.setattr(conf, 'find_paired_reads', lambda fastq_directory, forward_id, reverse_id: [[str(forward), str(reverse)]])
    monkeypatch.setattr(conf, 'find_unpaired_reads', lambda *args, **kwargs: [])
    monkeypatch.setattr(conf, 'find_contamination', lambda *args, **kwargs: None)

    caplog.set_level('DEBUG')
    conf.confindr(args=args)

    assert 'first unique coords' in caplog.text


def test_database_setup_main_removes_existing_dir(monkeypatch, tmp_path):
    output_folder = tmp_path / 'db'
    output_folder.mkdir()
    (output_folder / 'old.txt').write_text('old', encoding='utf-8')
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    called = {'cgmlst': False, 'setup': False, 'mash': False}

    monkeypatch.setattr(dbsetup, 'download_cgmlst_derived_data', lambda output_folder: called.__setitem__('cgmlst', True))
    monkeypatch.setattr(dbsetup, 'setup_confindr_database', lambda **kwargs: called.__setitem__('setup', True))
    monkeypatch.setattr(dbsetup, 'download_mash_sketch', lambda output_folder: called.__setitem__('mash', True))
    monkeypatch.setattr(sys, 'argv', [
        'database_setup.py',
        '-o', str(output_folder),
        '-s', str(secret_file),
        '-u'
    ])

    original_context = dbsetup.ssl._create_default_https_context
    try:
        dbsetup.main()
    finally:
        dbsetup.ssl._create_default_https_context = original_context

    assert called['cgmlst'] is True
    assert called['setup'] is True
    assert called['mash'] is True
    assert not (output_folder / 'old.txt').exists()
    assert (output_folder / 'download_date.txt').exists()


def test_rmlstrest_get_loci_and_scheme_url_falls_back_to_text(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers, text):
            self.status_code = status_code
            self.headers = headers
            self.text = text
        def json(self):
            raise AssertionError('json should not be called')

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(200, {'content-type': 'text/plain'}, '{"loci": "http://example.com/loci", "schemes": "http://example.com/profile"}')

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'

    rmlst.get_loci_and_scheme_url()
    assert rmlst.loci == 'http://example.com/loci'
    assert rmlst.profile == 'http://example.com/profile'


def test_rmlstrest_get_session_token_success(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, json_data):
            self.status_code = status_code
            self._json_data = json_data
        def json(self):
            return self._json_data

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(200, {
                'oauth_token': 'sample_token',
                'oauth_token_secret': 'sample_secret'
            })

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.access_token = 'token'
    rmlst.access_secret = 'secret'

    rmlst.get_session_token()
    assert rmlst.session_token == 'sample_token'
    assert rmlst.session_secret == 'sample_secret'


def test_rmlstrest_download_loci_success_writes_files(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers=None, text=''):
            self.status_code = status_code
            self.headers = headers or {'content-type': 'text/plain'}
            self.text = text
        def json(self):
            try:
                return json.loads(self.text)
            except ValueError:
                return self.text

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            if url.endswith('/alleles_fasta'):
                return FakeResponse(200, {'content-type': 'text/plain'}, text='>locus_1\nACGT\n')
            return FakeResponse(200, {'content-type': 'application/json'}, text='{"loci": ["http://example.com/locus1"]}')

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.loci = 'http://example.com/locus1'
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'

    rmlst.download_loci()
    expected_file = tmp_path / 'locus1.tfa'
    assert expected_file.exists()
    assert expected_file.read_text(encoding='utf-8') == '>locus_1\nACGT\n'


def test_create_gene_allele_file_generates_genus_mapping(tmp_path):
    profiles = tmp_path / 'profiles.tsv'
    profiles.write_text(
        'genus\tBACT000001\tBACT000002\n'
        'Escherichia\t1\t2\n'
        'Escherichia/Shigella\t3\tN\n',
        encoding='utf-8'
    )
    output_file = tmp_path / 'genes.txt'

    genera = dbsetup.create_gene_allele_file(
        profiles_file=str(profiles),
        gene_allele_file=str(output_file)
    )

    assert genera == {'Escherichia'}
    contents = output_file.read_text(encoding='utf-8')
    assert 'Escherichia:BACT000001_1,BACT000002_2,BACT000001_3,' in contents


def test_confindr_nanopore_logs_warning(monkeypatch, caplog, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    forward = input_dir / 'reads_R1.fastq.gz'
    reverse = input_dir / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='4g',
        cgmlst=False,
        data_type='Nanopore',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)
    monkeypatch.setattr(conf, 'check_for_databases_and_download', lambda database_location: None)
    monkeypatch.setattr(conf, 'find_paired_reads', lambda fastq_directory, forward_id, reverse_id: [[str(forward), str(reverse)]])
    monkeypatch.setattr(conf, 'find_unpaired_reads', lambda *args, **kwargs: [])
    monkeypatch.setattr(conf, 'find_contamination', lambda *args, **kwargs: None)

    caplog.set_level('WARNING')
    conf.confindr(args=args)

    assert 'Nanopore contamination detection is highly experimental' in caplog.text


def test_confindr_removes_tmp_directory_on_completion(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()
    tmp_dir = tmp_path / 'tmp'
    tmp_dir.mkdir()
    (tmp_dir / 'temp.txt').write_text('tmp', encoding='utf-8')

    forward = input_dir / 'reads_R1.fastq.gz'
    reverse = input_dir / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_dir),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='4g',
        cgmlst=None,
        data_type='Illumina',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)
    monkeypatch.setattr(conf, 'check_for_databases_and_download', lambda database_location: None)
    monkeypatch.setattr(conf, 'find_paired_reads', lambda fastq_directory, forward_id, reverse_id: [[str(forward), str(reverse)]])
    monkeypatch.setattr(conf, 'find_unpaired_reads', lambda *args, **kwargs: [])
    monkeypatch.setattr(conf, 'find_contamination', lambda *args, **kwargs: None)

    conf.confindr(args=args)

    assert not tmp_dir.exists()


def test_confindr_exits_on_invalid_base_fraction(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=1.5,
        Xmx='4g',
        cgmlst=None,
        data_type='Illumina',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=None,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: False)

    with pytest.raises(SystemExit) as excinfo:
        conf.confindr(args=args)
    assert excinfo.value.code == 1


def test_confindr_main_installs_coloredlogs_when_available(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()

    monkeypatch.setattr(sys, 'argv', [
        'confindr.py',
        '-i', str(input_dir),
        '-o', str(output_dir),
        '-d', str(tmp_path / 'db')
    ])

    installed = {'called': False}
    class FakeColoredLogs:
        @staticmethod
        def install(**kwargs):
            installed['called'] = True
    monkeypatch.setattr(conf, 'coloredlogs', FakeColoredLogs)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)
    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_for_databases_and_download', lambda database_location: None)
    monkeypatch.setattr(conf, 'find_paired_reads', lambda fastq_directory, forward_id, reverse_id: [])
    monkeypatch.setattr(conf, 'find_unpaired_reads', lambda *args, **kwargs: [])
    monkeypatch.setattr(conf, 'find_contamination', lambda *args, **kwargs: None)

    conf.main()
    assert installed['called'] is True


def test_rmlstrest_get_loci_and_scheme_url_handles_json(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers, json_data):
            self.status_code = status_code
            self.headers = headers
            self._json_data = json_data
        def json(self):
            return self._json_data

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(200, {'content-type': 'application/json'}, {'loci': ['http://example.com/locus'], 'schemes': 'http://example.com/profile'})

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'

    rmlst.get_loci_and_scheme_url()
    assert rmlst.loci == ['http://example.com/locus']
    assert rmlst.profile == 'http://example.com/profile'


def test_rmlstrest_download_loci_failure_does_not_write(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers, text=None, json_data=None):
            self.status_code = status_code
            self.headers = headers
            self.text = text
            self._json_data = json_data
        def json(self):
            if self._json_data is not None:
                return self._json_data
            if self.text is not None:
                return self.text
            return {}

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            if url.endswith('/alleles_fasta'):
                return FakeResponse(404, {'content-type': 'text/plain'}, text='Not found')
            return FakeResponse(200, {'content-type': 'application/json'}, json_data={'loci': ['http://example.com/locus']})

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)
    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.loci = 'http://example.com/locus'
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'

    rmlst.download_loci()
    assert not (tmp_path / 'locus.tfa').exists()


def test_rmlstrest_download_profile_handles_json(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers, json_data):
            self.status_code = status_code
            self.headers = headers
            self._json_data = json_data
        def json(self):
            return self._json_data

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(200, {'content-type': 'application/json'}, 'genus\tBACT000001\n')

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)
    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.profile = 'http://example.com/profile'
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'

    rmlst.download_profile()
    assert (tmp_path / 'profiles.txt').exists()
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code):
            self.status_code = status_code
        def json(self):
            return {}

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(400)

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )

    with pytest.raises(SystemExit):
        rmlst.get_session_token()


def test_rmlstrest_get_request_and_access_token_workflow(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, json_data=None, text=''):
            self.status_code = status_code
            self._json_data = json_data
            self.text = text
        def json(self):
            return self._json_data

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True, params=None):
            if url.endswith('/oauth/get_session_token'):
                return FakeResponse(200, {'oauth_token': 'session', 'oauth_token_secret': 'secret'})
            if url.endswith('/oauth/get_access_token'):
                return FakeResponse(200, {'oauth_token': 'access', 'oauth_token_secret': 'access_secret'})
            return FakeResponse(200, {'content-type': 'application/json'}, text='{}')
        def request(self, method, url, params=None):
            return FakeResponse(200, {'oauth_token': 'request', 'oauth_token_secret': 'request_secret'})

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)
    monkeypatch.setattr('builtins.input', lambda prompt: 'verifier')

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.get_request_token()
    assert rmlst.request_token == 'request'
    assert rmlst.request_secret == 'request_secret'

    rmlst.get_access_token()
    assert rmlst.access_token == 'access'
    assert rmlst.access_secret == 'access_secret'


def test_rmlstrest_download_profile_writes_output(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')
    profile_path = tmp_path / 'profiles.txt'

    class FakeResponse:
        def __init__(self, status_code, text):
            self.status_code = status_code
            self.headers = {'content-type': 'text/plain'}
            self.text = text
        def json(self):
            return self.text

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(200, 'genus\tBACT000001\nEscherichia\t1\n')

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.profile = 'http://example.com/profile'
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'

    rmlst.download_profile()
    assert profile_path.exists()
    assert 'genus' in profile_path.read_text(encoding='utf-8')


def test_download_loci_failure_exits(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers, text=None, json_data=None):
            self.status_code = status_code
            self.headers = headers
            self.text = text
            self._json_data = json_data
        def json(self):
            if self._json_data is not None:
                return self._json_data
            return json.loads(self.text)

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            if url.endswith('/alleles_fasta'):
                return FakeResponse(404, {'content-type': 'text/plain'}, text='Not found')
            return FakeResponse(200, {'content-type': 'application/json'}, json_data={'loci': ['http://example.com/locus1']})

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)
    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.loci = 'http://example.com/locus1'
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'

    rmlst.download_loci()
    # The locus download failed, so the output file should not exist.
    assert not (tmp_path / 'locus1.tfa').exists()


def test_blast_result_parses_line():
    br = cgdb.BlastResult('query subject 99.0 100 100 1 100 1 100 1e-10')
    assert br.query_name == 'query'
    assert br.subject_name == 'subject'
    assert br.percent_identity == 99.0
    assert br.query_coverage == pytest.approx(100.0)


def test_get_potential_genes_selects_by_proportion(tmp_path):
    report = tmp_path / 'gene_report.tsv'
    report.write_text(
        'Gene\tOneHitPerGenome\n'
        'a\t1.0\n'
        'b\t0.5\n'
        'c\t0.4\n',
        encoding='utf-8'
    )
    result = cgdb.get_potential_genes(gene_report=str(report), desired_genes=2)
    assert result == ['a', 'b']


def test_check_for_similar_genes_creates_confirmed_genes(monkeypatch, tmp_path):
    gene1 = tmp_path / 'gene1.fasta'
    gene2 = tmp_path / 'gene2.fasta'
    genome1 = tmp_path / 'genome1.fasta'
    genome2 = tmp_path / 'genome2.fasta'

    gene1.write_text('>g1\nACGT\n', encoding='utf-8')
    gene2.write_text('>g2\nTGCA\n', encoding='utf-8')
    genome1.write_text('>g1\nACGT\n', encoding='utf-8')
    genome2.write_text('>g2\nTGCA\n', encoding='utf-8')

    def fake_call(cmd, shell=True):
        if 'blastn' in cmd:
            out_file = cmd.split('-out ')[1].split()[0]
            Path(out_file).write_text('', encoding='utf-8')
        return 0

    monkeypatch.setattr(subprocess, 'call', fake_call)
    result = cgdb.check_for_similar_genes(
        potential_genes=[str(gene1), str(gene2)],
        genomes=[str(genome1), str(genome2)]
    )
    assert set(result) == {str(gene1), str(gene2)}


def test_download_refseq_genomes_creates_files(monkeypatch, tmp_path):
    assembly = tmp_path / 'assembly_summary_refseq.txt'
    assembly.write_text(
        '#comment\n'
        'field0\tfield1\tfield2\tfield3\tfield4\tfield5\tfield6\tEscherichia coli\tfield8\tfield9\tfield10\tComplete\tfield12\tfield13\tfield14\tfield15\tfield16\tfield17\tfield18\tftp://example.com/genome\n',
        encoding='utf-8'
    )

    def fake_urlretrieve(url, dest):
        Path(dest).write_bytes(b'fakegzip')
        return url, dest

    def fake_call(cmd, shell=True):
        gz_path = cmd.split(' ')[1]
        out_path = gz_path.replace('.gz', '')
        Path(out_path).write_bytes(b'A' * 2100000)
        return 0

    monkeypatch.setattr('urllib.request.urlretrieve', fake_urlretrieve)
    monkeypatch.setattr(subprocess, 'call', fake_call)

    cgdb.download_refseq_genomes(
        output_folder=str(tmp_path),
        assembly_summary=str(assembly),
        genus='Escherichia'
    )
    assert (tmp_path / 'genome_1.fasta').exists()


def test_download_refseq_summary_creates_file(monkeypatch, tmp_path):
    output_folder = tmp_path / 'db'
    output_folder.mkdir()

    def fake_urlretrieve(url, dest):
        Path(dest).write_text('dummy', encoding='utf-8')
        return url, dest

    monkeypatch.setattr(cgdb.urllib.request, 'urlretrieve', fake_urlretrieve)

    cgdb.download_refseq_summary(output_folder=str(output_folder))

    assert (output_folder / 'assembly_summary_refseq.txt').exists()
    assert (output_folder / 'assembly_summary_refseq.txt').read_text() == 'dummy'


def test_find_hits_per_genome_writes_reports(monkeypatch, tmp_path):
    gene_folder = tmp_path / 'genes'
    genome_folder = tmp_path / 'genomes'
    gene_folder.mkdir()
    genome_folder.mkdir()

    gene = gene_folder / 'g1.fasta'
    genome = genome_folder / 'genome1.fasta'
    gene.write_text('>g1\nACGT\n', encoding='utf-8')
    genome.write_text('>s1\nACGT\n', encoding='utf-8')

    def fake_call(cmd, shell=True):
        if 'makeblastdb' in cmd:
            return 0
        if 'blastn' in cmd:
            out_file = cmd.split('-out ')[1].split()[0]
            Path(out_file).write_text(
                'query subject 95.0 100 100 1 100 1 100 1e-10\n',
                encoding='utf-8'
            )
            return 0
        return 0

    monkeypatch.setattr(cgdb.subprocess, 'call', fake_call)

    cgdb.find_hits_per_genome(
        genes_folder=str(gene_folder),
        genomes_folder=str(genome_folder)
    )

    assert (genome_folder / 'gene_hit_report.tsv').exists()
    assert (genome_folder / 'genome_hit_report.tsv').exists()


def test_confindr_passes_max_expected_positions(monkeypatch, tmp_path):
    input_dir = tmp_path / 'input'
    input_dir.mkdir()
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    db_dir = tmp_path / 'db'
    db_dir.mkdir()

    forward = input_dir / 'reads_R1.fastq.gz'
    reverse = input_dir / 'reads_R2.fastq.gz'
    with gzip.open(forward, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/1\nACGT\n+\nIIII\n')
    with gzip.open(reverse, 'wt', encoding='utf-8') as handle:
        handle.write('@r1/2\nACGT\n+\nIIII\n')

    args = SimpleNamespace(
        input_directory=str(input_dir),
        output_name=str(output_dir),
        databases=str(db_dir),
        threads=1,
        tmp=str(tmp_path / 'tmp'),
        keep_files=False,
        quality_cutoff=20,
        min_quality=15,
        base_cutoff=3,
        base_fraction_cutoff=0.05,
        Xmx='4g',
        cgmlst=None,
        data_type='Illumina',
        rmlst=False,
        use_prob_scoring=False,
        score_threshold=None,
        min_matching_hashes=150,
        max_expected_positions=0.01,
        fasta=False,
        forward_id='_R1',
        reverse_id='_R2',
        verbosity='info',
        downsample_depth=None,
        subreplicates=1,
        subreplicate_seed=None,
        subreplicate_consensus=0.5,
    )

    captured = {'value': None}

    monkeypatch.setattr(conf, 'dependency_check', lambda dependency: True)
    monkeypatch.setattr(conf, 'check_valid_base_fraction', lambda base_fraction: True)
    monkeypatch.setattr(conf, 'check_for_databases_and_download', lambda database_location: None)
    monkeypatch.setattr(conf, 'find_paired_reads', lambda fastq_directory, forward_id, reverse_id: [[str(forward), str(reverse)]])
    monkeypatch.setattr(conf, 'find_unpaired_reads', lambda *args, **kwargs: [])

    def fake_find_contamination(*args, **kwargs):
        captured['value'] = kwargs.get('max_expected_positions')
        return None

    monkeypatch.setattr(conf, 'find_contamination', fake_find_contamination)

    conf.confindr(args=args)

    assert captured['value'] == 0.01


def test_create_genus_specific_db_main_calls_steps(monkeypatch, tmp_path):
    output_folder = tmp_path / 'out'
    input_folder = tmp_path / 'input'
    output_folder.mkdir()
    input_folder.mkdir()
    gene = input_folder / 'g1.fasta'
    gene.write_text('>g1\nACGT\n', encoding='utf-8')

    called = []

    monkeypatch.setattr(sys, 'argv', [
        'create_genus_specific_db.py',
        '-o', str(output_folder),
        '-i', str(input_folder),
        '-g', 'Escherichia'
    ])
    monkeypatch.setattr(cgdb, 'download_refseq_summary', lambda output_folder: called.append('summary'))
    monkeypatch.setattr(cgdb, 'download_refseq_genomes', lambda output_folder, assembly_summary, genus: called.append('genomes'))
    monkeypatch.setattr(cgdb, 'find_hits_per_genome', lambda genes_folder, genomes_folder: called.append('hits'))
    monkeypatch.setattr(cgdb, 'get_potential_genes', lambda gene_report, desired_genes: [str(gene)])
    monkeypatch.setattr(cgdb, 'check_for_similar_genes', lambda potential_genes, genomes: potential_genes)

    def fake_call(cmd, shell=True):
        if 'cat ' in cmd and ' >> ' in cmd:
            source, destination = cmd.split(' >> ')
            source = source.split('cat ')[1]
            Path(destination).write_text(Path(source).read_text(), encoding='utf-8')
        return 0

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(cgdb.subprocess, 'call', fake_call)

    cgdb.main()

    assert 'summary' in called
    assert 'genomes' in called
    assert 'hits' in called
    output_file = tmp_path / 'Escherichia_db_cgderived.fasta'
    assert output_file.exists()
    output_file.unlink()


def test_database_setup_main_with_secret_file_and_unverified(monkeypatch, tmp_path):
    output_folder = tmp_path / 'db'
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    called = {'cgmlst': False, 'setup': False, 'mash': False}

    monkeypatch.setattr(dbsetup, 'download_cgmlst_derived_data', lambda output_folder: called.__setitem__('cgmlst', True))
    monkeypatch.setattr(dbsetup, 'setup_confindr_database', lambda **kwargs: called.__setitem__('setup', True))
    monkeypatch.setattr(dbsetup, 'download_mash_sketch', lambda output_folder: called.__setitem__('mash', True))
    monkeypatch.setattr(sys, 'argv', [
        'database_setup.py',
        '-o', str(output_folder),
        '-s', str(secret_file),
        '-u'
    ])

    dbsetup.main()

    assert called['cgmlst'] is True
    assert called['setup'] is True
    assert called['mash'] is True
    assert (output_folder / 'download_date.txt').exists()


def test_database_setup_main_without_secret_file_skips_setup(monkeypatch, tmp_path):
    output_folder = tmp_path / 'db'
    called = {'cgmlst': False, 'setup': False, 'mash': False}

    monkeypatch.setattr(dbsetup, 'download_cgmlst_derived_data', lambda output_folder: called.__setitem__('cgmlst', True))
    monkeypatch.setattr(dbsetup, 'setup_confindr_database', lambda **kwargs: called.__setitem__('setup', True))
    monkeypatch.setattr(dbsetup, 'download_mash_sketch', lambda output_folder: called.__setitem__('mash', True))
    monkeypatch.setattr(sys, 'argv', [
        'database_setup.py',
        '-o', str(output_folder)
    ])

    dbsetup.main()

    assert called['cgmlst'] is True
    assert called['mash'] is True
    assert called['setup'] is False
    assert (output_folder / 'download_date.txt').exists()


def test_setup_confindr_database_creates_combined_files(monkeypatch, tmp_path):
    output_folder = tmp_path / 'db'
    output_folder.mkdir()
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeRmlstRest:
        def __init__(self, consumer_secret_file, output_folder, unverified=False):
            assert consumer_secret_file == str(secret_file)
            assert output_folder == str(output_folder)
            self.output_folder = output_folder
            self.unverified = unverified
        def get_request_token(self):
            pass
        def get_access_token(self):
            pass
        def get_session_token(self):
            pass
        def get_loci_and_scheme_url(self):
            pass
        def download_loci(self):
            for idx in [1, 2]:
                path = Path(self.output_folder) / f'BACT00000{idx}.tfa'
                path.write_text(f'>BACT00000{idx}_1\nACGT\n', encoding='utf-8')
        def download_profile(self):
            headers = ['genus'] + [f'BACT00000{i}' for i in range(1, 66)]
            row = ['Escherichia'] + [str(i) for i in range(1, 66)]
            profile = Path(self.output_folder) / 'profiles.txt'
            profile.write_text('\t'.join(headers) + '\n' + '\t'.join(row) + '\n', encoding='utf-8')

    called = {'index': False}
    monkeypatch.setattr(dbsetup, 'RmlstRest', FakeRmlstRest)
    monkeypatch.setattr(dbsetup, 'index', lambda **kwargs: called.__setitem__('index', True))

    dbsetup.setup_confindr_database(
        output_folder=str(output_folder),
        consumer_secret=str(secret_file),
        index_databases=True,
        unverified=True
    )

    assert (output_folder / 'rMLST_combined.fasta').exists()
    assert (output_folder / 'gene_allele.txt').exists()
    assert called['index']


def test_rmlstrest_download_profile_writes_profiles(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers, text):
            self.status_code = status_code
            self.headers = headers
            self.text = text
    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return FakeResponse(200, {'content-type': 'text/plain'}, 'genus\tBACT000001\n')
    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'
    rmlst.profile = 'http://example.com/profile'

    rmlst.download_profile()
    assert (tmp_path / 'profiles.txt').read_text() == 'genus\tBACT000001\n'


def test_rmlstrest_get_request_token_and_access_token(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')
    responses = [
        SimpleNamespace(status_code=200, json=lambda: {'oauth_token': 'token', 'oauth_token_secret': 'secret'}),
        SimpleNamespace(status_code=200, json=lambda: {'oauth_token': 'token2', 'oauth_token_secret': 'secret2'})
    ]

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def request(self, method, url, params=None):
            return responses.pop(0)
        def get(self, url, verify=True, params=None):
            return responses.pop(0)

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)
    monkeypatch.setattr('builtins.input', lambda prompt='': 'verifier')

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.get_request_token()
    assert rmlst.request_token == 'token'
    assert rmlst.request_secret == 'secret'
    rmlst.get_access_token()
    assert rmlst.access_token == 'token2'
    assert rmlst.access_secret == 'secret2'


def test_rmlstrest_get_session_token_sets_tokens(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, json_data):
            self.status_code = status_code
            self._json_data = json_data
        def json(self):
            return self._json_data

    responses = [
        FakeResponse(200, {'oauth_token': 'token', 'oauth_token_secret': 'secret'})
    ]

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return responses.pop(0)

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )

    rmlst.get_session_token()

    assert rmlst.session_token == 'token'
    assert rmlst.session_secret == 'secret'


def test_rmlstrest_get_loci_and_scheme_url_decodes_json(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers, json_data):
            self.status_code = status_code
            self.headers = headers
            self._json_data = json_data
        def json(self):
            return self._json_data

    responses = [
        FakeResponse(200, {'content-type': 'application/json'}, {'loci': 'http://example.com/loci', 'schemes': 'http://example.com/profile'})
    ]

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True):
            return responses.pop(0)

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'

    rmlst.get_loci_and_scheme_url()

    assert rmlst.loci == 'http://example.com/loci'
    assert rmlst.profile == 'http://example.com/profile'


def test_rmlstrest_download_loci_and_profile_writes_files(monkeypatch, tmp_path):
    secret_file = tmp_path / 'secret.txt'
    secret_file.write_text('key\nsecret\n', encoding='utf-8')

    class FakeResponse:
        def __init__(self, status_code, headers, text=None, json_data=None):
            self.status_code = status_code
            self.headers = headers
            self.text = text
            self._json_data = json_data
        def json(self):
            return self._json_data

    responses = [
        FakeResponse(200, {'content-type': 'application/json'}, json_data={'loci': ['http://example.com/locus1']}),
        FakeResponse(200, {'content-type': 'text/plain'}, text='>allele1\nACGT\n'),
        FakeResponse(200, {'content-type': 'text/plain'}, text='genus\tBACT000001\n')
    ]

    class FakeSession:
        def __init__(self, *args, **kwargs):
            pass
        def get(self, url, verify=True, params=None):
            return responses.pop(0)

    monkeypatch.setattr(dbsetup, 'OAuth1Session', FakeSession)

    rmlst = dbsetup.RmlstRest(
        consumer_secret_file=str(secret_file),
        output_folder=str(tmp_path),
        unverified=False
    )
    rmlst.session_token = 'token'
    rmlst.session_secret = 'secret'
    rmlst.loci = 'http://example.com/locus1'
    rmlst.profile = 'http://example.com/profile1'

    rmlst.download_loci()
    assert (tmp_path / 'locus1.tfa').exists()
    assert '>allele1' in (tmp_path / 'locus1.tfa').read_text()

    rmlst.download_profile()
    assert (tmp_path / 'profiles.txt').exists()
    assert 'genus' in (tmp_path / 'profiles.txt').read_text()


def test_create_gene_allele_file_writes_mapping(tmp_path):
    profiles = tmp_path / 'profiles.txt'
    profiles.write_text(
        'genus\tBACT000001\tBACT000002\n'
        'Escherichia\t1\tN\n'
        'Listeria/Shigella\t2\t3\n',
        encoding='utf-8'
    )
    out_file = tmp_path / 'gene_allele.txt'
    genera = dbsetup.create_gene_allele_file(
        profiles_file=str(profiles),
        gene_allele_file=str(out_file)
    )
    assert 'Escherichia' in genera
    assert (out_file).exists()
    contents = out_file.read_text(encoding='utf-8')
    assert 'Escherichia:BACT000001_1,' in contents


def test_setup_confindr_database_combines_loci_and_profiles(monkeypatch, tmp_path):
    output_folder = tmp_path / 'db'
    output_folder.mkdir()
    # Create fake locus files
    for locus in ['BACT000001.tfa', 'BACT000002.tfa']:
        Path(output_folder / locus).write_text('>allele1\nACGT\n', encoding='utf-8')
    profiles = output_folder / 'profiles.txt'
    profiles.write_text('genus\tBACT000001\tBACT000002\nEscherichia\t1\t1\n', encoding='utf-8')

    class FakeRmlstRest:
        def __init__(self, consumer_secret_file, output_folder, unverified=False):
            pass
        def get_request_token(self):
            pass
        def get_access_token(self):
            pass
        def get_session_token(self):
            pass
        def get_loci_and_scheme_url(self):
            pass
        def download_loci(self):
            pass
        def download_profile(self):
            pass

    monkeypatch.setattr(dbsetup, 'RmlstRest', FakeRmlstRest)
    monkeypatch.setattr(dbsetup, 'index', lambda *args, **kwargs: None)

    dbsetup.setup_confindr_database(
        output_folder=str(output_folder),
        consumer_secret=str(tmp_path / 'secret.txt'),
        index_databases=False,
        unverified=False
    )

    assert (output_folder / 'rMLST_combined.fasta').exists()
    assert (output_folder / 'gene_allele.txt').exists()
