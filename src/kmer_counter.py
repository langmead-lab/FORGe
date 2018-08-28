#!/usr/bin/env python

"""
Interfaces to various k-mer counters
"""

from __future__ import print_function
import re
import os
import pytest
import logging
import shutil
import random
import tempfile
import subprocess
from util import revcomp, dump_to_dict
from collections import Counter


_non_acgt = re.compile(b'[^ACGTacgt]')


def slice_canonical(st, r):
    return map(lambda i: min(st[i:i+r], revcomp(st[i:i+r])) if not _non_acgt.search(st[i:i + r]) else None,
               range(len(st) - r + 1))


def slice(st, r):
    return map(lambda i: st[i:i+r] if not _non_acgt.search(st[i:i + r]) else None,
               range(len(st) - r + 1))


def run(cmd, stdout, stderr):
    logging.info(cmd)
    ret = os.system(cmd + ' >%s 2>%s' % (stdout, stderr))
    if ret != 0:
        msg = ('Error running command "%s"\n' % cmd) + \
              ('Exitlevel: %d\n' % ret) + \
               'Output:\n' + \
               open(stdout, 'r').read() + '\n' + \
               'Error:\n' + \
               open(stderr, 'r').read() + '\n'
        raise RuntimeError(msg)
    return ret


class SimpleKmerCounter(object):
    """
    Use collections.Counter and do exact counting
    """

    def __init__(self, name, r, dont_care_below=1,
                 temp=None, cache_from=None, cache_to=None):
        self.name = name
        self.r = r
        self.counter = Counter()
        self.dont_care_below = dont_care_below
        self.tmp_dir = tempfile.mkdtemp(dir=temp)

    def add(self, s):
        """ Add canonicalized version of each k-mer substring """
        assert isinstance(s, bytes)
        kmers = list(filter(lambda x: x is not None, slice_canonical(s, self.r)))
        self.counter.update(kmers)
        return len(kmers)

    def flush(self):
        pass  # nop

    def finalize(self):
        pass  # nop

    def cache_to(self):
        pass  # nop

    def query(self, query, threads=0, procid=0):
        """ Query with each k-mer substring """
        return list(map(lambda kmer: -1 if kmer is None else self.counter.get(kmer, 0), slice_canonical(query, self.r)))

    def query_batch(self, queries, threads=0, procid=0):
        """ Query with each k-mer substring of each query string """
        output_fn = os.path.join(self.tmp_dir, 'output.tsv')
        with open(output_fn, 'wb') as fh:
            for query in queries:
                for result in self.query(query):
                    fh.write(b'%d\n' % result)
        return output_fn


class KMC3KmerCounter(object):

    def __init__(self,  name, r, threads=1, gb=4, batch_size=-1,
                 dont_care_below=1, temp=None, cache_from=None, cache_to=None):
        self.name = name
        self.r = r
        self.tmp_dir = tempfile.mkdtemp(dir=temp)
        self.batch_mfa = os.path.join(self.tmp_dir, 'batch.mfa')
        self.working_dir = os.path.join(self.tmp_dir, 'working')
        assert not os.path.exists(self.working_dir)
        os.mkdir(self.working_dir)
        self.batch_db = os.path.join(self.tmp_dir, 'batch.db')
        self.combined_db = os.path.join(self.tmp_dir, name + '.db')
        self.combined_tmp_db = os.path.join(self.tmp_dir, 'combined_tmp.db')
        self.tmp_fh = open(self.batch_mfa, 'wb')
        self.query_mfa = os.path.join(self.tmp_dir, 'query.mfa')
        self.query_counts_db = os.path.join(self.tmp_dir, 'query_counts.db')
        self.buffered_bytes = 0
        self.batch_size = batch_size
        self.dont_care_below = dont_care_below
        self.first_flush = True
        self.stdout_log = os.path.join(self.tmp_dir, 'stdout.log')
        self.stderr_log = os.path.join(self.tmp_dir, 'stderr.log')
        self.output_fn = os.path.join(self.tmp_dir, 'query_output.log')
        self.gb = gb
        self.threads = threads
        self.cache_from_dir = cache_from
        self.cache_to_dir = cache_to
        self.cached_to = False
        self.batch_no = 0
        if cache_from is not None:
            for suf in ['.kmc_pre', '.kmc_suf']:
                fn = name + '.db' + suf
                if not os.path.exists(os.path.join(cache_from, fn)):
                    logging.warn('    cache-from enabled but could not find file "%s"' % fn)
                    self.cache_from_dir = None
                    break
            self.combined_db = os.path.join(cache_from, name + '.db')
        self.finalized = self.cache_from_dir is not None
        logging.info('    Initialized KMC3 k-mer counter with k=%d, gb=%d, threads=%d and batch_size=%d' %
                     (self.r, self.gb, self.threads, self.batch_size))

    def __del__(self):
        self.tmp_fh.close()
        to_delete = [self.batch_mfa, self.query_mfa, self.batch_db,
                     self.combined_tmp_db + '.kmc_pre',
                     self.combined_tmp_db + '.kmc_suf',
                     self.query_counts_db + '.kmc_pre',
                     self.query_counts_db + '.kmc_suf',
                     self.stdout_log,
                     self.stderr_log]
        if self.cached_to is None:
            to_delete.append(self.combined_db + '.kmc_pre')
            to_delete.append(self.combined_db + '.kmc_suf')
        for fn in to_delete:
            try:
                os.remove(fn)
            except OSError:
                pass
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir, ignore_errors=True)

    def add(self, s):
        if self.cache_from_dir is not None:
            raise RuntimeError('Cannot add to counter taken from cache')
        if self.finalized:
            raise RuntimeError('Cannot add after first query')
        assert isinstance(s, bytes)
        assert not self.tmp_fh.closed
        st = b'>r\n%s\n' % s
        self.buffered_bytes += len(st)
        assert self.buffered_bytes > 0
        self.tmp_fh.write(st)
        # negative batch size means no limit
        if self.buffered_bytes > self.batch_size >= 0:
            self.flush()
            assert self.buffered_bytes == 0

    def flush(self):
        threads = ('-t%d' % self.threads) if self.threads > 0 else ''
        if self.cache_from_dir is not None:
            raise RuntimeError('Cannot flush counter taken from cache')
        self.tmp_fh.close()
        logging.info('Counting batch of ~%d bytes' % self.buffered_bytes)
        run('kmc -sm -m%d %s -k%d -ci1 -fm %s %s %s' %
            (self.gb, threads, self.r, self.batch_mfa, self.batch_db, self.working_dir),
            self.stdout_log, self.stderr_log)
        if self.first_flush:
            self.first_flush = False
            os.system('mv %s.kmc_pre %s.kmc_pre' % (self.batch_db, self.combined_db))
            os.system('mv %s.kmc_suf %s.kmc_suf' % (self.batch_db, self.combined_db))
        else:
            assert self.batch_size >= 0
            run('kmc_tools %s, simple %s %s union %s' %
                (threads, self.batch_db, self.combined_db, self.combined_tmp_db),
                self.stdout_log, self.stderr_log)
            os.system('mv %s.kmc_pre %s.kmc_pre' % (self.combined_tmp_db, self.combined_db))
            os.system('mv %s.kmc_suf %s.kmc_suf' % (self.combined_tmp_db, self.combined_db))
            os.remove('%s.kmc_pre' % self.batch_db)
            os.remove('%s.kmc_suf' % self.batch_db)
        assert self.tmp_fh.closed
        os.remove(self.batch_mfa)
        self.tmp_fh = open(self.batch_mfa, 'wb')
        self.buffered_bytes = 0

    def finalize(self):
        threads = ('-t%d' % self.threads) if self.threads > 0 else ''
        if not self.finalized:
            if self.buffered_bytes > 0:
                self.flush()
            if self.cache_from_dir is not None:
                raise RuntimeError('Cannot finalize counter taken from cache')
            assert self.buffered_bytes == 0
            if self.dont_care_below > 1:
                os.system('mv %s.kmc_pre %s.kmc_pre' % (self.combined_db, self.combined_tmp_db))
                os.system('mv %s.kmc_suf %s.kmc_suf' % (self.combined_db, self.combined_tmp_db))
                run('kmc_tools %s transform %s reduce %s -ci%d' %
                    (threads, self.combined_tmp_db, self.combined_db, self.dont_care_below),
                    self.stdout_log, self.stderr_log)
                os.remove('%s.kmc_pre' % self.combined_tmp_db)
                os.remove('%s.kmc_suf' % self.combined_tmp_db)
            self.finalized = True
            self._cache_to()
        else:
            assert self.buffered_bytes == 0

    def _cache_to(self):
        if self.cache_to_dir is not None and not self.cached_to:
            name = os.path.basename(self.combined_db)
            full = os.path.join(self.cache_to_dir, name)
            for suf in ['.kmc_pre', '.kmc_suf']:
                shutil.move(self.combined_db + suf, full + suf)
                assert not os.path.exists(self.combined_db + suf)
                assert os.path.exists(full + suf)
            self.combined_db = full
            self.cached_to = True

    def query(self, query, threads=0, procid=0):
        assert isinstance(query, bytes)
        output_fn = self.query_batch([query], threads=threads, procid=procid)
        with open(output_fn, 'rb') as fh:
            return list(map(int, fh.readlines()))

    def query_batch(self, queries, threads=0, procid=0):
        """ Query with a batch of long strings """
        self.finalize()
        self.batch_no += 1
        if threads == 0:
            threads = self.threads
        threads = '-t%d' % threads
        # Write queries to multi-fasta
        query_mfa, batch_db = self.query_mfa, self.batch_db
        query_counts_db, working_dir = self.query_counts_db, self.working_dir
        output_fn = self.output_fn
        if procid > 0:
            query_mfa = os.path.join(self.tmp_dir,
                                     'query_thread%02d_batch%d.mfa' % (procid, self.batch_no))
            batch_db = os.path.join(self.tmp_dir,
                                    'batch_thread%02d_batch%d.db' % (procid, self.batch_no))
            query_counts_db = os.path.join(self.tmp_dir,
                                           'query_counts_thread%02d_batch%d.db' % (procid, self.batch_no))
            working_dir = os.path.join(self.tmp_dir,
                                       'working_thread%02d_batch%d' % (procid, self.batch_no))
            output_fn = os.path.join(self.tmp_dir,
                                     'query_output%02d_batch%d' % (procid, self.batch_no))
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        else:
            if not os.path.isdir(working_dir):
                raise RuntimeError('Working dir "%s" not a directory' % working_dir)
        with open(query_mfa, 'wb') as fh:
            for query in queries:
                st = b'>r\n%s\n' % query
                fh.write(st)
        # Count k-mers (counting is redundant; we just want a set)
        run('kmc -sm -m%d %s -k%d -ci1 -fm %s %s %s' %
            (self.gb, threads, self.r, query_mfa, batch_db, working_dir),
            self.stdout_log, self.stderr_log)
        os.remove(query_mfa)
        # Perform intersection with the big count file
        run('kmc_tools %s simple %s %s intersect %s -ocleft' %
            (threads, self.combined_db, batch_db, query_counts_db),
            self.stdout_log, self.stderr_log)
        # Can remove query set now
        os.remove('%s.kmc_pre' % batch_db)
        os.remove('%s.kmc_suf' % batch_db)
        # Dump the intersection file
        dump_cmd = 'kmc_dump %s /dev/stdout' % query_counts_db
        logging.info(dump_cmd)
        p = subprocess.Popen(dump_cmd, stdout=subprocess.PIPE, shell=True)
        h = dump_to_dict(p.stdout, chunksize=10000)
        p.wait()
        if p.returncode != 0:
            raise RuntimeError('Dump command failed: "%s"' % dump_cmd)
        # Can remove query counts now
        os.remove('%s.kmc_pre' % query_counts_db)
        os.remove('%s.kmc_suf' % query_counts_db)
        with open(output_fn, 'wb') as fh:
            for query in queries:
                for result in map(lambda kmer: -1 if kmer is None else h.get(kmer, 0),
                                  slice_canonical(query, self.r)):
                    fh.write(b'%d\n' % result)
        return output_fn


def pytest_generate_tests(metafunc):
    if 'counter_class' in metafunc.fixturenames:
        counter_classes = [SimpleKmerCounter, KMC3KmerCounter]
        metafunc.parametrize("counter_class", counter_classes, indirect=True)


@pytest.fixture(scope='session')
def counter_class(request):
    return request.param


def test_revcomp():
    assert b'ACGT' == revcomp(b'ACGT')
    assert b'A' == revcomp(b'T')
    assert b'' == revcomp(b'')
    assert b'N' == revcomp(b'N')
    assert b'ANG' == revcomp(b'CNT')


def test_slice_1():
    lst = list(slice(b'ACGT', 1))
    assert [b'A', b'C', b'G', b'T'] == lst


def test_slice_2():
    lst = list(slice(b'ACGT', 4))
    assert [b'ACGT'] == lst
    lst = list(slice(b'ACNT', 4))
    assert [ None ] == lst


def test_slice_and_canonicalize():
    lst = list(slice_canonical(b'ACGT', 1))
    assert [b'A', b'C', b'C', b'A'] == lst
    lst = list(slice_canonical(b'ACGT', 2))
    assert [b'AC', b'CG', b'AC'] == lst


def test_counter_1(counter_class):
    cnt = counter_class('Test', 4)
    cnt.add(b'ACGT')
    res = cnt.query(b'ACGT')
    assert 1 == len(res)
    assert 1 == res[0]

    cnt = counter_class('Test', 4)
    cnt.add(b'ACGT')
    cnt.add(b'ACNT')
    res = cnt.query(b'ACGT')
    assert 1 == len(res)
    assert 1 == res[0]


def test_counter_2(counter_class):
    cnt = counter_class('Test', 8)
    cnt.add(b'ACGTACGTNACGTACGT')
    #         ACGTACGT
    #                  ACGTACGT
    res = cnt.query(b'ACGTACGTN')
    assert 2 == len(res)
    assert -1 == res[1]
    assert 2 == res[0]


def test_counter_3(counter_class):
    cnt = counter_class('Test', 4)
    cnt.add(b'ACGTACGTNACGTACGT')
    #        0000    x5555      ACGT x 2
    #         1111   x 6666     CGTA x 2
    #          2222  x  7777    GTAC x 2
    #           3333 x   8888   TACG -> CGTA x 2
    #            4444x    9999  ACGT x 2
    res = cnt.query(b'ACGTACG')
    #                0000
    #                 1111
    #                  2222
    #                   3333
    assert 4 == len(res)
    assert res[0] == 4  # ACGT
    assert res[1] == 4  # CGTA
    assert res[2] == 2  # GTAC
    assert res[3] == 4  # TACG -> CGTA


def test_counter_4(counter_class):
    cnt = counter_class('Test', 4)
    cnt.add(b'ACGT')
    cnt.add(b'ACGT')
    cnt.add(b'ACGT')
    cnt.add(b'CCGT')
    cnt.add(b'CCGT')
    cnt.add(b'TTTT')
    cnt.add(b'AAAA')
    res_fn = cnt.query_batch([b'ACGT', b'CCGT', b'AAAAA'])
    #                        3        2        22
    with open(res_fn, 'rb') as fh:
        res = list(map(int, fh.readlines()))
    assert 4 == len(res)
    assert 3 == res[0]
    assert 2 == res[1]
    assert 2 == res[2]
    assert 2 == res[3]


def _test_random(cnt, k=4, stlen = 10000):
    random.seed(21361)
    st1 = b''.join([random.choice([b'A', b'C', b'G', b'T']) for _ in range(stlen)])
    st2 = b''.join([random.choice([b'A', b'C', b'G', b'T']) for _ in range(stlen)])
    results_len = stlen - k + 1
    cnt.add(st1)
    results = list(cnt.query(st2))
    assert results_len == len(results)
    for i in range(results_len-k+1):
        km_orig = st2[i:i+k]
        km = min(km_orig, revcomp(km_orig))
        true_count = st1.count(km)
        if true_count < 8:
            assert results[i] >= true_count, (results[i], km_orig, km, st1.count(km_orig), st1.count(revcomp(km_orig)))


def test_random_k4_10000(counter_class, k=4):
    _test_random(counter_class('Test', k), k)


def test_random_k100_100000(counter_class, k=100):
    _test_random(counter_class('Test', k), k)
