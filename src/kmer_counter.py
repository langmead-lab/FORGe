#!/usr/bin/env python

"""
Interfaces to various k-mer counters
"""

from __future__ import print_function
import re
import pytest
from collections import Counter
import logging
import random
from util import revcomp
from _api import ffi, lib


# at module level
malloc = ffi.new_allocator(should_clear_after_alloc=False)

_non_acgt = re.compile(b'[^ACGTacgt]')


def slice_canonical(st, r):
    return map(lambda i: min(st[i:i+r], revcomp(st[i:i+r])) if not _non_acgt.search(st[i:i + r]) else None,
               range(len(st) - r + 1))


def slice(st, r):
    return map(lambda i: st[i:i+r] if not _non_acgt.search(st[i:i + r]) else None,
               range(len(st) - r + 1))


class SimpleKmerCounter(object):
    """
    Use collections.Counter and do exact counting
    """

    def __init__(self, name, r):
        self.name = name
        self.r = r
        self.counter = Counter()

    def add(self, s):
        """ Add canonicalized version of each k-mer substring """
        assert isinstance(s, bytes)
        kmers = list(filter(lambda x: x is not None, slice_canonical(s, self.r)))
        self.counter.update(kmers)
        return len(kmers)

    def query(self, s):
        """ Query with each k-mer substring """
        assert isinstance(s, bytes)
        result = list(map(lambda kmer: -1 if kmer is None else self.counter.get(kmer, 0), slice_canonical(s, self.r)))
        assert all(map(lambda x: x is not None, result)), (s, result)
        return result, len(result)


class BounterKmerCounter(object):
    def __init__(self,  name, r, size_mb=128, log_counting=8):
        if log_counting not in [None, 8, 1024]:
            raise ValueError('log_counting must be one of: None (for non-log), '
                             '8 (for 3-bit exact), 1024 (for 10-bit exact)')
        if log_counting != 8:
            raise ValueError('Only 8-bit log counting supported for now')
        if log_counting == 8:
            self.cell_size = cell_size = 1
        elif log_counting == 1024:
            self.cell_size = cell_size = 2
        else:
            self.cell_size = cell_size = 4
        self.width = 1 << (size_mb * (2 ** 20) // (cell_size * 8 * 2)).bit_length()
        self.depth = (size_mb * (2 ** 20)) // (self.width * cell_size)
        self.name = name
        self.r = r
        self.result_len = 1024
        self.results = malloc("int64_t[]", self.result_len)
        self.sketch = lib.bounter_new(self.width, self.depth)

    def close(self):
        lib.bounter_delete(self.sketch)

    def add(self, s):
        """ Add canonicalized version of each k-mer substring """
        assert isinstance(s, bytes)
        lib.bounter_string_injest(self.sketch, self.r, s, len(s))

    def query(self, s):
        """ Query with each k-mer substring """
        assert isinstance(s, bytes)
        result_len = len(s)-self.r+1
        if result_len > self.result_len:
            del self.results
            self.result_len = result_len
            self.results = malloc("int64_t[]", self.result_len)
        lib.bounter_string_query(self.sketch, self.r, s, len(s), self.results, result_len)
        return self.results, result_len


class SqueakrKmerCounter(object):
    """
    What is the size of the CQF?
    """
    def __init__(self, name, r, qbits=10):
        self.name = name
        self.r = r
        self.qbits = qbits
        self.result_len = 1024
        self.results = malloc("int64_t[]", self.result_len)
        self.db = lib.cqf_new(r, qbits)

    def close(self):
        lib.cqf_delete(self.db)

    def add(self, s):
        """ Add canonicalized version of each k-mer substring """
        assert isinstance(s, bytes)
        lib.cqf_string_injest(s, len(s), self.db)

    def query(self, s):
        """ Query with each k-mer substring """
        assert isinstance(s, bytes)
        result_len = len(s)-self.r+1
        if result_len > self.result_len:
            del self.results
            self.result_len = result_len
            self.results = malloc("int64_t[]", self.result_len)
        lib.cqf_string_query(s, len(s), self.results, result_len, self.db)
        return self.results, result_len

    def report_fpr(self):
        results, results_len = lib.cqf_est_fpr(self.db)
        for i in range(0, results_len, 2):
            nobs, ntrue = results[i], results[i+1]
            logging.info('%d %d %d' % (nobs-ntrue, nobs, ntrue))


class JellyfishKmerCounter(object):
    def __init__(self,  name, r, initial_size=1024, bits_per_value=5):
        import jellyfish
        self.name = name
        self.r = r
        jellyfish.MerDNA.k(r)
        self.counter = jellyfish.HashCounter(initial_size, bits_per_value)

    def add(self, s):
        """ Add canonicalized version of each k-mer substring """
        for mer in jellyfish.string_canonicals(s):
            self.counter.add(mer, 1)

    def query(self, s):
        res = []
        for mer in jellyfish.string_canonicals(s):
            cnt = self.counter[mer]
            assert cnt is not None
            res.append(cnt)
        return res


def pytest_generate_tests(metafunc):
    if 'counter_class' in metafunc.fixturenames:
        # Add classes to list as they are added
        #counter_classes = [SimpleKmerCounter, BounterKmerCounter, SqueakrKmerCounter]
        #counter_classes = [SimpleKmerCounter, BounterKmerCounter]
        counter_classes = [SimpleKmerCounter, BounterKmerCounter]
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
    res, reslen = cnt.query(b'ACGT')
    assert reslen == 1
    assert res[0] == 1

    cnt.add(b'ACNT')
    res, reslen = cnt.query(b'ACGT')
    assert reslen == 1
    assert res[0] == 1  # cumulative


def test_counter_2(counter_class):
    cnt = counter_class('Test', 8)
    cnt.add(b'ACGTACGTNACGTACGT')
    #         ACGTACGT
    #                  ACGTACGT
    res, reslen = cnt.query(b'ACGTACGTN')
    assert reslen == 2
    assert res[1] == -1
    assert res[0] == 2


def test_counter_3(counter_class):
    cnt = counter_class('Test', 4)
    cnt.add(b'ACGTACGTNACGTACGT')
    #        0000    x5555      ACGT x 2
    #         1111   x 6666     CGTA x 2
    #          2222  x  7777    GTAC x 2
    #           3333 x   8888   TACG -> CGTA x 2
    #            4444x    9999  ACGT x 2
    res, reslen = cnt.query(b'ACGTACG')
    #                0000
    #                 1111
    #                  2222
    #                   3333
    assert reslen == 4
    assert res[0] == 4  # ACGT
    assert res[1] == 4  # CGTA
    assert res[2] == 2  # GTAC
    assert res[3] == 4  # TACG -> CGTA


def _test_random(cnt, k=4, len=1000):
    stlen = 10000
    random.seed(21361)
    st1 = b''.join([random.choice([b'A', b'C', b'G', b'T']) for _ in range(stlen)])
    st2 = b''.join([random.choice([b'A', b'C', b'G', b'T']) for _ in range(stlen)])
    results_len = stlen - k + 1
    cnt.add(st1)
    results, results_len2 = cnt.query(st2)
    assert results_len == results_len2
    for i in range(results_len-k+1):
        km_orig = st2[i:i+k]
        km = min(km_orig, revcomp(km_orig))
        true_count = st1.count(km)
        if true_count < 8:
            assert results[i] >= true_count, (results[i], km_orig, km, st1.count(km_orig), st1.count(revcomp(km_orig)))


def test_random_k4_10000(counter_class, k=4):
    _test_random(counter_class('Test', k), k, 10000)


def test_random_k100_100000(counter_class, k=100):
    _test_random(counter_class('Test', k), k, 100000)
