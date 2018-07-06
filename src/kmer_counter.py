#!/usr/bin/env python

"""
Interfaces to various k-mer counters
"""

from __future__ import print_function
import re
import string
import pytest
from collections import Counter
import logging
from _api import ffi, lib


# at module level
malloc = ffi.new_allocator(should_clear_after_alloc=False)

try:
    _revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
except AttributeError:
    _revcomp_trans = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")

_non_acgt = re.compile(b'[^ACGTacgt]')


def revcomp(x):
    return x[::-1].translate(_revcomp_trans)


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

    def __init__(self,  name, r):
        self.name = name
        self.r = r
        self.counter = Counter()

    def add(self, s):
        """ Add canonicalized version of each k-mer substring """
        assert isinstance(s, bytes)
        self.counter.update(filter(lambda x: x is not None, slice_canonical(s, self.r)))

    def query(self, s):
        """ Query with each k-mer substring """
        assert isinstance(s, bytes)
        result = list(map(lambda kmer: -1 if kmer is None else self.counter.get(kmer), slice_canonical(s, self.r)))
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
        counter_classes = [SimpleKmerCounter, SqueakrKmerCounter, BounterKmerCounter]
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


if __name__ == '__main__':
    if False:
        bk = BounterKmerCounter('blah', 8, size_mb=128)
        bk.add(b'ACGTACGTNACGTACGT')
        res = bk.query(b'ACGTACGTN')
        print(res)
        bk.close()
    else:
        bk = SqueakrKmerCounter('blah', 4, 10)
        bk.add(b'ACGTACGTNACGTACGT')
        res = bk.query(b'ACGTACGTN')
        print(res)
        bk.close()
