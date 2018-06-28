#!/usr/bin/env python

"""
Interfaces to various k-mer counters
"""

from __future__ import print_function
import squeakr
import re
import string
import pytest
from bounter import bounter
from collections import Counter
import logging


_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
_non_acgt = re.compile('[^ACGTacgt]')


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
        self.counter.update(filter(lambda x: x is not None, slice_canonical(s, self.r)))

    def query(self, s):
        """ Query with each k-mer substring """
        return map(lambda kmer: -1 if kmer is None else self.counter.get(kmer), slice_canonical(s, self.r))


class BounterKmerCounter(object):
    def __init__(self,  name, r, size_mb=128, log_counting=8):
        if log_counting not in [None, 8, 1024]:
            raise ValueError('log_counting must be one of: None (for non-log), '
                             '8 (for 3-bit exact), 1024 (for 10-bit exact)')
        self.name = name
        self.r = r
        self.counter = bounter(size_mb=size_mb, need_iteration=False, need_counts=True, log_counting=log_counting)

    def add(self, s):
        """ Add canonicalized version of each k-mer substring """
        s = bytes(s)
        self.counter.update(filter(lambda x: x is not None, slice_canonical(s, self.r)))

    def query(self, s):
        """ Query with each k-mer substring """
        s = bytes(s)
        return map(lambda kmer: -1 if kmer is None else self.counter[kmer], slice_canonical(s, self.r))


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


class SqueakrKmerCounter(object):
    """
    What is the size of the CQF?
    """
    def __init__(self, name, r, qbits=10, seed=777, create_local=False):
        self.name = name
        self.r = r
        self.qbits = qbits
        self.seed = seed
        self.db = squeakr.cqf_new_db(r, qbits, seed, create_local=create_local)
        self.nadded = 0

    def add(self, s):
        self.nadded += squeakr.cqf_injest(self.db, s)

    def query(self, s):
        return squeakr.cqf_query(self.db, s)

    def report_fpr(self):
        results, results_len = squeakr.cqf_est_fpr(self.db)
        for i in range(0, results_len, 2):
            nobs, ntrue = results[i], results[i+1]
            logging.info('%d %d %d' % (nobs-ntrue, nobs, ntrue))


def pytest_generate_tests(metafunc):
    if 'counter_class' in metafunc.fixturenames:
        # Add classes to list as they are added
        counter_classes = [SimpleKmerCounter, BounterKmerCounter]
        metafunc.parametrize("counter_class", counter_classes, indirect=True)


@pytest.fixture(scope='session')
def counter_class(request):
    return request.param


def test_revcomp():
    assert 'ACGT' == revcomp('ACGT')
    assert 'A' == revcomp('T')
    assert '' == revcomp('')
    assert 'N' == revcomp('N')
    assert 'ANG' == revcomp('CNT')


def test_slice_1():
    lst = slice('ACGT', 1)
    assert list('ACGT') == lst


def test_slice_2():
    lst = slice('ACGT', 4)
    assert ['ACGT'] == lst
    lst = slice('ACNT', 4)
    assert [ None ] == lst


def test_slice_and_canonicalize():
    lst = slice_canonical('ACGT', 1)
    assert list('ACCA') == lst
    lst = slice_canonical('ACGT', 2)
    assert ['AC', 'CG', 'AC'] == lst


def test_counter_1(counter_class):
    cnt = counter_class('Test', 4)
    cnt.add('ACGT')
    res = cnt.query('ACGT')
    assert len(res) == 1
    assert res[0] == 1

    cnt.add('ACNT')
    res = cnt.query('ACGT')
    assert len(res) == 1
    assert res[0] == 1  # cumulative


def test_counter_2(counter_class):
    cnt = counter_class('Test', 8)
    cnt.add('ACGTACGTNACGTACGT')
    res = cnt.query('ACGTACGTN')
    assert len(res) == 2
    assert res[0] == 2
    assert res[1] == -1


def test_counter_3(counter_class):
    cnt = counter_class('Test', 4)
    cnt.add('ACGTACGTNACGTACGT')
    #        0000    x5555      ACGT x 2
    #         1111   x 6666     CGTA x 2
    #          2222  x  7777    GTAC x 2
    #           3333 x   8888   TACG -> CGTA x 2
    #            4444x    9999  ACGT x 2
    res = cnt.query('ACGTACG')
    #                0000
    #                 1111
    #                  2222
    #                   3333
    assert len(res) == 4
    assert res[0] == 4  # ACGT
    assert res[1] == 4  # CGTA
    assert res[2] == 2  # GTAC
    assert res[3] == 4  # TACG -> CGTA
