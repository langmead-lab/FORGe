#!/usr/bin/env python

"""
Interfaces to various k-mer counters
"""

from __future__ import print_function
from abc import ABCMeta, abstractmethod
import squeakr


class KmerCounter(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def add(self, s):
        pass

    @abstractmethod
    def query(self, s):
        pass


class JellyfishKmerCounter(KmerCounter):
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


class SqueakrKmerCounter(KmerCounter):
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
