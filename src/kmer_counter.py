#!/usr/bin/env python

"""
Interfaces to various k-mer counters
"""

import jellyfish
import squeakr_query


class JellyfishKmerCounter(object):
    def __init__(self, r, initial_size=1024, bits_per_value=5):
        self.r = r
        jellyfish.MerDNA.k(r)
        self.counter = jellyfish.HashCounter(initial_size, bits_per_value)

    def add_string(self, s):
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
    def __init__(self, r, qbits, seed=777, create_local=False):
        self.r = r
        self.qbits = qbits
        self.seed = seed
        self.db = squeakr_query.cqf_new_db(r, qbits, seed, create_local=create_local)

    def add_string(self, s):
        squeakr_query.cqf_injest(self.db, s)

    def query(self, s):
        return squeakr_query.cqf_query(self.db, s)
