#!/usr/bin/env python

"""
Interfaces to various k-mer counters
"""

import jellyfish


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
    def __init__(self, r):
        pass

    def add_string(self, s):
        pass

    def query(self, s):
        pass
