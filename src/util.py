#!/usr/bin/env python

"""
Utility functions
"""

from future.utils import implements_iterator
import re
import string
import cProfile
from variant import VariantSet
from collections import OrderedDict


try:
    _revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
except AttributeError:
    _revcomp_trans = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")

_non_acgt = re.compile(b'[^ACGTacgt]')


def revcomp(x):
    return x[::-1].translate(_revcomp_trans)


def profileit(func):
    def wrapper(*args, **kwargs):
        datafn = func.__name__ + ".profile" # Name the data file sensibly
        prof = cProfile.Profile()
        retval = prof.runcall(func, *args, **kwargs)
        prof.dump_stats(datafn)
        return retval

    return wrapper


class Quantiler(object):

    def __init__(self):
        self.h = OrderedDict()

    def add(self, val):
        self.h[val] = self.h.get(val, 0) + 1

    def quantiles(self):
        size = sum(self.h.values())
        idx = {'min': 0,
               'q1': size / 4,
               'q2': size / 2,
               'q3': size * 3 / 4,
               'max': size-1}

        def accumulate(it):
            tot = 0
            for x in it:
                tot += x
                yield tot

        def get_val(d, i):
            return next(k for k, x in zip(d, accumulate(d.values())) if i < x)

        return {k: get_val(self.h, v) for k, v in idx.items()}

    def __str__(self):
        qu = self.quantiles()
        return '%d -- [%d, %d, %d] -- %d' % (qu['min'], qu['q1'], qu['q2'], qu['q3'], qu['max'])


def get_next_vector(k, counts, v):
    """
    Loop through all allele vectors in order
    counts: Array of length k storing the number of alternate alleles for each variant
    If v is the last array, return None
    """
    if v is None:
        v = [1] + ([0] * (k - 1))
        return v, False
    j = k-1
    while v[j] == counts[j]:
        # carry
        if j == 0:
            return None, True
        v[j] = 0
        j -= 1
    v[j] += 1
    return v, False


def vec_to_id(v, counts):
    """
    Convert an id vector, indicating which ALTs are present, into an integer
    that can be used as a concise key
    """
    ident = 0
    assert len(v) == len(counts)
    for i in range(len(v)):
        if v[i] > counts[i]:
            msg = 'Error in allele vector! Vector is %s but counts are %d' % (str(v[i]), counts[i])
            raise RuntimeError(msg)
        ident = ident * (counts[i]+1) + v[i]
    return ident


all_acgt_re = re.compile(b'^[ACGTacgt]*$')


@implements_iterator
class PseudocontigIterator(object):
    """ Loop over all pseudocontigs that contain a certain set of variants """

    def __init__(self, chrom_seq, variants, ids, r):
        assert isinstance(chrom_seq, bytes)
        self.seq = chrom_seq
        self.vars = variants
        self.ids = ids
        self.r = r
        self.k = len(ids)
        self.counts = [variants.num_alts(i) for i in ids]
        self.vec = None
        self.read_chunks = []
        for i in range(1, len(ids)):
            assert ids[i] > ids[i-1]
            assert variants.poss[ids[i]] > variants.poss[ids[i - 1]]
            # start from just after the previous REF
            st = variants.poss[ids[i - 1]] + len(variants.origs[ids[i - 1]])
            # take everything up to but not including this REF/ALT
            en = variants.poss[ids[i]]
            self.read_chunks.append(self.seq[st:en])

    def __iter__(self):
        return self

    def __next__(self):
        """
        Return current read and iterate to the next one
        """
        while True:
            self.vec, done = get_next_vector(self.k, self.counts, self.vec)
            if done:
                raise StopIteration
            variants, ids, vec = self.vars, self.ids, self.vec
            first_alt_base = 0
            read = variants.get_alt(ids[0], vec[0]-1)
            last_alt_base = len(read)
            for i in range(1, self.k):
                if vec[i] == 0:
                    read += self.read_chunks[i-1] + variants.origs[ids[i]]
                else:
                    if first_alt_base < 0:
                        first_alt_base = len(read)
                    read += self.read_chunks[i-1] + variants.get_alt(ids[i], vec[i]-1)
                    last_alt_base = len(read)

            # trim/extend ends of read to the right length
            end_offset = first_alt_base + self.r
            if end_offset <= last_alt_base:
                end_offset = last_alt_base+1
            rlen = len(read)
            if end_offset > rlen:
                pos = variants.poss[ids[-1]] + len(variants.origs[ids[-1]])
                suffix = self.seq[pos : pos+end_offset-rlen]
                read += suffix
            else:
                read = read[:end_offset]

            start_offset = last_alt_base - self.r
            if start_offset < 0:
                pos = variants.poss[ids[0]]
                assert pos >= 0
                prefix = self.seq[max(pos+start_offset, 0):pos]
                read = prefix + read
            else:
                read = read[start_offset:]

            if all_acgt_re.match(read):
                return read


def test_get_next_vector():
    vec = [0, 0]
    counts = [1, 1]
    vec, done = get_next_vector(2, counts, vec)
    assert not done
    assert [0, 1] == vec
    vec, done = get_next_vector(2, counts, vec)
    assert not done
    assert [1, 0] == vec
    vec, done = get_next_vector(2, counts, vec)
    assert not done
    assert [1, 1] == vec
    vec, done = get_next_vector(2, counts, vec)
    assert done


def test_pc_iter_1():
    seq = b'AAAAAAAAA'
    #      012345678
    #          T

    variants = VariantSet()
    variants.add_var(4, b'A', b'T', 0.25)

    pi = PseudocontigIterator(seq, variants, [0], 4)
    pcs = [x for x in pi]
    assert 1 == len(pcs)
    assert b'AAATAAA' == pcs[0]

    pi = PseudocontigIterator(seq, variants, [0], 5)
    pcs = [x for x in pi]
    assert 1 == len(pcs)
    assert b'AAAATAAAA' == pcs[0]

    pi = PseudocontigIterator(seq, variants, [0], 6)
    pcs = [x for x in pi]
    assert 1 == len(pcs)
    assert b'AAAATAAAA' == pcs[0]

    pi = PseudocontigIterator(seq, variants, [0], 9)
    pcs = [x for x in pi]
    assert 1 == len(pcs)
    assert b'AAAATAAAA' == pcs[0]


def test_pc_iter_2():
    seq = b'AAAAAAAAAA'
    #      0123456789
    #          CG

    variants = VariantSet()
    variants.add_var(4, b'A', b'C', 0.25)
    variants.add_var(5, b'A', b'G', 0.25)

    pi = PseudocontigIterator(seq, variants, [0, 1], 4)
    pcs = [x for x in pi]
    assert 2 == len(pcs)
    pcs.sort()
    assert b'AAACAAA' == pcs[0]
    assert b'AACGAA' == pcs[1]

    pi = PseudocontigIterator(seq, variants, [0, 1], 5)
    pcs = [x for x in pi]
    assert 2 == len(pcs)
    pcs.sort()
    assert b'AAAACAAAA' == pcs[0]
    assert b'AAACGAAA' == pcs[1]


def test_pc_iter_3():
    seq = b'AAAAAAAAAAA'
    #      01234567890
    #          CGT

    variants = VariantSet()
    variants.add_var(4, b'A', b'C', 0.25)
    variants.add_var(5, b'A', b'G', 0.25)
    variants.add_var(6, b'A', b'T', 0.25)

    pi = PseudocontigIterator(seq, variants, [0, 1, 2], 4)
    pcs = [x for x in pi]
    assert 4 == len(pcs)
    pcs.sort()
    assert b'AAACAAA' == pcs[0]
    assert b'AACGAA' == pcs[1]
    assert b'ACATA' == pcs[2]
    assert b'ACGTA' == pcs[3]


def test_pc_iter_4():
    seq = b'AAANAAAAA'
    #      012345678
    #          T

    variants = VariantSet()
    variants.add_var(4, b'A', b'T', 0.25)

    pi = PseudocontigIterator(seq, variants, [0], 4)
    pcs = [x for x in pi]
    assert 0 == len(pcs)


def test_quantiler_1():
    q = Quantiler()
    for i in range(100):
        q.add(i)
    qu = q.quantiles()
    assert 0 == qu['min']
    assert 25 == qu['q1']
    assert 50 == qu['q2']
    assert 75 == qu['q3']
    assert 99 == qu['max']
