#!/usr/bin/enb python

"""
Utility functions
"""

import re


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
    for i in range(len(v)):
        if v[i] > counts[i]:
            msg = 'Error in allele vector! Vector is %s but counts are %d' % (str(v[i]), counts[i])
            raise RuntimeError(msg)
        ident = ident * (counts[i]+1) + v[i]
    return ident


all_acgt_re = re.compile("^[ACGTacgt]*$")


class PseudocontigIterator:
    """ Loop over all pseudocontigs that contain a certain set of variants """

    def __init__(self, chrom_seq, variants, ids, r):
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

    def next(self):
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


if __name__ == '__main__':
    import sys
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    from variant import VariantSet

    class TestGenNextVector(unittest.TestCase):

        def test_simple1(self):
            vec = [0, 0]
            counts = [1, 1]
            vec, done = get_next_vector(2, counts, vec)
            self.assertFalse(done)
            self.assertEqual([0, 1], vec)
            vec, done = get_next_vector(2, counts, vec)
            self.assertFalse(done)
            self.assertEqual([1, 0], vec)
            vec, done = get_next_vector(2, counts, vec)
            self.assertFalse(done)
            self.assertEqual([1, 1], vec)
            vec, done = get_next_vector(2, counts, vec)
            self.assertTrue(done)

    class TestPseudocontigIterator(unittest.TestCase):

        def test_simple1(self):
            seq = 'AAAAAAAAA'
            #      012345678
            #          T

            variants = VariantSet()
            variants.add_var(4, 'A', 'T', 0.25)

            pi = PseudocontigIterator(seq, variants, [0], 4)
            pcs = [x for x in pi]
            self.assertEqual(1, len(pcs))
            self.assertEqual('AAATAAA', pcs[0])

            pi = PseudocontigIterator(seq, variants, [0], 5)
            pcs = [x for x in pi]
            self.assertEqual(1, len(pcs))
            self.assertEqual('AAAATAAAA', pcs[0])

            pi = PseudocontigIterator(seq, variants, [0], 6)
            pcs = [x for x in pi]
            self.assertEqual(1, len(pcs))
            self.assertEqual('AAAATAAAA', pcs[0])

            pi = PseudocontigIterator(seq, variants, [0], 9)
            pcs = [x for x in pi]
            self.assertEqual(1, len(pcs))
            self.assertEqual('AAAATAAAA', pcs[0])

        def test_double1(self):
            seq = 'AAAAAAAAAA'
            #      0123456789
            #          CG

            variants = VariantSet()
            variants.add_var(4, 'A', 'C', 0.25)
            variants.add_var(5, 'A', 'G', 0.25)

            pi = PseudocontigIterator(seq, variants, [0, 1], 4)
            pcs = [x for x in pi]
            self.assertEqual(2, len(pcs))
            pcs.sort()
            self.assertEqual('AAACAAA', pcs[0])
            self.assertEqual('AACGAA', pcs[1])

            pi = PseudocontigIterator(seq, variants, [0, 1], 5)
            pcs = [x for x in pi]
            self.assertEqual(2, len(pcs))
            pcs.sort()
            self.assertEqual('AAAACAAAA', pcs[0])
            self.assertEqual('AAACGAAA', pcs[1])

        def test_triple1(self):
            seq = 'AAAAAAAAAAA'
            #      01234567890
            #          CGT

            variants = VariantSet()
            variants.add_var(4, 'A', 'C', 0.25)
            variants.add_var(5, 'A', 'G', 0.25)
            variants.add_var(6, 'A', 'T', 0.25)

            pi = PseudocontigIterator(seq, variants, [0, 1, 2], 4)
            pcs = [x for x in pi]
            self.assertEqual(4, len(pcs))
            pcs.sort()
            self.assertEqual('AAACAAA', pcs[0])
            self.assertEqual('AACGAA', pcs[1])
            self.assertEqual('ACATA', pcs[2])
            self.assertEqual('ACGTA', pcs[3])

        def test_invalid(self):
            seq = 'AAANAAAAA'
            #      012345678
            #          T

            variants = VariantSet()
            variants.add_var(4, 'A', 'T', 0.25)

            pi = PseudocontigIterator(seq, variants, [0], 4)
            pcs = [x for x in pi]
            self.assertEqual(0, len(pcs))

    unittest.main()
