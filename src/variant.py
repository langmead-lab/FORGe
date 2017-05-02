#! /usr/bin/env python2.7

class Variant:
    def __init__(self, name, chrom, pos, orig, alts, probs):
        self.name = name
        self.chrom = chrom
        self.pos = pos
        self.orig = orig
        self.alts = alts
        self.probs = probs
        self.num_alts = len(alts)

    def add_alt(self, alt, prob):
        self.alts.append(alt)
        self.probs.append(prob)
        self.num_alts += 1
