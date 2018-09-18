#!/usr/bin/env python

from __future__ import print_function
import os
import gzip
import sys

rs_to_off = {}

with open('sm.1ksnp', 'rt') as fh:
    for ln in fh:
        toks = ln.rstrip().split()
        rs = toks[-1]
        off = toks[1]
        rs_to_off[rs] = off

with gzip.open('ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', 'rt') as fh:
    for ln in fh:
        if ln[0] == '#':
            print(ln, end='')
        else:
            toks = ln.rstrip().split()
            rs = toks[2]
            if rs in rs_to_off:
                toks[1] = rs_to_off[rs]
                print('\t'.join(toks))
                len_before = len(rs_to_off)
                del rs_to_off[rs]
                len_after = len(rs_to_off)
                print('len decreased from %d to %d' % (len_before, len_after), file=sys.stderr)
            if len(rs_to_off) == 0:
                print('saw all rss', file=sys.stderr)
                break
