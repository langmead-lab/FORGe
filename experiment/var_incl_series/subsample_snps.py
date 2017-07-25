#!/usr/bin/env python

from __future__ import print_function
import sys
import random

random.seed(2677)

if len(sys.argv) < 2:
    raise RuntimeError('Specify a percentage as first arg')
pct = float(sys.argv[1])
pct /= 100.0
assert 0.0 <= pct <= 1.0
last_ln = None

# Sample consecutive pairs of lines
kept, total = 0, 0
for ln in sys.stdin:
    total += 1
    if last_ln is not None:
        if random.random() < pct:
            print(last_ln, end='')
            print(ln, end='')
            kept += 2
        last_ln = None
    else:
        last_ln = ln

print('Kept %d lines out of %d (%0.02f%%)' % (kept, total, kept*100.0/total), file=sys.stderr)
