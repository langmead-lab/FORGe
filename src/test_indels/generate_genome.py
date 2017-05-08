#! /usr/bin/env python2.7

import sys
import random

l = 10000
G = ''
nts = ['A', 'C', 'G', 'T']
for _ in range(l):
    G += random.choice(nts)

with open(sys.argv[1], 'w') as f:
    f.write('>1\n')
    for i in range(0, l, 60):
        f.write(G[i:i+60] + '\n')
