#! /usr/bin/env python2.7

import sys
import random

indel_freq = 0.005

# Indel lengths will be chosen uniformly at random from this list
indel_lens = [1,2,3,4,5]
nts = ['A', 'C', 'G', 'T']
def sim_insert(l):
    ins = ''
    for _ in range(l):
        ins += random.choice(nts)
    return ins

G = ''
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == '>':
            continue
        G += line.rstrip()

l = len(G)

indels = []
for i in range(int(l * indel_freq)):
    added = False
    while not added:
        pos = random.randint(10, l-10)
        if not indels:
            indels = [pos]
            added = True
        for i in range(len(indels)):
            if pos > indels[i]-10 and pos < indels[i]:
                break
            elif pos < indels[i]-10:
                added = True
                indels = indels[:i] + [pos] + indels[i:]
                break
        if not added and pos > indels[-1] + 10:
            indels += [pos]
            added = True

with open('indels.txt', 'w') as f:
    offset = 0
    for pos in indels:
        if random.randint(1,2) == 1:
            # Deletion
            l = random.choice(indel_lens)
            f.write('1\t' + str(pos) + '\t-' + str(l) + '\n')
            G = G[:pos+offset] + G[pos+offset+l:]
            offset -= l

        else:
            # Insertion
            l = random.choice(indel_lens)
            f.write('1\t' + str(pos) + '\t' + str(l) + '\n')
            seq = sim_insert(l)
            G = G[:pos+offset] + seq + G[pos+offset:]
            offset += l

# Write genome to file
with open(sys.argv[2], 'w') as f:
    f.write('>1\n')
    for i in range(0, len(G), 60):
        f.write(G[i:i+60] + '\n')

