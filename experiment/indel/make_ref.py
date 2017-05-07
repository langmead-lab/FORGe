#!/usr/bin/env python

"""
python make_ref.py [snp or indel] [FASTA] [stride] [length]

Extract some sequence from an existing [FASTA], then create
an alternate version of that sequence by adding alternating
insertions and deletions every [stride] characters.  The
final genome, and its alternate version, will both be
[length] characters long.
"""

import os
import sys
import random

random.seed(463)

if len(sys.argv) < 2:
    raise RuntimeError('Must specify indel or snp as first argument')
type = sys.argv[1]

if len(sys.argv) < 3 or not os.path.exists(sys.argv[2]):
    raise RuntimeError('Must specify existing FASTA file as first argument')
fa = sys.argv[2]

if len(sys.argv) < 4:
    raise RuntimeError('Must specify stride length as second argument')
stride = int(sys.argv[3])

if len(sys.argv) < 5:
    raise RuntimeError('Must specify genome length as third argument')
length = int(sys.argv[4])

lines = []
lines_ln = 0
with open(fa) as fh:
    for ln in fh:
        if ln[0] == '>':
            continue
        if 'N' in ln:
            continue
        ln = ln.rstrip()
        lines.append(ln)
        lines_ln += len(ln)
        if lines_ln >= length:
            break

genome = ''.join(lines)[:length]

alt_genome = []
vars = []  # [id], [single, deletion, or insertion], [chromosome], [0-based off], [alt]

i = 0
genome_name = '%s_stride_%d' % (type, stride)
while i + 2*stride < length:
    alt_genome.append(genome[i:i+stride])
    # add insertion
    if type == 'indel':
        in_base = random.choice('ACGT')
        alt_genome.append(in_base)
        vars.append(['ins_%d' % (i+stride), 'insertion', genome_name, i+stride, in_base])
    else:
        in_base = random.choice('ACGT')
        while in_base == alt_genome[-1][-1]:
            in_base = random.choice('ACGT')
        alt_genome[-1] = alt_genome[-1][:-1] + in_base
        vars.append(['snp_%d' % (i+stride), 'single', genome_name, i+stride-1, in_base])

    alt_genome.append(genome[i+stride:i+stride+stride])
    # add deletion
    if type == 'indel':
        del_base = alt_genome[-1][-1]
        alt_genome[-1] = alt_genome[-1][:-1]
        vars.append(['del_%d' % (i+stride+stride-1), 'deletion', genome_name, i+stride+stride-1, '1'])
    else:
        in_base = random.choice('ACGT')
        while in_base == alt_genome[-1][-1]:
            in_base = random.choice('ACGT')
        alt_genome[-1] = alt_genome[-1][:-1] + in_base
        vars.append(['snp_%d' % (i+stride), 'single', genome_name, i+stride+stride-1, in_base])

    i += 2*stride

alt_genome.append(genome[i:i+(length-i)])

# Base reference sequence
with open('%s_base_stride_%d.fa' % (type, stride), 'w') as ofh:
    ofh.write('>%s_stride_%d\n' % (type, stride))
    ofh.write(genome)
    ofh.write('\n')

# Alt sequence
with open('%s_alt_stride_%d.fa' % (type, stride), 'w') as ofh:
    ofh.write('>%s_stride_%d\n' % (type, stride))
    ofh.write(''.join(alt_genome))
    ofh.write('\n')

# Variants with respect to base
with open('%s_base_stride_%d.snp' % (type, stride), 'w') as ofh:
    for var in vars:
        ofh.write('\t'.join(map(str, var)) + '\n')
