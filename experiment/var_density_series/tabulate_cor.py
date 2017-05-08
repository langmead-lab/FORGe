#!/usr/bin/env python

from __future__ import print_function

cols = ['aligner', 'type', 'alts', 'stride', 'flag', 'cor']
print(','.join(cols))

for aligner in ['hisat2']:  # TODO: add ERG
    for type in ['indel', 'snp']:
        for alts in ['', '_alts']:
            for stride in [50, 100, 200, 400, 800]:
                cor_fn = '%s_base_stride_%d%s.cor' % (type, stride, alts)
                with open(cor_fn) as fh:
                    for ln in fh:
                        toks = ln.rstrip().split()
                        if (int(toks[0]) & 256) != 0:
                            continue
                        print(','.join([aligner, type, 'T' if alts == '_alts' else 'F', str(stride)] + toks))
