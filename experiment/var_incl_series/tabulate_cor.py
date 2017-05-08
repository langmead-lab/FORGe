#!/usr/bin/env python

from __future__ import print_function

cols = ['aligner', 'type', 'pct', 'stride', 'flag', 'cor']
pcts = [''] + list(map(lambda x: '_pct_%d_alts' % x, [10, 20, 30, 40, 60, 80, 100]))
strides = [450, 900]
print(','.join(cols))

for aligner in ['hisat2']:  # TODO: add ERG
    for type in ['indel', 'snp']:
        for pct in pcts:
            for stride in strides:
                cor_fn = '%s_base_stride_%d%s.cor' % (type, stride, pct)
                with open(cor_fn) as fh:
                    for ln in fh:
                        toks = ln.rstrip().split()
                        
                        if (int(toks[0]) & 256) != 0:
                            continue
                        pct_str = '0' if len(pct) == 0 else pct[5:-5]
                        print(','.join([aligner, type, pct_str, str(stride)] + toks))
