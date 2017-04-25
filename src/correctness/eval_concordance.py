"""
Turn ZT:Zs and MAPQs into percentiles for further analysis.
"""

from __future__ import print_function
from collections import defaultdict
import sys
import re
import logging


"""
First pass:
- Collect MAPQ distributions for tandem & input
- Collect ZT:Z distributions for tandem & input

Second pass:
- For tandem and input:
  + For each read, check for incorrectness.  If incorrect:
    - Record vector of percentiles for MAPQ and ZT:Zs

Finally:
- For tandem and input:
  + Sort the percentile vector by MAPQ percentile and record

Example usage:
pypy eval_concordance.py $HOME/input.sam $HOME/predictions.csv dists.txt incor.mapq incor.mapq_orig incor.both

TODO:
- Are we doing paired-end correctly for the various name formats?

# Load data and give appropriate column names
setwd('~/git/qsim-experiments/bin')
mq <- read.table('incor.mapq', sep=',', header=F, quote="", comment.char="")
mo <- read.table('incor.mapq_orig', sep=',', header=F, quote="", comment.char="")
mb <- read.table('incor.both', sep=',', header=F, quote="", comment.char="")
ztn <- c('ztz0', 'ztz1', 'ztz2', 'ztz3', 'ztz4', 'ztz5', 'ztz6', 'ztz7', 'ztz8', 'ztz9')
cn <- c('mapq', ztn)
colnames(mq) <- cn
colnames(mo) <- cn
colnames(mb) <- c('mapq', 'mapq_orig', ztn)

# Plot some basics
plot(jitter(mq$ztz1), jitter(mq$mapq))
points(jitter(mo$ztz1), jitter(mo$mapq), col='blue')

# Plot mapq difference as a function of diff percentile
library(ggplot2)
library(ggExtra)
p <- ggplot(mb, aes(x=ztz1, y=mapq - mapq_orig)) + geom_point() + theme_classic()
ggExtra::ggMarginal(p)
"""


def pass1_fh(fh_sam, fh_pred):
    """ Histograms all the alignments """
    mapq_dist, mapq_orig_dist = defaultdict(int), defaultdict(int)
    ztz_dist = defaultdict(lambda: defaultdict(int))
    line = 0
    for ln in fh_sam:
        line += 1
        if ln[0] == '@':
            continue
        toks = ln.rstrip().split('\t')
        flags = int(toks[1])
        if (flags & 4) != 0 or (flags & 2048) != 0:
            continue
        pred = fh_pred.readline()
        assert len(pred) > 0
        predline, pred = pred.split(',')
        assert int(predline) == line, (predline, line)
        mapq_dist[int(round(float(pred)))] += 1
        mapq_orig_dist[int(toks[4])] += 1
        for tok in toks[12:]:
            if tok.startswith('ZT:Z:'):
                for i, ztok in enumerate(tok[5:].split(',')):
                    if ztok == 'NA':
                        pass
                    else:
                        ztz_dist[i][float(ztok)] += 1  # ZT:Z
                break
    return mapq_dist, mapq_orig_dist, ztz_dist


def pass1_fn(fn_sam, fn_pred):
    """ Histograms all the input and tandem alignments """
    with open(fn_sam) as fh_sam:
        with open(fn_pred) as fh_pred:
            return pass1_fh(fh_sam, fh_pred)

"""
Example: 10_26049747_26049846_0:0:0_0:0:0_100_100_0_3999999
offset is 0-based?
"""
_wgsimex_re = re.compile('(.+)_([^_]+)_([^_]+)_([^:]+):([^:]+):([^_]+)_([^:]+):([^:]+):([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^/]+).*')
#                           1   2       3       4       5       6       7       8       9       10      11      12      13


def name_is_extended_wgsim(nm):
    return _wgsimex_re.match(nm) is not None


def pos_from_extended_wgsim(name, mate2=False):
    res = _wgsimex_re.match(name)
    refid, fragst1, fragen1 = res.group(1), int(res.group(2))-1, int(res.group(3))-1
    len1, len2 = int(res.group(10)), int(res.group(11))
    flip = res.group(12) == '1'
    ln = len2 if mate2 else len1
    if flip == mate2:
        return refid, fragst1, True
    else:
        return refid, fragen1 - (ln-1), False

"""
Example: qsim!:GL000229.1:+:6005:100:u
pretty sure offset is 0-based
"""
_qsim_re = re.compile('qsim!:([^:]+):([+-]):([^:]+):.*')


def name_is_qsim(nm):
    return _qsim_re.match(nm) is not None


def pos_from_qsim(name):
    res = _qsim_re.match(name)
    return res.group(1), int(res.group(3)), res.group(2) == '+'


"""
Example: !h!chr9!118085975!+!50!0
offset is 0-based
"""
_hint_re = re.compile('!h!([^!]+)!([0-9]+)!([+-])!([0-9]+)!([0-9]+)')


def name_is_hint(name):
    return _hint_re.match(name) is not None


def pos_from_hint(name):
    res = _hint_re.match(name)
    return res.group(1), int(res.group(2)), res.group(3) == '+'


"""
Example:
hg38_50nt_mason1_unp.fastq.000999999 contig=chr9 haplotype=1 length=50 orig_begin=118085976 orig_end=118086026 snps=0 indels=0 haplotype_infix=ATGACTCTTGAAGCTGGGCGCAGTGGCTCATGCCTGTAATCCTAGCACTT edit_string=MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM strand=forward
"""
_mason_re = re.compile('[^ ]+ contig=([^ ]+) .* orig_begin=([^ ]+) .* strand=([fr]).*')


def name_is_mason1(name):
    return _mason_re.match(name) is not None


def pos_from_mason1(name):
    res = _hint_re.match(name)
    return res.group(1), int(res.group(2)), res.group(3) == 'f'


def same_pos(pos1, pos2, wiggle=30):
    """ Returns true when the two positions are basically the same """
    refid1, pos1, strand1 = pos1
    refid2, pos2, strand2 = pos2
    if refid1 != refid2 or strand1 != strand2:
        return False
    return abs(pos1 - pos2) < wiggle


def is_correct(toks, wiggle=30):
    """ Checks whether alignment, tokenized in toks, is correct """
    flags = int(toks[1])
    aligned_pos = (toks[2], int(toks[3])-1, (flags & 16) == 0)
    paired = (flags & 1) != 0
    mate2 = paired and (flags & 128) != 0
    if name_is_extended_wgsim(toks[0]):
        true_pos = pos_from_extended_wgsim(toks[0], mate2)
    elif name_is_qsim(toks[0]):
        true_pos = pos_from_qsim(toks[0])
    elif name_is_mason1(toks[0]):
        true_pos = pos_from_mason1(toks[0])
    elif name_is_hint(toks[0]):
        true_pos = pos_from_hint(toks[0])
    else:
        raise RuntimeError('Name was not formatted as expected: "%s"' % toks[0])
    return same_pos(true_pos, aligned_pos, wiggle=wiggle)


def percentileize(dist):
    last_v = 0
    cum = {}
    for k, v in sorted(dist.items()):
        assert k not in cum
        cum[k] = last_v + v / 2.0
        last_v += v
    pct_dict = {}
    for k in dist.keys():
        pct_dict[k] = (cum[k], float(cum[k])/last_v, dist[k])
    return pct_dict


def pass2_fh(fh_sam, fh_preds, mapq_dist, mapq_orig_dist, ztz_dist, ofh, ofh_orig, ofh_both):
    logging.info('  Calculating percentiles')
    mapq_pctile_dict = percentileize(mapq_dist)
    mapq_orig_pctile_dict = percentileize(mapq_orig_dist)
    ztz_pctile_dict = {x: percentileize(v) for x, v in ztz_dist.items()}
    incor_mapq, incor_mapq_orig, incor_both = [], [], []
    line = 0
    for ln in fh_sam:
        line += 1
        if ln[0] == '@':
            continue
        toks = ln.rstrip().split('\t')
        flags = int(toks[1])
        if (flags & 4) != 0 or (flags & 2048) != 0:
            continue
        pred = fh_preds.readline()
        assert len(pred) > 0
        predline, pred = pred.split(',')
        assert int(predline) == line
        if not is_correct(toks):
            mapq = int(round(float(pred)))
            mapq_orig = int(toks[4])  # TODO: make mapq_orig dist
            assert mapq in mapq_pctile_dict
            mapq_pctile = mapq_pctile_dict[mapq]
            mapq_orig_pctile = mapq_orig_pctile_dict[mapq_orig]
            ls = [mapq_pctile[1]]
            ls_orig = [mapq_orig_pctile[1]]
            ls_both = [mapq_pctile[1], mapq_orig_pctile[1]]
            for tok in toks[12:]:
                if tok.startswith('ZT:Z:'):
                    for i, ztok in enumerate(tok[5:].split(',')):
                        if ztok == 'NA':
                            ls.append('NA')
                            ls_orig.append('NA')
                            ls_both.append('NA')
                        else:
                            iztok = float(ztok)
                            assert iztok in ztz_pctile_dict[i]
                            ztz_pctile = ztz_pctile_dict[i][iztok]
                            ls.append(ztz_pctile[1])
                            ls_orig.append(ztz_pctile[1])
                            ls_both.append(ztz_pctile[1])
                    break
            incor_mapq.append(ls)
            incor_mapq_orig.append(ls_orig)
            incor_both.append(ls_both)
    logging.info('  Sorting and printing')
    for ls in sorted(incor_mapq, reverse=True):
        ofh.write(','.join(map(str, ls)) + '\n')
    for ls in sorted(incor_mapq_orig, reverse=True):
        ofh_orig.write(','.join(map(str, ls)) + '\n')
    for ls in sorted(incor_both, reverse=True):
        ofh_both.write(','.join(map(str, ls)) + '\n')


def pass2_fn(sam_fn, predictions_fn, mapq_dist, mapq_orig_dist, ztz_dist, ofn, ofn_orig, ofn_both):
    with open(ofn, 'w') as ofh:
        with open(ofn_orig, 'w') as ofh_orig:
            with open(ofn_both, 'w') as ofh_both:
                with open(sam_fn) as fh_sam:
                    with open(predictions_fn) as fh_preds:
                        pass2_fh(fh_sam, fh_preds, mapq_dist, mapq_orig_dist, ztz_dist, ofh, ofh_orig, ofh_both)


def write_dists_fh(mapq_dist, mapq_orig_dist, ztz_dist, ofh):
    for k, v in sorted(mapq_dist.items()):
        ofh.write('%d:%d ' % (int(k), v))
    ofh.write('\n')
    for k, v in sorted(mapq_orig_dist.items()):
        ofh.write('%d:%d ' % (k, v))
    ofh.write('\n')
    for k_outer, v_outer in sorted(ztz_dist.items()):
        ofh.write('ztz%d: ' % k_outer)
        for k, v in sorted(v_outer.items()):
            ofh.write('%d:%d ' % (k, v))
        ofh.write('\n')


def write_dists_fn(mapq_dist, mapq_orig_dist, ztz_dist, ofn):
    with open(ofn, 'w') as ofh:
        write_dists_fh(mapq_dist, mapq_orig_dist, ztz_dist, ofh)


def go():
    format_str = '%(asctime)s:%(levelname)s:%(message)s'
    level = logging.INFO
    logging.basicConfig(format=format_str, datefmt='%m/%d/%y-%H:%M:%S', level=level)

    logging.info('Pass 1')
    mapq_dist, mapq_orig_dist, ztz_dist = pass1_fn(sys.argv[1], sys.argv[2])
    logging.info('Writing distributions')
    write_dists_fn(mapq_dist, mapq_orig_dist, ztz_dist, sys.argv[3])
    logging.info('Pass 2')
    pass2_fn(sys.argv[1], sys.argv[2], mapq_dist, mapq_orig_dist, ztz_dist, sys.argv[4], sys.argv[5], sys.argv[6])
    logging.info('Done')


if __name__ == '__main__':
    go()
