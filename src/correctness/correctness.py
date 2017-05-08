#! /usr/bin/env python

"""
Given a SAM file with read names in Mason 1 format, Qsim tandem read format,
Extended wgsim format, or Bowtie 2 --hints format, measure how many are
correct.
"""

from __future__ import print_function
import sys
import re


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
_hint_seed_re = re.compile('!h!([^!]+)!([0-9]+)!([+-])!([0-9]+)!([0-9]+)')
_hint_ival_re = re.compile('!h!([^!]+)!([0-9]+)!([0-9]+)!([0-9]+)!([0-9]+)')


def name_is_hint(name):
    return _hint_seed_re.match(name) or _hint_ival_re.match(name)


def pos_from_hint(name):
    res = _hint_seed_re.match(name)
    if res is not None:
        return res.group(1), int(res.group(2)), res.group(3) == '+'
    raise RuntimeError("Can't really check interval hints since strand and exact offset are missing")


"""
Example:
hg38_50nt_mason1_unp.fastq.000999999 contig=chr9 haplotype=1 length=50 orig_begin=118085976 orig_end=118086026 snps=0 indels=0 haplotype_infix=ATGACTCTTGAAGCTGGGCGCAGTGGCTCATGCCTGTAATCCTAGCACTT edit_string=MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM strand=forward
"""
_mason_re = re.compile('[^ ]+ contig=([^ ]+) .* orig_begin=([^ ]+) .* strand=([fr]).*')


def name_is_mason1(name):
    return _mason_re.match(name) is not None


def pos_from_mason1(name):
    res = _mason_re.match(name)
    return res.group(1), int(res.group(2)), res.group(3) == 'f'


def same_pos(pos1, pos2, wiggle=0):
    """ Returns true when the two positions are basically the same """
    refid1, pos1, strand1 = pos1
    refid2, pos2, strand2 = pos2
    if refid1 != refid2 or strand1 != strand2:
        return False
    return abs(pos1 - pos2) < wiggle

def parse_label(label):
    # Parse label and return number of SNPs and their frequencies
    label = label.split(' ')
    snps = 0
    for l in label:
        if l[:6] == 'nsnps=':
            snps = int(l[6:])
        elif l[:6] == 'freqs=':
            freqs = l[6:]
            if len(freqs) > 0:
                freqs = [float(f) for f in l[6:].split(',')]
            else:
                freqs = []
    if snps == 0:
        return 0, []
    if not snps == len(freqs):
        print('Error!')
        print('%d SNPs' % snps)
        print('Freqs: ' + str(freqs))
        exit()
    return snps, freqs

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


def go():
    ncorrect, nincorrect = 0, 0
    ncorrect_by_snp, nincorrect_by_snp, ntotal_by_snp = [0], [0], [0]
    max_nsnps = 0

    rare_freq = 0.05
    ncorrect_by_rare, nincorrect_by_rare, ntotal_by_rare = [0], [0], [0]
    max_rare = 0

    for ln in sys.stdin:
        if ln[0] == '@':
            continue
        ln = ln.rstrip()
        toks = ln.split('\t')

        nsnps, freqs = parse_label(toks[0])
        if nsnps > max_nsnps:
            ncorrect_by_snp += [0] * (nsnps - max_nsnps)
            nincorrect_by_snp += [0] * (nsnps - max_nsnps)
            ntotal_by_snp += [0] * (nsnps - max_nsnps)
            max_nsnps = nsnps
        ntotal_by_snp[nsnps] += 1

        nrare = 0
        for f in freqs:
            if f < rare_freq:
                nrare += 1
        if nrare > max_rare:
            ncorrect_by_rare += [0] * (nrare - max_rare)
            nincorrect_by_rare += [0] * (nrare - max_rare)
            ntotal_by_rare += [0] * (nrare - max_rare)
            max_rare = nrare
        ntotal_by_rare[nrare] += 1

        if toks[1] == '4':
            #print(ln + '\tZC:i:-1')
            continue
        if is_correct(toks, 30):
            #print(ln + '\tZC:i:1')
            ncorrect += 1
            ncorrect_by_snp[nsnps] += 1
            ncorrect_by_rare[nrare] += 1
        else:
            #print(ln + '\tZC:i:0')
            nincorrect += 1
            nincorrect_by_snp[nsnps] += 1
            nincorrect_by_rare[nrare] += 1
    print('correct=%d (%0.4f%%)' % (ncorrect, ncorrect*100.0/(ncorrect+nincorrect)), file=sys.stderr)

    print('')
    print('Stratified by # SNPs:')
    print('# reads\tAligned\tCorrect\tOverall Correct')
    for i in range(max_nsnps+1):
        t = ntotal_by_snp[i]
        c = ncorrect_by_snp[i]
        i = nincorrect_by_snp[i]
        if t == 0:
            acc = 0
            corr_overall = 0
        else:
            acc = float(c+i)/t
            corr_overall = float(c)/t
        if (c+i) == 0:
            corr = 0
        else:
            corr = float(c)/(c+i)

        print('%d\t%f\t%f\t%f' % (t, acc, corr, corr_overall))

    '''
    print('')
    print('Stratified by # rare (< %f) SNPs:' % rare_freq)
    print('# reads\tAligned\tCorrect\tOverall Correct')
    for i in range(max_rare+1):
        t = ntotal_by_rare[i]
        c = ncorrect_by_rare[i]
        i = nincorrect_by_rare[i]
        if t == 0:
            acc = 0
            corr_overall = 0
        else:
            acc = float(c+i)/t
            corr_overall = float(c)/t
        if (c+i) == 0:
            corr = 0
        else:
            corr = float(c)/(c+i)

        print('%d\t%f\t%f\t%f' % (t, acc, corr, corr_overall))
    '''

if __name__ == '__main__':
    go()
