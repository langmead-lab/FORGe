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


def same_pos(pos1, pos2, wiggle=30):
    """ Returns true when the two positions are basically the same """
    refid1, pos1, strand1 = pos1
    refid2, pos2, strand2 = pos2
    if refid1 != refid2 or strand1 != strand2:
        return False
    return abs(pos1 - pos2) < wiggle

def parse_label(label, rare_thresh=0.05):
    # Parse label and return # SNPs, # rare SNPs (frequency < rare_thresh), and # deleterious SNPs
    label = label.split(' ')
    snps, rare, deleterious, exome, conf, rep, alu = 0,0,0,0,0,0,0
    for l in label:
        if l[:6] == 'nsnps=':
            snps = int(l[6:])
        elif l[:6] == 'freqs=':
            freqs = l[6:]
            rare = 0
            if len(freqs) > 0:
                freqs = [float(f) for f in l[6:].split(',')]
                for f in freqs:
                    if f < rare_thresh:
                        rare += 1
        elif l[:4] == 'del=':
            deleterious = int(l[4:])
        elif l[:5] == 'conf=':
            c = int(l[5:])
            if c == 2:
                conf = 1
            else:
                conf = 0
        elif l[:6] == 'exome=':
            e = int(l[6:])
            if e == 2:
                exome = 1
            else:
                exome = 0
        elif l[:4] == 'rep=':
            r = int(l[4:])
            if r == 2:
                rep = 1
            else:
                rep = 0
        elif l[:4] == 'alu=':
            r = int(l[4:])
            if r == 2:
                alu = 1
            else:
                alu = 0
    return snps, rare, deleterious, exome, conf, rep, alu

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


def go(results_prefix, pct, desc):
    ncorrect, nincorrect, ntotal = 0, 0, 0
    ncorrect_snp, nincorrect_snp, ntotal_snp = [0], [0], [0]
    max_snps = 0

    rare_freq = 0.05
    ncorrect_rare, nincorrect_rare, ntotal_rare = [0,0], [0,0], [0,0]
    #max_rare = 0

    ncorrect_del, nincorrect_del, ntotal_del = [0], [0], [0]
    max_del = 0

    ncorrect_exome, nincorrect_exome, ntotal_exome = [0,0], [0,0], [0,0]
    ncorrect_conf, nincorrect_conf, ntotal_conf = [0,0], [0,0], [0,0]
    ncorrect_rep, nincorrect_rep, ntotal_rep = [0,0], [0,0], [0,0]
    ncorrect_alu, nincorrect_alu, ntotal_alu = [0,0], [0,0], [0,0]
    ncorrect_cent, nincorrect_cent, ntotal_cent, cent_bounds = 0, 0, 0, [39000000,71000000]

    #f_incorrect0 = open('incorrect0.sam', 'w')
    #f_incorrect1 = open('incorrect1.sam', 'w')

    for ln in sys.stdin:
        if ln[0] == '@':
            continue
        ln = ln.rstrip()
        toks = ln.split('\t')

        # Skip unaligned and secondary alignments
        if int(toks[1]) & 256:
            #print(ln + '\tZC:i:-1')
            continue

        ntotal += 1
        nsnps, nrare, ndel, exome, conf, rep, alu = parse_label(toks[0], rare_freq)
        if nsnps > max_snps:
            ncorrect_snp += [0] * (nsnps - max_snps)
            nincorrect_snp += [0] * (nsnps - max_snps)
            ntotal_snp += [0] * (nsnps - max_snps)
            max_snps = nsnps

        rare_id = -1
        if nsnps == 1:
            if nrare == 1:
                rare_id = 1
            else:
                rare_id = 0

        '''
        if nrare > max_rare:
            ncorrect_rare += [0] * (nrare - max_rare)
            nincorrect_rare += [0] * (nrare - max_rare)
            ntotal_rare += [0] * (nrare - max_rare)
            max_rare = nrare
        '''

        if ndel > max_del:
            ncorrect_del += [0] * (ndel - max_del)
            nincorrect_del += [0] * (ndel - max_del)
            ntotal_del += [0] * (ndel - max_del)
            max_del = ndel

        ntotal_snp[nsnps] += 1
        if rare_id >= 0:
            ntotal_rare[rare_id] += 1
        ntotal_del[ndel] += 1
        ntotal_exome[exome] += 1
        ntotal_conf[conf] += 1
        ntotal_rep[rep] += 1
        ntotal_alu[alu] += 1

        true_pos = pos_from_mason1(toks[0])
        if true_pos >= cent_bounds[0] and true_pos < cent_bounds[1]:
            ntotal_cent += 1

        if toks[1] == '4':
            continue

        if is_correct(toks, 30):
            #print(ln + '\tZC:i:1')
            ncorrect += 1
            ncorrect_snp[nsnps] += 1
            if rare_id >= 0:
                ncorrect_rare[rare_id] += 1
            ncorrect_del[ndel] += 1
            ncorrect_exome[exome] += 1
            ncorrect_conf[conf] += 1
            ncorrect_rep[rep] += 1
            ncorrect_alu[alu] += 1
            if true_pos >= cent_bounds[0] and true_pos < cent_bounds[1]:
                ncorrect_cent += 1
        else:
            #print(ln + '\tZC:i:0')
            nincorrect += 1
            #if nsnps == 0:
            #    f_incorrect0.write(ln+'\n')
            #elif nsnps == 1:
            #    f_incorrect1.write(ln+'\n')

            nincorrect_snp[nsnps] += 1
            if rare_id >= 0:
                nincorrect_rare[rare_id] += 1
            nincorrect_del[ndel] += 1
            nincorrect_exome[exome] += 1
            nincorrect_conf[conf] += 1
            nincorrect_rep[rep] += 1
            nincorrect_alu[alu] += 1
            if true_pos >= cent_bounds[0] and true_pos < cent_bounds[1]:
                nincorrect_cent += 1

    #f_incorrect0.close()
    #f_incorrect1.close()

    aligned = 100 * float(ncorrect + nincorrect) / ntotal
    correct = 100 * float(ncorrect) / (ncorrect + nincorrect)
    overall = 100 * float(ncorrect) / ntotal
    with open(results_prefix+'.tsv', 'a') as f:
        f.write('%s\t%s\t%f\t%f\t%f\n' % (pct, desc, aligned, correct, overall))

    with open(results_prefix+'.strat_snp.tsv', 'a') as f:
        for i in range(max_snps+1):
            tot = ntotal_snp[i]
            cor = ncorrect_snp[i]
            inc = nincorrect_snp[i]
            if tot > 0:
                aligned = float(cor+inc)/tot
                overall = float(cor)/tot
                if (cor+inc) == 0:
                    correct = -1
                else:
                    correct = float(cor)/(cor+inc)

                f.write('%s\t%d\t%d\t%s\t%f\t%f\t%f\n' % (pct, i, tot, desc, aligned, correct, overall))

    with open(results_prefix+'.strat_rare.tsv', 'a') as f:
        for i in range(2):
            tot = ntotal_rare[i]
            cor = ncorrect_rare[i]
            inc = nincorrect_rare[i]
            if tot > 0:
                aligned = float(cor+inc)/tot
                overall = float(cor)/tot
                if (cor+inc) == 0:
                    correct = -1
                else:
                    correct = float(cor)/(cor+inc)

                f.write('%s\t%d\t%d\t%s\t%f\t%f\t%f\n' % (pct, i, tot, desc, aligned, correct, overall))

    with open(results_prefix+'.strat_del.tsv', 'a') as f:
        for i in range(max_del+1):
            tot = ntotal_del[i]
            cor = ncorrect_del[i]
            inc = nincorrect_del[i]
            if tot > 0:
                aligned = float(cor+inc)/tot
                overall = float(cor)/tot
                if (cor+inc) == 0:
                    correct = -1
                else:
                    correct = float(cor)/(cor+inc)

                f.write('%s\t%d\t%d\t%s\t%f\t%f\t%f\n' % (pct, i, tot, desc, aligned, correct, overall))

    with open(results_prefix+'.strat_region.tsv', 'a') as f:
        tot = ntotal
        cor = ncorrect
        inc = nincorrect
        if tot > 0:
            aligned = float(cor+inc)/tot
            overall = float(cor)/tot
            if (cor+inc) == 0:
                correct = -1
            else:
                correct = float(cor)/(cor+inc)
            f.write('%s\t%s\t%d\t%s\t%f\t%f\t%f\n' % (pct, 'Total', tot, desc, aligned, correct, overall))

        i = 1
        tot = ntotal_exome[i]
        cor = ncorrect_exome[i]
        inc = nincorrect_exome[i]
        if tot > 0:
            aligned = float(cor+inc)/tot
            overall = float(cor)/tot
            if (cor+inc) == 0:
                correct = -1
            else:
                correct = float(cor)/(cor+inc)
            f.write('%s\t%s\t%d\t%s\t%f\t%f\t%f\n' % (pct, 'Exome', tot, desc, aligned, correct, overall))

        i = 0
        tot = ntotal_conf[i]
        cor = ncorrect_conf[i]
        inc = nincorrect_conf[i]
        if tot > 0:
            aligned = float(cor+inc)/tot
            overall = float(cor)/tot
            if (cor+inc) == 0:
                correct = -1
            else:
                correct = float(cor)/(cor+inc)
            f.write('%s\t%s\t%d\t%s\t%f\t%f\t%f\n' % (pct, 'NonConf', tot, desc, aligned, correct, overall))

        i = 1
        tot = ntotal_rep[i]
        cor = ncorrect_rep[i]
        inc = nincorrect_rep[i]
        if tot > 0:
            aligned = float(cor+inc)/tot
            overall = float(cor)/tot
            if (cor+inc) == 0:
                correct = -1
            else:
                correct = float(cor)/(cor+inc)
            f.write('%s\t%s\t%d\t%s\t%f\t%f\t%f\n' % (pct, 'Rep', tot, desc, aligned, correct, overall))

        i = 1
        tot = ntotal_alu[i]
        cor = ncorrect_alu[i]
        inc = nincorrect_alu[i]
        if tot > 0:
            aligned = float(cor+inc)/tot
            overall = float(cor)/tot
            if (cor+inc) == 0:
                correct = -1
            else:
                correct = float(cor)/(cor+inc)
            f.write('%s\t%s\t%d\t%s\t%f\t%f\t%f\n' % (pct, 'Alu', tot, desc, aligned, correct, overall))

        tot = ntotal_cent
        cor = ncorrect_cent
        inc = nincorrect_cent
        if tot > 0:
            aligned = float(cor+inc)/tot
            overall = float(cor)/tot
            if (cor+inc) == 0:
                correct = -1
            else:
                correct = float(cor)/(cor+inc)
            f.write('%s\t%s\t%d\t%s\t%f\t%f\t%f\n' % (pct, 'Cent', tot, desc, aligned, correct, overall))

if __name__ == '__main__':
    results_prefix = sys.argv[1]
    pct = sys.argv[2]
    desc = sys.argv[3]
    go(results_prefix, pct, desc)
