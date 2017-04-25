#!/usr/bin/env python

"""
rep.py

Given a SAM file and a repeatmasker database (in the .fa.out format that
the repeatmasker.org site provides), write a new version of the SAM file
where each read is annotated according to which repeat families/
subfamilies it overlaps and by how many bases.

Look here for repeatmasker databases:
http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html

Depends on bx-tools Python module by James Taylor, available through
bioconda.
https://bioconda.github.io/recipes/bx-python/README.html

For GRCh37, the MHC region, containing HLA genes, is:
chr6:28,477,797-33,448,354
https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37

For GRCh37, the LRC region, containing KIR genes, is:
chr19:54,528,888-55,595,686
https://www.ncbi.nlm.nih.gov/grc/human/regions/LRC?asm=GRCh37.p13

Nextera exome capture probes:
http://support.illumina.com/sequencing/sequencing_kits/nextera-rapid-capture-exome-kit/downloads.html
http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed

Author: Ben Langmead
Date: 7/25/2012
Contact: langmea@cs.jhu.edu
"""

from __future__ import print_function
import argparse
import sys
import string
import re
import time
from collections import defaultdict
from bx.intervals.intersection import Interval
from bx.intervals.intersection import IntervalTree


def openex(fn, mode="rb"):
    if fn.endswith(".gz"):
        import gzip
        return gzip.open(fn, mode)
    elif fn.endswith(".bz2"):
        import bz2
        return bz2.BZ2File(fn, mode)
    else:
        return open(fn, mode)

parser = argparse.ArgumentParser(
    description='Analyze errors and read-level evidence in a bisulfite bam file.')

# Basic inputs
parser.add_argument(
    '--basename', metavar='path', type=str,
    required=True, help='Base name of sam file')
parser.add_argument(
    '--repeat-masker', metavar='path', type=str,
    required=True, help='RepeatMasker file(s) containing repeat intervals')
parser.add_argument(
    '--exome-bed', metavar='path', type=str,
    required=True, help='BED file giving exome capture regions')

# Annotate SAM
parser.add_argument(
    '--sam-input', metavar='path', dest='sam_in', type=str,
    required=False, help='.sam/.bam files containing alignments to annotate.  '
                         '--sam-input must also be set.')
parser.add_argument(
    '--sam-output', metavar='path', dest='sam_out', type=str,
    required=False, help='Output for annotated .sam files.  --sam-input must '
                         'also be set.')

# Basics
parser.add_argument(
    '--test', dest='test', action='store_const', const=True, default=False,
    help='Do unit tests')
parser.add_argument(
    '--sanity', dest='sanity', action='store_const', const=True, default=False,
    help='Perform extra sanity checks')
parser.add_argument(
    '--verbose', dest='verbose', action='store_const', const=True,
    default=False, help='Be talkative')
parser.add_argument(
    '--profile', dest='profile', action='store_const', const=True,
    default=False, help='Profile the mode using python cProfile module')

args = parser.parse_args()

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
_revcomp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'Z': 'Z'}


def revcomp(x): return x[::-1].translate(_revcomp_trans)


class Repeat:
    def __init__(self, ln):
        (self.swsc, self.pctdiv, self.pctdel, self.pctins, self.refid,
         self.ref_i, self.ref_f, self.ref_remain, self.orient, self.rep_nm,
         self.rep_cl, self.rep_prior, self.rep_i, self.rep_f, self.unk) = ln.split()
        self.ref_i, self.ref_f = int(self.ref_i), int(self.ref_f)

_basename_re = re.compile('hap([AB])_([a-z_]*)([0-9]+)(_r100)?')
_basename_re2 = re.compile('hap([AB])*')

def parse_basename(s):
    if 'auto_haps' in s:
        ma = _basename_re2.match(s)
        hap = ma.group(1)
        return hap, None, None, None
    else:
        ma = _basename_re.match(s)
        hap, rank, pct, r100 = ma.group(1), ma.group(2), ma.group(3), ma.group(4)
        return hap, rank, pct, r100

"""
Example:
hg38_50nt_mason1_unp.fastq.000999999 contig=chr9 haplotype=1 length=50 orig_begin=118085976 orig_end=118086026 snps=0 indels=0 haplotype_infix=ATGACTCTTGAAGCTGGGCGCAGTGGCTCATGCCTGTAATCCTAGCACTT edit_string=MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM strand=forward
"""
_mason_re = re.compile('[^ ]+ contig=([^ ]+) .* orig_begin=([^ ]+) orig_end=([^ ]+) .* strand=([fr]).*')


def name_is_mason1(name):
    return _mason_re.match(name) is not None


def pos_from_mason1(name):
    assert name is not None
    res = _mason_re.match(name)
    assert res is not None, name
    return res.group(1), int(res.group(2)), int(res.group(3)), res.group(4) == 'f'


def same_pos(pos1, pos2, wiggle=30):
    """ Returns true when the two positions are basically the same """
    refid1, pos1, _, strand1 = pos1
    refid2, pos2, strand2 = pos2
    if refid1 != refid2 or strand1 != strand2:
        return False
    return abs(pos1 - pos2) < wiggle


def is_correct(toks, wiggle=30):
    """ Checks whether alignment, tokenized in toks, is correct """
    flags = int(toks[1])
    aligned_pos = (toks[2], int(toks[3])-1, (flags & 16) == 0)
    true_pos = pos_from_mason1(toks[0])
    return same_pos(true_pos, aligned_pos, wiggle=wiggle)


def go():
    nrep_fns, nreps, nrepnt, nsamrecs, nsamolap = 0, 0, 0, 0, 0
    rep_itree, km_itree, ex_itree = {}, {}, defaultdict(IntervalTree)
    recs = defaultdict(int)

    # MHC: chr6:28,477,797-33,448,354
    km_itree['6'] = IntervalTree()
    km_itree['6'].insert_interval(Interval(28477797, 33448354, value='MHC'))

    # KIR: chr19:54,528,888-55,595,686
    km_itree['19'] = IntervalTree()
    km_itree['19'].insert_interval(Interval(54528888, 55595686, value='KIR'))

    #
    # First, we scan through all the RepeatMasker files and extract the
    # information we need for annotating SAM and/or reporting repeat
    # sequences.
    #
    hd_re = re.compile('^\s*[0-9]')
    with openex(args.repeat_masker) as ifh:
        nrep_fns += 1

        print("Processing RepeatMasker file %s ..." % args.repeat_masker, file=sys.stderr)
        print(time.ctime(), file=sys.stderr)
        for ln in ifh:
            # Skip header lines
            if not hd_re.match(ln):
                continue

            # Parse fields; see:
            #  http://www.repeatmasker.org/webrepeatmaskerhelp.html
            # 1. SW score
            # 2. Percent divergence
            # 3. Percent deletions
            # 4. Percent insertions
            # 5. Query sequence (reference genome sequence)
            # 6. Position begin
            # 7. Position end
            # 8. No of bases in query past end of match
            # 9. C = match is with the complement of the consensus
            # 10. Name of matching interspersed repeat
            # 11. Class of the repeat
            # 12. # bases in repeat prior to beginning of the match
            # 13. Starting position of match in repeat consensus
            # 14. Ending position of match in repeat consensus
            # 15. ????
            if args.sanity:
                if not len(ln.split()) == 15:
                    raise ValueError("Expected 15 fields:\n" + ln)
            rp = Repeat(ln)
            nreps += 1
            nrepnt += (rp.ref_f - rp.ref_i)

            if rp.refid not in rep_itree:
                rep_itree[rp.refid] = IntervalTree()
            iv = Interval(rp.ref_i, rp.ref_f, value=rp.rep_cl)
            rep_itree[rp.refid].insert_interval(iv)

        with openex(args.exome_bed) as ifh:
            print("Processing Exome bed file %s ..." % args.exome_bed, file=sys.stderr)
            print(time.ctime(), file=sys.stderr)
            for ln in ifh:
                toks = ln.rstrip().split()
                iv = Interval(int(toks[1]), int(toks[2]), value="exome")
                ex_itree[toks[0]].insert_interval(iv)

    # Finished scanning RepeatMasker files

    # Scan through .sam/.bam files and annotate them with information about
    # what repeats they overlap
    print("Processing SAM/BAM file %s ..." % args.sam_in, file=sys.stderr)
    print(time.ctime(), file=sys.stderr)
    ifh = sys.stdin
    if args.sam_in is not None:
        ifh = openex(args.sam_in)
    for ln in ifh:
        if ln[0] == '@':
            continue
        ln = ln.rstrip()
        toks = ln.split('\t')
        refid, ref_i, ref_f, _ = pos_from_mason1(toks[0])
        rec = []
        for tree in [rep_itree, ex_itree]:
            max_olap, max_lab = 0, None
            if refid in tree:
                for olap in tree[refid].find(ref_i, ref_f):
                    nsamolap += 1
                    en, st = min(olap.end, ref_f), max(olap.start, ref_i)
                    if en - st > max_olap:
                        max_olap = min(en - st, len(toks[9]))
                        max_lab = olap.value
            rec.append(str(max_lab))
        correctness = -1
        if toks[1] != '4':
            correctness = 1 if is_correct(toks) else 0
        rec.append(str(correctness))
        recs[','.join(rec)] += 1

    ifh.close()

    hap, rank, pct, r100 = parse_basename(args.basename)
    if rank is None:
        prefix = '%s,auto_hap,100,%s' % (hap, '35' if r100 is None else '100')
    else:
        prefix = '%s,%s,%s,%s' % (hap, rank, pct, '35' if r100 is None else '100')

    ofh = sys.stdout
    if args.sam_out is not None:
        ofh = open(args.sam_out, 'wb')

    for k, v in recs.items():
        ofh.write('%s,%s,%d\n' % (prefix, k, v))

    if args.sam_out is not None:
        ofh.close()

    print("# repeat files: %d" % nrep_fns, file=sys.stderr)
    print("# repeats: %d" % nreps, file=sys.stderr)
    print("# repeat nucleotides: %d" % nrepnt, file=sys.stderr)
    print("# SAM records: %d" % nsamrecs, file=sys.stderr)
    print("# repeat/SAM overlaps: %d" % nsamolap, file=sys.stderr)


if not args.test:
    if args.profile:
        import cProfile

        cProfile.run('go()')
    else:
        go()

if args.test:
    import unittest


    class TestCases(unittest.TestCase):
        def test_blah(self):
            pass


    unittest.main(argv=[sys.argv[0]])
    sys.exit()
