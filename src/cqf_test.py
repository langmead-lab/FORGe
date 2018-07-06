#!/usr/bin/env python

"""squeakr

Usage:
  squeakr fpr-test <keys> [options]

Options:
  --seed <int>     pseudo-random seed [default: 777].
  -q <int>         # quotient bits [default: 10].
  -k <int>         k-mer length [default: 20].
  -h, --help       Show this screen.
"""

"""
Author: Ben Langmead
Contact: ben.langmead@gmail.com

Example Python script with functions for allocating, populating and querying a
CQF.  Uses CFFI.  Run squeakr_build.py first.
"""

from _api import ffi, lib
import math
import random
import logging
from collections import Counter
from docopt import docopt


def cqf_est_fpr(db):
    """
    Query the given CQF with k-mers from the given string.  Return a C array
    of int64_t giving the count for each k-mer.  The C array is not the same
    as a Python list, but could be immediately converted to one using
    list(...) if desired.
    """
    results_len = 100 # lib.FP_BUF_ELTS*2
    results = ffi.new('int64_t[]', results_len)
    lib.cqf_est_fpr(results, results_len, db)
    return results, results_len


def test_simple_query_1():
    db = lib.cqf_new(4, 10)
    lib.cqf_string_injest(b'ACGT', 4, db)
    results = ffi.new('int64_t[]', 1)
    nresults = lib.cqf_string_query(b'ACGT', 4, results, 1, db)
    assert 1 == nresults
    assert 1 == results[0]
    lib.cqf_delete(db)


def test_query_with_ns():
    k = 4
    db = lib.cqf_new(k, 10)
    to_injest = b'ACGTACGT'
    lib.cqf_string_injest(to_injest, len(to_injest), db)
    to_query = b'ACGTNACGTNACGT'
    #           01234567890123
    stlen = len(to_query)
    results_len = stlen - k + 1
    results = ffi.new('int64_t[]', results_len)
    nresults = lib.cqf_string_query(to_query, stlen, results, results_len, db)
    assert nresults == results_len
    assert 2 == results[0]
    assert -1 == results[1]
    assert -1 == results[2]
    assert -1 == results[3]
    assert -1 == results[4]
    assert 2 == results[5]
    assert -1 == results[6]
    assert -1 == results[7]
    assert -1 == results[8]
    assert -1 == results[9]
    assert 2 == results[10]
    lib.cqf_delete(db)


def fpr_test(nkeys, seed, qbits, ksize, nchunks=None):
    db = lib.cqf_new(ksize, qbits)
    print(db)
    if nchunks is None:
        nchunks = int(math.sqrt(nkeys))+1
    nchunks = min(nchunks, nkeys)
    chunk_sz = int(nkeys / nchunks)
    so_far = 0
    random.seed(seed)
    fpr_results = ffi.new('int64_t[]', 200)
    for i in range(nchunks):
        st = b''.join([random.choice([b'A', b'C', b'G', b'T']) for _ in range(ksize + chunk_sz - 1)])
        print(st)
        lib.cqf_string_injest(st, len(st), db)
        print('Injested keys %d -- %d' % (so_far, so_far + chunk_sz))
        so_far += chunk_sz
        lib.est_fpr(fpr_results, 100, db)
        cnt = Counter(map(lambda x: x[1] - x[0], ((fpr_results[i], fpr_results[i+1]) for i in range(0, 100, 2))))
        print(cnt)
    lib.cqf_delete(db)


if __name__ == '__main__':
    args = docopt(__doc__)
    seed = int(args['--seed'])
    qbits = int(args['-q'])
    ksize = int(args['-k'])
    logging.basicConfig()
    if args['fpr-test']:
        fpr_test(int(args['<keys>']), seed, qbits, ksize)
