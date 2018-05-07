#!/usr/bin/env python

"""
Author: Ben Langmead
Contact: ben.langmead@gmail.com

Example Python script with functions for allocating, populating and querying a
CQF.  Uses CFFI.  Run squeakr_build.py first.
"""

from _api import ffi, lib


def cqf_new_filter(qbits, seed):
    """
    Returns a new counting quotient filter with the given number of quotient
    bits and the given seed. 
    """
    qf = ffi.new('struct quotient_filter*')
    lib.qf_init(qf, 1 << qbits, qbits+8, 0, True, "", seed)
    return qf


def cqf_new_db(ksize, qbits, seed, create_local=False):
    """
    Returns a new "db" (basically the `flush_object` structure from Squeakr)
    which can be used to represent a CQF, its metadata & state.  Though we
    don't use it for the (yet), it could also be used to hold per-thread state
    (i.e. a local CQF) which could be useful in scenarios where a writer is
    running concurrently with other readers or writers.
    """
    qf = cqf_new_filter(qbits, seed)
    qf_local = ffi.NULL
    if create_local:
        qf_local = cqf_new_filter(qbits, seed)
    flush_object = ffi.new('struct flush_object*')
    flush_object.ksize = ksize
    flush_object.count = 0
    flush_object.main_qf = qf
    flush_object.local_qf = qf_local
    return flush_object


def cqf_injest(db, s):
    lib.string_injest(s, len(s), db, 0)


def cqf_query(db, st):
    """
    Query the given CQF with k-mers from the given string.  Return a C array
    of int64_t giving the count for each k-mer.  The C array is not the same
    as a Python list, but could be immediately converted to one using
    list(...) if desired.
    """
    stlen = len(st)
    results_len = stlen - db.ksize + 1
    results = ffi.new('int64_t[]', results_len)
    nresults = lib.string_query(st, stlen, results, results_len, db)
    assert nresults == results_len
    return results


if __name__ == '__main__':
    seed = 777
    qbits = 10
    ksize = 4
    db = cqf_new_db(ksize, qbits, seed, False)
    to_injest = 'ACGTACGT'
    lib.string_injest(to_injest, len(to_injest), db, 0)
    to_query = 'ACGTNACGTNACGT'
    #           01234567890123
    results = cqf_query(db, to_query)
    for i, res in enumerate(results):
        print((i, res))
