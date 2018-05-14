#!/usr/bin/env python

"""
Author: Ben Langmead
Contact: ben.langmead@gmail.com

Example Python script with functions for allocating, populating and querying a
CQF.  Uses CFFI.  Run squeakr_build.py first.
"""

from _api import ffi, lib


def cqf_new_db(ksize, qbits, seed, create_local=False):
    """
    Returns a new "db" (basically the `flush_object` structure from Squeakr)
    which can be used to represent a CQF, its metadata & state.  Though we
    don't use it for the (yet), it could also be used to hold per-thread state
    (i.e. a local CQF) which could be useful in scenarios where a writer is
    running concurrently with other readers or writers.
    """
    db = ffi.new('struct flush_object*')
    db.main_qf = lib.qf_init_simple(1 << qbits, qbits + 8, 0, seed)
    db.local_qf = ffi.NULL
    if create_local:
        db.local_qf = lib.qf_init_simple(1 << qbits, qbits + 8, 0, seed)
    db.ksize = ksize
    db.count = 0
    return db


def cqf_injest(db, s):
    return lib.string_injest(s, len(s), db, 0)


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
    if False:
        db = cqf_new_db(ksize, qbits, seed, False)
        to_injest = 'ACGTACGT'
        lib.string_injest(to_injest, len(to_injest), db, 0)
        to_query = 'ACGTNACGTNACGT'
        #           01234567890123
        results = cqf_query(db, to_query)
        #for i, res in enumerate(results):
        #    print((i, res))
    else:
        db = ffi.new('struct flush_object*')
        db.main_qf = lib.qf_init_simple(1 << qbits, qbits + 8, 0, seed)
        db.local_qf = ffi.NULL
        db.ksize = ksize
        db.count = 0

        to_injest = 'ACGTACGT'
        lib.string_injest(to_injest, len(to_injest), db, 0)

        to_query = 'ACGTNACGTNACGT'
        stlen = len(to_query)
        results_len = stlen - db.ksize + 1
        results = ffi.new('int64_t[]', results_len)
        nresults = lib.string_query(to_query, stlen, results, results_len, db)
        assert nresults == results_len
