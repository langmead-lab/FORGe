#!/usr/bin/env python

"""
Author: Ben Langmead
Contact: ben.langmead@gmail.com

CFFI generator file.  Used to build object files that can then be used from
python after "from _api import ffi, lib".
"""

from cffi import FFI
ffibuilder = FFI()

ffibuilder.set_source("_api", r"""
    #include "bounter/bounter.h"
    """,
    libraries=[],
    sources=['bounter/bounter.c', 'bounter/murmur3.c'],
    # need -std=gnu99 for inline asm & loops w/ variables declared inside
    extra_compile_args=["-std=gnu99"])

ffibuilder.cdef("""
    void *bounter_new(
        int width,
        int depth);

    void bounter_delete(
        void *sketch);

    int bounter_string_ingest(
        void *sketch,
        int k,
        const char *read,
        size_t read_len);

    int bounter_string_query(
        void *sketch,
        int k,
        const char *read,
        size_t read_len,
        int64_t *count_array,
        size_t count_array_len);

    int bounter_tsv_ingest(
        void *sketch,
        int k,
        const char *filename);
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
