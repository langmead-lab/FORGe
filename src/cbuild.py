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
    #include "cqf/cqf.h"
    #include "cqf/gqf.h"
    #include "bounter/bounter.h"
    """,
    libraries=[],
    sources=['cqf/cqf.c', 'cqf/gqf.c', 'bounter/bounter.c', 'bounter/murmur3.c'],
    # need -std=gnu99 for inline asm & loops w/ variables declared inside
    extra_compile_args=["-std=gnu99"])

ffibuilder.cdef("""
    struct flush_object * cqf_new(
        int ksize,
        int qbits);

    void cqf_delete(
        struct flush_object *o);

    int cqf_string_injest(
        const char *read,
        size_t read_len,
        struct flush_object *obj);

    int cqf_string_query(
        const char *read,
        size_t read_len,
        int64_t *count_array,
        size_t count_array_len,
        struct flush_object *obj);

    int cqf_est_fpr(
        int64_t *count_array,
        size_t count_array_len,
        struct flush_object *obj);

    void *bounter_new(
        int width,
        int depth);

    void bounter_delete(
        void *sketch);

    int bounter_string_injest(
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

    void qf_destroy(
        struct quotient_filter *qf,
        bool mem);
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
