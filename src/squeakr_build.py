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
    #include "threadsafe-gqf/gqf.h"
    #include "squeakr_c_api.h"
    """,
    libraries=[],
    sources=['squeakr_c_api.c', 'threadsafe-gqf/gqf.c'],
    # need -std=gnu99 for inline asm & loops w/ variables declared inside
    extra_compile_args=["-std=gnu99"])

ffibuilder.cdef("""
    struct quotient_filter {
        ...;
    };

    struct flush_object {
        struct quotient_filter* local_qf;
        struct quotient_filter* main_qf;
        uint32_t count;
        uint32_t ksize;
        uint64_t fp_buf[100];
    };

    struct flush_object * cqf_new(int ksize, int qbits);

    void cqf_delete(struct flush_object *o);

    int string_injest(
        const char *read,
        size_t read_len,
        struct flush_object *obj);

    int string_query(
        const char *read,
        size_t read_len,
        int64_t *count_array,
        size_t count_array_len,
        struct flush_object *obj);

    int est_fpr(
        int64_t *count_array,
        size_t count_array_len,
        struct flush_object *obj);

    void qf_destroy(struct quotient_filter *qf, bool mem);
""")

if __name__ == "__main__":
    import sys
    ffibuilder.compile(verbose=True, debug=True if '--debug' in sys.argv else None)
