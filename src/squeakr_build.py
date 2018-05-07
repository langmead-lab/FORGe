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
    sources=['squeakr_c_api.c', 'threadsafe-gqf/gqf.c'])

ffibuilder.cdef("""
    struct quotient_filter {
        ...;
    };

    struct flush_object {
        struct quotient_filter* local_qf;
        struct quotient_filter* main_qf;
        uint32_t count;
        uint32_t ksize;
    };

    void string_injest(const char *read,
                       size_t read_len,
                       struct flush_object *obj,
                       int get_lock);

    int string_query(const char *read,
                     size_t read_len,
                     int64_t *count_array,
                     size_t count_array_len,
                     struct flush_object *obj);

    void qf_init(struct quotient_filter *qf,
                 uint64_t nslots,
                 uint64_t key_bits,
                 uint64_t value_bits,
                 bool mem,
                 const char *path,
                 uint32_t seed);
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
