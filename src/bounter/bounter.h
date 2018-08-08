#include <stdint.h>

// TODO: avoid void *; here it's used because bounter doesn't have header
// files we can include to get the proper struct definitions

extern void *bounter_new(int width, int depth);

extern void bounter_delete(void *sketch);

/**
 * Given a string, extract every k-mer, canonicalize, and add to the QF.
 */
extern int bounter_string_injest(
    void *sketch,            // bounter sketch
    int k,                   // k-mer length
    const char *read,        // string whose k-mers to add
    size_t read_len);        // length of string

/**
 * Assumes no locking is needed, i.e. that no other thread
 * is trying to update the CQF.
 */
int bounter_string_query(
    void *sketch,            // bounter sketch
    int k,                   // k-mer length
    const char *read_orig,   // query string
    size_t read_len_orig,    // length of query
    int64_t *count_array,    // result array
    size_t count_array_len); // length of result array
