/*
 * Author: Ben Langmead
 * Contact: ben.langmead@gmail.com
 *
 * A C API for the CQF.
 */

#ifndef _SQUEAKR_C_API_H_
#define _SQUEAKR_C_API_H_

#include "gqf.h"

#define FP_BUF_ELTS 100

struct flush_object {
	QF *local_qf;
	QF *main_qf;
	uint32_t count;
	uint32_t ksize;

	// Holds true counts for elements that hash to numbers < FP_BUF_ELTS.  For
	// estimating false positive rate.  It's not quite right to have this in
	// the flush_object, since this is per-thread.  Better to add this to the
	// QF structure itself.
	uint64_t fp_buf[FP_BUF_ELTS];
};

/**
 * Allocate and initialize a CQF and flush_object struct
 */
struct flush_object * cqf_new(int ksize, int qbits);

/**
 * Free the CQF and flush_object struct.
 */
void cqf_delete(struct flush_object *o);

/**
 * Injest all the k-mers in a string into a CQF.
 */
int cqf_string_injest(
    const char *read,
    size_t read_len,
    struct flush_object *obj);

/**
 * Populate int array with the counts of all the k-mers in the given string.
 */
int cqf_string_query(
    const char *read,
    size_t read_len,
    int64_t *count_array,
    size_t count_array_len,
    struct flush_object *obj);

/**
 * Put false positive rate information into the given count_array.  Returns up
 * to min(FP_BUF_ELTS, range) pairs of uint64_ts.  Each pair consists of an
 * observed count and a true count.
 */
int cqf_est_fpr(
    int64_t *count_array,
    size_t count_array_len,
    struct flush_object *obj);

#endif /* _SQUEAKR_C_API_H_ */
