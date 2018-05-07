/*
 * Author: Ben Langmead
 * Contact: ben.langmead@gmail.com
 *
 * A C API for the CQF.
 */

#ifndef _SQUEAKR_C_API_H_
#define _SQUEAKR_C_API_H_

#include "threadsafe-gqf/gqf.h"

struct flush_object {
	QF *local_qf;
	QF *main_qf;
	uint32_t count;
	uint32_t ksize;
};

/**
 * Injest all the k-mers in a string into a CQF.
 */
void string_injest(const char *read,
                   size_t read_len,
                   struct flush_object *obj,
                   int get_lock);

/**
 * Populate int array with the counts of all the k-mers in the given string.
 */
int string_query(const char *read,
                 size_t read_len,
                 int64_t *count_array,
                 size_t count_array_len,
                 struct flush_object *obj);

#endif /* _SQUEAKR_C_API_H_ */
