/*
 * Author: Ben Langmead
 * Contact: ben.langmead@gmail.com
 *
 * A C API for the CQF.  Required rewriting some C++ functions from Squeakr in
 * C, but changes are minor overall.  Intended for use with CFFI.  As of now,
 * only supports creating an empty CQF (qf_init), populating it with k-mers
 * from a longer string (string_injest), and querying it with k-mers from a
 * string (string_query).  It does not (yet) allow reading or writing CQFs to
 * disk, or injesting/querying from a file. 
 *
 * See:
 *
 * squeakr_c_api.h -- prototypes for string_injest and string_query
 *
 * squeakr_build.py -- builds the C code & C/Pyhton interface using CFFI
 * squeakr_query.py -- example of how to allocate, populate and query a CQF
 *                     with this API
 *
 * If you build with -DSQUEAKR_TEST_MAIN (use "make squeakr-test" target) you
 * can use the main function defined here.  Usage:
 *
 * ./squeakr-test <k-mer length> <qbits> <string to injest> \
 *                <query string 1> [<query string 2>...]
 *
 */

#include "squeakr_c_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : \
                                        (1ULL << (nbits)) - 1ULL)
#define QBITS_LOCAL_QF 16

enum DNA_MAP {DNA_C, DNA_A, DNA_T, DNA_G};  // A=1, C=0, T=2, G=3

/*return the integer representation of the base */
static inline uint8_t map_base(char base) {
	switch(toupper(base)) {
		case 'A': { return DNA_A; }
		case 'T': { return DNA_T; }
		case 'C': { return DNA_C; }
		case 'G': { return DNA_G; }
		default:  { return DNA_G+1; }
	}
}

/* Return the reverse complement of a base */
static inline int reverse_complement_base(int x) { return 3 - x; }

/* Calculate the revsese complement of a kmer */
static inline uint64_t reverse_complement(uint64_t kmer, uint32_t K) {
	uint64_t rc = 0;
	uint8_t base = 0;
	for (int i=0; i<K; i++) {
		base = kmer & 3ULL;
		base = reverse_complement_base(base);
		kmer >>= 2;
		rc |= base;
		rc <<= 2;
	}
	rc >>=2;
	return rc;
}

/* Compare the kmer and its reverse complement and return the result
 * Return true if the kmer is greater than or equal to its
 * reverse complement.
 * */
static inline bool compare_kmers(uint64_t kmer, uint64_t kmer_rev) {
	return kmer >= kmer_rev;
}

/* dump the contents of a local QF into the main QF */
static void dump_local_qf_to_main(struct flush_object *obj) {
	QFi local_cfi;
	if (qf_iterator(obj->local_qf, &local_cfi, 0)) {
		do {
			uint64_t key = 0, value = 0, count = 0;
			qfi_get(&local_cfi, &key, &value, &count);
			qf_insert(obj->main_qf, key, 0, count, true, true);
		} while (!qfi_next(&local_cfi));
		qf_reset(obj->local_qf);
	}
}

/**
 * Converts a string of "ATCG" to a uint64_t
 * where each character is represented by using only two bits
 */
static uint64_t str_to_int(const char *str, size_t len) {
	uint64_t strint = 0;
	for(size_t i = 0; i < len; i++) {
		uint8_t curr = 0;
		switch (str[i]) {
			case 'A': { curr = DNA_A; break; }
			case 'T': { curr = DNA_T; break; }
			case 'C': { curr = DNA_C; break; }
			case 'G': { curr = DNA_G; break; }
		}
		strint = strint | curr;
		strint = strint << 2;
	}
	return strint >> 2;
}

static uint64_t MurmurHash64A(const void * key,
                              int len,
                              unsigned int seed)
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;
	
	uint64_t h = seed ^ (len * m);
	
	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);
	
	while(data != end) {
		uint64_t k = *data++;
		
		k *= m;
		k ^= k >> r;
		k *= m;
		
		h ^= k;
		h *= m;
	}
	
	const unsigned char * data2 = (const unsigned char*)data;
	
	switch(len & 7) {
		case 7: h ^= (uint64_t)data2[6] << 48;
		case 6: h ^= (uint64_t)data2[5] << 40;
		case 5: h ^= (uint64_t)data2[4] << 32;
		case 4: h ^= (uint64_t)data2[3] << 24;
		case 3: h ^= (uint64_t)data2[2] << 16;
		case 2: h ^= (uint64_t)data2[1] << 8;
		case 1: h ^= (uint64_t)data2[0];
			h *= m;
	};
	
	h ^= h >> r;
	h *= m;
	h ^= h >> r;
	
	return h;
}

void string_injest(const char *read,
                   size_t read_len,
                   struct flush_object *obj,
                   int get_lock)
{
	const uint32_t seed = obj->main_qf->metadata->seed;
	const __uint128_t range = obj->main_qf->metadata->range;
	assert(obj->local_qf == NULL || obj->local_qf->metadata->range == range);
	assert(obj->local_qf == NULL || obj->local_qf->metadata->seed == seed);
	size_t left_chop = 0;
	do {
		if (read_len - left_chop < obj->ksize) {
			return; // start with the next read if length is smaller than K
		}
		uint64_t first = 0;
		uint64_t first_rev = 0;
		uint64_t item = 0;
		
		for(int i=0; i<obj->ksize; i++) { //First kmer
			uint8_t curr = map_base(read[left_chop + i]);
			if (curr > DNA_G) { // 'N' is encountered
				left_chop += (i+1);
				continue; // BTL: does it properly check if there's <k left?
			}
			first = first | curr;
			first = first << 2;
		}
		first = first >> 2;
		first_rev = reverse_complement(first, obj->ksize);
		item = compare_kmers(first, first_rev) ? first : first_rev;
		item = MurmurHash64A(((void*)&item), sizeof(item), seed);
		item %= range;

		/*
		 * first try and insert in the main QF.
		 * If lock can't be accuired in the first attempt then
		 * insert the item in the local QF.
		 */
		if (!qf_insert(obj->main_qf, item, 0, 1, get_lock, false)) {
			assert(obj->local_qf != NULL);
			qf_insert(obj->local_qf, item, 0, 1, false, false);
			obj->count++;
			// check of the load factor of the local QF is more than 50%
			if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
				dump_local_qf_to_main(obj);
				obj->count = 0;
			}
		}
		
		uint64_t next = (first << 2) & BITMASK(2*obj->ksize);
		uint64_t next_rev = first_rev >> 2;
		
		for(uint32_t i=obj->ksize; i < (read_len - left_chop); i++) { //next kmers
			uint8_t curr = map_base(read[i + left_chop]);
			if (curr > DNA_G) { // 'N' is encountered
				left_chop += (i+1);
				continue;
			}
			next |= curr;
			uint64_t tmp = reverse_complement_base(curr);
			tmp <<= (obj->ksize*2-2);
			next_rev = next_rev | tmp;
			item = compare_kmers(next, next_rev) ? next : next_rev;
			item = MurmurHash64A(((void*)&item), sizeof(item), seed);
			item %= range;
			
			/*
			 * first try and insert in the main QF.
			 * If lock can't be accuired in the first attempt then
			 * insert the item in the local QF.
			 */
			if (!qf_insert(obj->main_qf, item, 0, 1, get_lock, false)) {
				qf_insert(obj->local_qf, item, 0, 1, false, false);
				obj->count++;
				// check of the load factor of the local QF is more than 50%
				if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
					dump_local_qf_to_main(obj);
					obj->count = 0;
				}
			}
			next = (next << 2) & BITMASK(2*obj->ksize);
			next_rev = next_rev >> 2;
		}
	} while(false);
}

/**
 * Assumes no locking is needed, i.e. that no other thread
 * is trying to update the CQF.
 */
int string_query(const char *read_orig,
                 size_t read_len_orig,
                 int64_t *count_array,
                 size_t count_array_len,
                 struct flush_object *obj)
{
	const uint32_t seed = obj->main_qf->metadata->seed;
	const __uint128_t range = obj->main_qf->metadata->range;
	assert(obj->local_qf == NULL || obj->local_qf->metadata->range == range);
	assert(obj->local_qf == NULL || obj->local_qf->metadata->seed == seed);
	int64_t *cur_count = count_array;
	const char *read = read_orig;
	size_t read_len = read_len_orig;
	do {
		if(read_len < obj->ksize) {
			return (int)(cur_count - count_array);
		}
		uint64_t first = 0;
		uint64_t first_rev = 0;
		uint64_t item = 0;
		
		int do_continue = 0;
		for(int i = 0; i<obj->ksize; i++) {
			uint8_t curr = map_base(read[i]);
			if (curr > DNA_G) { // 'N' is encountered
				// append -1s for k-mers that include a non-ACGT 
				for(int j = 0; j <= i && j + obj->ksize <= read_len; j++) {
					*cur_count++ = -1;
					assert(cur_count - count_array <= count_array_len);
				}
				read += (i+1);
				read_len -= (i+1);
				do_continue = 1;
				break;
			}
			first = first | curr;
			first = first << 2;
		}
		if(do_continue) {
			continue;
		}
		first = first >> 2;
		first_rev = reverse_complement(first, obj->ksize);
		item = compare_kmers(first, first_rev) ? first : first_rev;
		item = MurmurHash64A(((void*)&item), sizeof(item), seed);
		item %= range;
		uint64_t count = qf_count_key_value(obj->main_qf, item, 0);
		*cur_count++ = count;
		assert(cur_count - count_array <= count_array_len);
		
		uint64_t next = (first << 2) & BITMASK(2*obj->ksize);
		uint64_t next_rev = first_rev >> 2;
		read += obj->ksize;
		assert(read_len >= obj->ksize);
		read_len -= obj->ksize;

		for(uint32_t i = 0; i < read_len; i++) { //next kmers
			uint8_t curr = map_base(read[i]);
			if (curr > DNA_G) { // 'N' is encountered
				for(int j = 0; j < obj->ksize; j++) {
					*cur_count++ = -1;
				}
				read++;
				read_len--;
				do_continue = 1;
				break;
			}
			next |= curr;
			uint64_t tmp = reverse_complement_base(curr);
			tmp <<= (obj->ksize*2-2);
			next_rev = next_rev | tmp;
			item = compare_kmers(next, next_rev) ? next : next_rev;
			item = MurmurHash64A(((void*)&item), sizeof(item), seed);
			item %= range;
			count = qf_count_key_value(obj->main_qf, item, 0);
			*cur_count++ = count;
			assert(cur_count - count_array <= count_array_len);
			next = (next << 2) & BITMASK(2*obj->ksize);
			next_rev = next_rev >> 2;
		}
		if(!do_continue) {
			break;
		}
	} while(1);
	return (int)(cur_count - count_array);
}

#ifdef SQUEAKR_TEST_MAIN
/* main method */
int main(int argc, char *argv[]) {
	QF cf;
	
	if(argc < 5) {
		fprintf(stderr, "Need at least 4 args\n");
		return 1;
	}
	
	int ksize = atoi(argv[1]);
	int qbits = atoi(argv[2]);
	
	fprintf(stderr, "Using ksize: %d\n", ksize);
	fprintf(stderr, "Using qbits: %d\n", qbits);

	struct flush_object obj;
	obj.local_qf = NULL;
	obj.main_qf = &cf;
	obj.ksize = ksize;
	obj.count = 0;
	
	int num_hash_bits = qbits+8;	// we use 8 bits for remainders in the main QF
	uint32_t seed = 2038074761;
	qf_init(&cf, (1ULL<<qbits), num_hash_bits, 0, true, "", seed);
	
	fprintf(stderr, "Building reference from: %s\n", argv[3]);
	string_injest(argv[3], strlen(argv[3]), &obj, false); // no locking
	
	for(size_t i = 4; i < argc; i++) {
		size_t len = strlen(argv[i]);
		size_t results_len = len - ksize + 1;
		int64_t *results = (int64_t*)malloc(sizeof(int64_t) * results_len);
		if(results == NULL) {
			fprintf(stderr, "ERROR: could not allocate results\n");
			return 1;
		}
		int nres = string_query(argv[i], len, results, results_len, &obj);
		assert(nres == (int)results_len);
		for(int j = 0; j < nres; j++) {
			fprintf(stderr, "[%d]: %lld\n", j, results[j]);
		}
		free(results);
	}
}
#endif
