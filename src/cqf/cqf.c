/*
 * Author: Ben Langmead
 * Contact: ben.langmead@gmail.com
 *
 * A C API for the CQF.  Required rewriting some C++ functions from Squeakr in
 * C, but changes are minor overall.  Intended for use with CFFI.  As of now,
 * only supports creating an empty CQF (qf_init), populating it with k-mers
 * from a longer string (cqf_string_injest), and querying it with k-mers from a
 * string (cqf_string_query).  It does not (yet) allow reading or writing CQFs to
 * disk, or injesting/querying from a file. 
 *
 * See:
 *
 * squeakr_c_api.h -- prototypes for cqf_string_injest and cqf_string_query
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

#include "cqf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : \
                                        (1ULL << (nbits)) - 1ULL)
#define QBITS_LOCAL_QF 16

/*return the integer representation of the base */
static inline uint8_t map_base(char base) {
	switch(toupper(base)) {
		case 'A': { return 0; }
		case 'T': { return 3; }
		case 'C': { return 1; }
		case 'G': { return 2; }
		default:  { return 4; }
	}
}

#define NON_ACGT(b) (b < 0 || b > 3)

/* Return the reverse complement of a base */
static inline int reverse_complement_base(int x) {
	assert(x < 4);
	assert(x >= 0);
	return 3 ^ x;
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

/**
 * Allocate and initialize a CQF and flush_object struct
 */
struct flush_object * cqf_new(int ksize, int qbits) {
	struct flush_object *obj = malloc(sizeof(struct flush_object));
	obj->local_qf = NULL;
	obj->main_qf = malloc(sizeof(QF));
	obj->ksize = ksize;
	obj->count = 0;
	memset(obj->fp_buf, 0, sizeof(uint32_t)*FP_BUF_ELTS);
	
	int num_hash_bits = qbits+8;	// we use 8 bits for remainders in the main QF
	uint32_t seed = 2038074761;
	qf_init(obj->main_qf, (1ULL << qbits), num_hash_bits, 0, true, "", seed);
	return obj;
}

/**
 * Allocate and initialize a CQF and flush_object struct
 */
void cqf_delete(struct flush_object *o) {
	if(o->local_qf != NULL) {
		free(o->local_qf);
		o->local_qf = NULL;
	}
	if(o->main_qf != NULL) {
		free(o->main_qf);
		o->main_qf = NULL;
	}
	free(o);
}

typedef struct {
	uint64_t u64s[4];
} b256;

/**
 * Set to all 0s.
 */
static void b256_clear(b256 *b) {
	memset(b, 0, 32);
}

/**
 * Left-shift by 2 bits and mask so that DNA string is no longer than k.  Put
 * result in dst.
 */
static void b256_lshift_and_mask(b256 *b, size_t k) {
	k *= 2;
	int word = (int)(k >> 6);
	for(int i = 3; i > 0; i--) {
		b->u64s[i] = (b->u64s[i] << 2) | (b->u64s[i-1] >> 62);
	}
	b->u64s[0] <<= 2;
	b->u64s[word] &= BITMASK(k & 63);
}

/**
 * Left-shift by 2 bits, in place.
 */
static void b256_lshift(b256 *b) {
	for(int i = 3; i > 0; i--) {
		b->u64s[i] = (b->u64s[i] << 2) | (b->u64s[i-1] >> 62);
	}
	b->u64s[0] <<= 2;
}

/**
 * Right-shift by 2 bits, in place.
 */
static void b256_rshift_2(b256 *b) {
	for(int i = 0; i < 3; i++) {
		b->u64s[i] = (b->u64s[i] >> 2) | (b->u64s[i+1] << 62);
	}
	b->u64s[3] >>= 2;
}

/**
 * OR the given value into the low bits of b, in place.
 */
static void b256_or_low(b256 *b, uint8_t val) {
	b->u64s[0] |= val;
}

/**
 * OR the given value into given bitpair of b, in place.
 */
static void b256_or_at(b256 *b, uint8_t val, int k) {
	assert(k < 128);
	assert(val < 4);
	k *= 2;
	size_t word = k >> 6, bitoff = k & 63;
	b->u64s[word] |= ((uint64_t)val << bitoff);
}

/**
 * Reverse complement the length-k DNA string in src, putting result in dst.
 */
static void b256_revcomp(const b256 *src, b256 *dst, size_t k) {
	// presumably *dst is all 0s
	assert(k > 0);
	size_t lo_word = 0, lo_bitoff = 0;
	size_t hi_word = ((k-1) * 2) >> 6;
	size_t hi_bitoff = ((k-1) * 2) & 63;
	while(1) {
		assert(lo_word >= 0);
		assert(lo_word < 4);
		assert(hi_word >= 0);
		assert(hi_word < 4);
		assert(lo_bitoff >= 0);
		assert(lo_bitoff < 64);
		assert(hi_bitoff >= 0);
		assert(hi_bitoff < 64);
		int base = (src->u64s[lo_word] >> lo_bitoff) & 3;
		uint64_t rcbase = reverse_complement_base(base);
		assert(rcbase < 4);
		dst->u64s[hi_word] |= rcbase << hi_bitoff;
		if(lo_bitoff == 62) {
			lo_bitoff = 0;
			lo_word++;
		} else {
			lo_bitoff += 2;
		}
		if(hi_bitoff == 0) {
			if(hi_word == 0) {
				break;
			}
			hi_word--;
			hi_bitoff = 62;
		} else {
			hi_bitoff -= 2;
		}
	}
}

/**
 * Return the pointer for whichever k-mer is minimal.
 */
static b256 *b256_min(b256 *a, b256 *b) {
	for(int i = 3; i >= 0; i--) {
		if(a->u64s[i] < b->u64s[i]) {
			return a;
		} else if(a->u64s[i] != b->u64s[i]) {
			return b;
		}
	}
	return a;
}

/**
 * Given a string, extract every k-mer, canonicalize, and add to the QF.
 */
int cqf_string_injest(
    const char *read, // string whose k-mers to add
    size_t read_len,  // length of string
    struct flush_object *obj)
{
	const uint32_t k = obj->ksize;
	const uint32_t seed = obj->main_qf->metadata->seed;
	const __uint128_t range = obj->main_qf->metadata->range;
	assert(obj->local_qf == NULL || obj->local_qf->metadata->range == range);
	assert(obj->local_qf == NULL || obj->local_qf->metadata->seed == seed);
	if(k > 127) {
		fprintf(stderr, "No support for k-mers sizes over 127; k-mer size "
		                "%u was specified\n", k);
		return -1;
	}
	int nadded = 0;
	int i = 0;
	for(; i < read_len; i++) {
		b256 first, first_rev, *item = NULL;
		b256_clear(&first);
		b256_clear(&first_rev);
		// Phase 1: get initial k-mer
		int do_continue = 0;
		for(int start = i; i < start + k; i++) { // First kmer
			uint8_t curr = map_base(read[i]);
			if(NON_ACGT(curr)) {
				do_continue = 1;
				break;
			}
			b256_or_low(&first, curr);
			b256_lshift(&first);
		}
		if(do_continue) {
			continue;
		}
		b256_rshift_2(&first);
		b256_revcomp(&first, &first_rev, k);
		item = b256_min(&first, &first_rev);
		uint64_t hash = MurmurHash64A((void*)item, 32, seed) % range;
		qf_insert(obj->main_qf, hash, 0, 1, false, false);
		nadded++;
		// Phase 2: get all subsequent k-mers
		b256 next = first, next_rev = first_rev;
		for(; i < read_len; i++) {
			uint8_t curr = map_base(read[i]);
			if(NON_ACGT(curr)) {
				break;
			}
			b256_lshift_and_mask(&next, k);
			b256_or_low(&next, curr);
			uint64_t tmp = reverse_complement_base(curr);
			b256_rshift_2(&next_rev);
			b256_or_at(&next_rev, tmp, k - 1);
			item = b256_min(&next, &next_rev);
			hash = MurmurHash64A((void*)item, 32, seed) % range;
			qf_insert(obj->main_qf, hash, 0, 1, false, false);
			nadded++;
		}
	}
	return nadded;
}

/**
 * Assumes no locking is needed, i.e. that no other thread
 * is trying to update the CQF.
 */
int cqf_string_query(
    const char *read_orig,
    size_t read_len_orig,
    int64_t *count_array,
    size_t count_array_len,
    struct flush_object *obj)
{
	const uint32_t k = obj->ksize;
	const uint32_t seed = obj->main_qf->metadata->seed;
	const __uint128_t range = obj->main_qf->metadata->range;
	assert(obj->local_qf == NULL || obj->local_qf->metadata->range == range);
	assert(obj->local_qf == NULL || obj->local_qf->metadata->seed == seed);
	int64_t *cur_count = count_array;
	const char *read = read_orig;
	size_t read_len = read_len_orig;
	do {
		if(read_len < k) {
			return (int)(cur_count - count_array);
		}
		b256 first, first_rev, *item = NULL;
		b256_clear(&first);
		b256_clear(&first_rev);
		int do_continue = 0;
		for(int i = 0; i < k; i++) {
			uint8_t curr = map_base(read[i]);
			if(NON_ACGT(curr)) {
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
			b256_or_low(&first, curr);
			b256_lshift(&first);
		}
		if(do_continue) {
			continue;
		}

		b256_rshift_2(&first);
		b256_revcomp(&first, &first_rev, k);
		item = b256_min(&first, &first_rev);
		uint64_t hash = MurmurHash64A((void*)item, 32, seed);
		hash %= range;
		if(hash < FP_BUF_ELTS) {
			obj->fp_buf[hash]++;
		}
		uint64_t count = qf_count_key_value(obj->main_qf, hash, 0);
		*cur_count++ = count;
		if(cur_count - count_array > count_array_len) {
			fprintf(stderr, "Wrote off end of count array\n");
			exit(1);
		}

		b256 next = first, next_rev = first_rev;
		read += k;
		if(read_len < k) {
			fprintf(stderr, "Query string became too short\n");
			exit(1);
		}
		read_len -= k;

		for(uint32_t i = 0; i < read_len; i++) { //next kmers
			uint8_t curr = map_base(read[i]);
			if(NON_ACGT(curr)) {
				for(int j = 0; j < obj->ksize && (cur_count - count_array) < count_array_len; j++) {
					*cur_count++ = -1;
				}
				read++;
				read_len--;
				do_continue = 1;
				break;
			}
			b256_lshift_and_mask(&next, k);
			b256_or_low(&next, curr);
			uint64_t tmp = reverse_complement_base(curr);
			b256_rshift_2(&next_rev);
			b256_or_at(&next_rev, tmp, k - 1);
			item = b256_min(&next, &next_rev);
			hash = MurmurHash64A((void*)item, 32, seed);
			hash %= range;
			if(hash < FP_BUF_ELTS) {
				obj->fp_buf[hash]++;
			}
			count = qf_count_key_value(obj->main_qf, hash, 0);
			*cur_count++ = count;
			if(cur_count - count_array > count_array_len) {
				fprintf(stderr, "Wrote off end of count array\n");
				exit(1);
			}
		}
		if(!do_continue) {
			break;
		}
	} while(1);
	return (int)(cur_count - count_array);
}

/**
 * Put false positive rate information into the given count_array.  Returns up
 * to min(FP_BUF_ELTS, range) pairs of uint64_ts.  Each pair consists of an
 * observed count and a true count.
 */
int cqf_est_fpr(
    int64_t *count_array,
    size_t count_array_len,
    struct flush_object *obj)
{
	const __uint128_t range = obj->main_qf->metadata->range;
	assert(obj->local_qf == NULL || obj->local_qf->metadata->range == range);
	int LIMIT = FP_BUF_ELTS < range ? FP_BUF_ELTS : range;
	LIMIT = LIMIT < ((int)count_array_len/2) ? LIMIT : ((int)count_array_len/2);
	for(size_t i = 0; i < LIMIT; i++) {
		uint64_t obs_count = qf_count_key_value(obj->main_qf, i, 0);
		uint64_t true_count = obj->fp_buf[i];
		assert(obs_count >= true_count);
		*count_array++ = obs_count;
		*count_array++ = true_count;
	}
	return LIMIT;
}

#ifdef SQUEAKR_TEST_MAIN
static void init(QF *cf, struct flush_object *obj, int ksize, int qbits) {
	obj->local_qf = NULL;
	obj->main_qf = cf;
	obj->ksize = ksize;
	obj->count = 0;
	memset(obj->fp_buf, 0, sizeof(uint32_t)*FP_BUF_ELTS);
	
	int num_hash_bits = qbits+8;	// we use 8 bits for remainders in the main QF
	uint32_t seed = 2038074761;
	qf_init(cf, (1ULL << qbits), num_hash_bits, 0, true, "", seed);
}

static void quick_tests(unsigned seed) {
	{
		int64_t results[4];
		const char *text = "ACGTACG";
		//                  0123
		for(int i = 0; i < 3; i++) {
			QF cf; struct flush_object obj; init(&cf, &obj, 4, 10);
			cqf_string_injest(text, i + 4, &obj);
			cqf_string_query(text, i+4, results, i+1, &obj);
			for(int j = 0; j <= i; j++) {
				assert(results[j] == 1);
			}
		}
		QF cf; struct flush_object obj; init(&cf, &obj, 4, 10);
		cqf_string_injest(text, 7, &obj);
		cqf_string_query(text, 7, results, 4, &obj);
		assert(results[0] == 1);
		assert(results[1] == 2);
		assert(results[2] == 1);
		assert(results[3] == 2);
	}

    {
        QF cf; struct flush_object obj; init(&cf, &obj, 4, 10);
        const char *text  = "TCCCGGGAGGGA";
        const char *query = "TCCCNGGGA";
        int64_t results[6];
        cqf_string_injest(text, 12, &obj);
        cqf_string_query(query, 9, results, 6, &obj);
        assert(results[0] == 3);
        assert(results[5] == 3);
    }

	{
		QF cf; struct flush_object obj; init(&cf, &obj, 60, 10);
		srand(seed);
		int textlen = 100000;
		int ksize = 60;
		char *text = (char*)malloc(textlen + ksize);
		char *cur = text;
		for(int i = 0; i < textlen + ksize - 1; i++) {
			*cur++ = "ACGT"[rand() % 4];
		}
		text[textlen+ksize-1] = '\0';
		cqf_string_injest(text, textlen + ksize - 1, &obj);
	}
}

int main(int argc, char *argv[]) {

	if(argc < 5) {
		fprintf(stderr, "Need at least 4 args\n");
		return 1;
	}

	int ksize = atoi(argv[1]), qbits = atoi(argv[2]);
	fprintf(stderr, "Using ksize: %d\n", ksize);
	fprintf(stderr, "Using qbits: %d\n", qbits);
	unsigned seed = 777;

	quick_tests(seed);

	QF cf;
	struct flush_object obj;
	init(&cf, &obj, ksize, qbits);

	fprintf(stderr, "Building reference from: %s\n", argv[3]);
	cqf_string_injest(argv[3], strlen(argv[3]), &obj);

	for(size_t i = 4; i < argc; i++) {
		size_t len = strlen(argv[i]);
		size_t results_len = len - ksize + 1;
		int64_t *results = (int64_t*)malloc(sizeof(int64_t) * results_len);
		if(results == NULL) {
			fprintf(stderr, "ERROR: could not allocate results\n");
			return 1;
		}
		int nres = cqf_string_query(argv[i], len, results, results_len, &obj);
		assert(nres == (int)results_len);
		for(int j = 0; j < nres; j++) {
			fprintf(stderr, "[%d]: %lld\n", j, results[j]);
		}
		free(results);
	}
}
#endif
