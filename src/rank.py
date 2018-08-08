#!/usr/bin/env python

"""
Rank a set of variants for inclusion in a graph genome, from highest to lowest priority
"""

from __future__ import print_function
import argparse
from iohelp import HaplotypeParser, read_genome, parse_1ksnp
import logging
import kmer_counter
import sys
from operator import itemgetter
from util import PseudocontigIterator, vec_to_id, get_next_vector
import re
import os


VERSION = '0.0.4'


class VarRanker:

    def __init__(self, genome, variants, r, phasing, max_v, counter_type, temp):
        logging.info('Creating ranker')
        self.genome = genome
        self.chrom_lens = dict()
        for chrom, seq in genome.items():
            self.chrom_lens[chrom] = len(seq)
        self.variants = variants
        self.num_v = sum(map(len, variants.values()))
        self.r = r
        self.phasing = phasing
        self.hap_parser = HaplotypeParser(phasing) if phasing else None
        self.max_v_in_window = max_v
        self.h_ref = self.h_added = None
        self.wgt_ref = self.wgt_added = None
        self.curr_vars = None
        self.freqs = {}
        self.counter_type = counter_type
        self.temp = temp

    def counter_maker(self, name, dont_care_below=1, cache_from=None, cache_to=None):
        if self.counter_type == 'Simple':
            return kmer_counter.SimpleKmerCounter(name, self.r,
                                                  dont_care_below=dont_care_below,
                                                  temp=self.temp,
                                                  cache_from=cache_from,
                                                  cache_to=cache_to)
        elif self.counter_type.startswith('KMC3'):
            toks = self.counter_type.split(',')
            threads, gb, batch_sz = 1, 4, -1
            if len(toks) > 1:
                threads = int(toks[1])
            if len(toks) > 2:
                gb = int(toks[2])
            if len(toks) > 3:
                batch_sz = int(toks[3])
            return kmer_counter.KMC3KmerCounter(name, self.r,
                                                threads=threads, gb=gb,
                                                batch_size=batch_sz,
                                                dont_care_below=dont_care_below,
                                                temp=self.temp,
                                                cache_from=cache_from,
                                                cache_to=cache_to)

    def avg_read_prob(self, cache_from=None, cache_to=None):
        logging.info('  Calculating average read probabilities')
        if self.wgt_ref and self.wgt_added:
            return

        if cache_from is not None:
            fn = os.path.join(cache_from, 'weights.csv')
            logging.info('    retrieving from cache file "%s"' % fn)
            with open(fn, 'rb') as fh:
                toks = fh.read().strip().split(',')
                self.wgt_ref = float(toks[0])
                self.wgt_added = float(toks[1])
        else:
            # Average probability (weighted by genome length) of a specific read from the linear genome being chosen 
            total_prob_ref, count_ref = 0, 0

            # Average probability (weighted by genome length) of a specific read from the added pseudocontigs being chosen 
            total_prob_added, count_added = 0, 0

            if self.hap_parser is not None:
                self.hap_parser.reset_chunk()

            r = self.r

            all_acgt_re = re.compile("^[ACGTacgt]*$")

            last_i, last_j, last_pref, last_added = -1, -1, -1, -1
            for chrom, seq in self.genome.items():
                variants = self.variants[chrom]
                num_v = len(variants)
                var_i = 0
                logging.info('    Processing chrom %s' % chrom)
                num_reads = self.chrom_lens[chrom] - r + 1
                count_ref += num_reads
                for i in range(num_reads):
                    read = seq[i:i+r]
                    if not all_acgt_re.match(read):
                        count_ref -= 1
                        continue

                    # Set [var_i, var_j) to the range of variants contained in the current read
                    while var_i < num_v and variants.poss[var_i] < i:
                        var_i += 1
                    var_j = var_i
                    while var_j < num_v and variants.poss[var_j] < i+r:
                        var_j += 1
                    num_vars = var_j - var_i

                    if num_vars == 0:
                        total_prob_ref += 1
                        continue

                    '''
                    counts = [variants[n].num_alts for n in range(var_i, var_j)]
                    p = 1 - self.prob_read(variants, range(var_i, var_j), [0]*num_vars)
                    total_prob_ref -= p
                    total_prob_added += p
                    num_pcs = 1
                    for c in range(var_i, var_j):
                        num_pcs *= (variants[c].num_alts + 1)
                    count_added += num_pcs-1

                    curr_count_added = 0
                    counts = [variants[n].num_alts for n in range(var_i, var_j)]
                    total_prob_ref2 += self.prob_read(variants, range(var_i, var_j), [0]*num_vars)

                    curr_p = self.prob_read(variants, range(var_i, var_j), [0]*num_vars)

                    vec = get_next_vector(num_vars, counts, [0]*num_vars)
                    while vec:
                        p = self.prob_read(variants, range(var_i, var_j), vec)
                        curr_p += p
                        total_prob_ref -= p
                        total_prob_added += p
                        count_added += 1
                        curr_count_added += 1

                        vec = get_next_vector(num_vars, counts, vec)

                    if abs(1 - curr_p) > 0.001:
                        print('%d (%d - %d)' % (num_vars, var_i, var_j))
                        print('Total prob: %f' % curr_p)
                    '''

                    if var_i == last_i and var_j == last_j:
                        p_ref = last_pref
                        p_alts = 1 - p_ref
                        count_added += last_added-1
                    else:
                        last_added = 1
                        for c in range(var_i, min(var_j, var_i+self.max_v_in_window)):
                            last_added *= (variants.num_alts(c) + 1)
                        count_added += last_added-1

                        p_ref = self.prob_read_ref(variants, range(var_i, var_j))
                        p_alts = 1 - p_ref
                    total_prob_ref += p_ref
                    total_prob_added += p_alts

                    last_i = var_i
                    last_j = var_j
                    last_pref = p_ref

            self.wgt_ref = float(total_prob_ref) / count_ref
            self.wgt_added = float(total_prob_added) / count_added

        if cache_to is not None:
            fn = os.path.join(cache_to, 'weights.csv')
            logging.info('  Caching average probabilities in "%s"' % fn)
            with open(fn, 'wb') as fh:
                fh.write(b'%f,%f\n' % (self.wgt_ref, self.wgt_added))

        logging.info('  Avg probability of reads in ref:  %f' % self.wgt_ref)
        logging.info('  Avg probability of added reads:   %f' % self.wgt_added)

    def count_kmers_ref(self, cache_from=None, cache_to=None):
        logging.info('  Counting reference k-mers')
        self.h_ref = self.counter_maker('Ref', dont_care_below=1,
                                        cache_from=cache_from,
                                        cache_to=cache_to)
        if cache_from is not None:
            logging.info('    retrieving from cache file "%s"' % self.h_ref.combined_db)
            return

        n_pcs, tot_pc_len = 0, 0
        for chrom in self.genome.values():
            self.h_ref.add(chrom)
            n_pcs += 1
            tot_pc_len += len(chrom)
            logging.info('    %d contigs, %d bases' % (n_pcs, tot_pc_len))
        self.h_ref.finalize()

    def count_kmers_added(self, cache_from=None, cache_to=None):
        logging.info('  Counting augmented k-mers')
        self.h_added = self.counter_maker('Aug', dont_care_below=2,
                                          cache_from=cache_from,
                                          cache_to=cache_to)

        if cache_from is not None:
            logging.info('    retrieving from cache file "%s"' % self.h_added.combined_db)
            self.h_added.finalize()
            return

        n_pcs, n_vars, max_vars = 0, 0, 0
        incr = 1.1
        bar = 1000

        for chrom, vars in self.variants.items():
            num_v = len(vars)
            chrom_seq = self.genome[chrom]
            for i in range(num_v):
                # Number of variants in window starting at this one
                k = 1
                while i+k < num_v and vars.poss[i+k] < vars.poss[i]+self.r:
                    k += 1

                # If necessary, select a subset of ALTs to include based on
                # allele frequency
                ids = range(i, i+k)
                if k > self.max_v_in_window:
                    alt_freqs = [(vars.sum_probs(i+j), i+j) for j in range(1, k)]
                    sortlist = list(sorted(alt_freqs, reverse=True))
                    ids = list(sorted([i] + [f[1] for f in sortlist[:self.max_v_in_window-1]]))

                n_vars += 1
                max_vars = max(max_vars, len(ids))

                for pc in PseudocontigIterator(chrom_seq, vars, ids, self.r):
                    n_pcs += 1
                    if n_pcs >= bar:
                        bar = max(bar+1, int(bar*incr))
                        logging.info('    %d contigs, %d vars; max vars=%d; %0.3f%% done' %
                                     (n_pcs, n_vars, max_vars, 100.0*float(n_vars)/self.num_v))
                    self.h_added.add(pc)
        self.h_added.finalize()

    def prob_read(self, variants, var_ids, vec):
        """
        Probability that a read contains the allele vector vec
        Simultaneously computes probabilities for all haplotypes of the given variants to save time
        """
        if not self.curr_vars or not (self.curr_vars == var_ids):
            self.curr_vars = var_ids
            num_v = len(var_ids)
            self.counts = counts = [variants.num_alts(v) for v in var_ids]
            if self.hap_parser is not None:
                # Initialize freqs based on population haplotype data
                self.freqs = self.hap_parser.get_freqs(var_ids, counts)
            else:
                # Inititalize freqs based on independently-assumed allele frequencies
                num_vecs = 1
                for c in counts:
                    num_vecs *= (c+1)
                vtmp = [0] * num_v
                done = False
                while not done:
                    p = 1
                    for i in range(num_v):
                        if vtmp[i]:
                            p *= variants.get_prob(var_ids[i], vtmp[i]-1)
                        else:
                            p *= (1 - variants.sum_probs(var_ids[i]))
                    self.freqs[vec_to_id(vtmp, counts)] = p
                    vtmp, done = get_next_vector(num_v, counts, vtmp)

        f = self.freqs[vec_to_id(vec, self.counts)]
        return f

    def prob_read_ref(self, variants, var_ids):
        """
        Probability that a read is from the reference genome, plus-one smoothed
        Faster than prob_read() when other haplotype probs are unneeded
        """
        if self.hap_parser is not None:
            self.curr_vars = var_ids
            self.counts = [variants[v].num_alts for v in var_ids]
            return self.hap_parser.get_ref_freq(var_ids, self.counts)
        else:
            p = 1
            for vi in var_ids:
                p *= (1 - variants.sum_probs(vi))
            return p

    def rank(self, method, out_file, query_batch_size=1000000,
             cache_from=None, cache_to=None):
        ordered_blowup = None
        if method == 'popcov':
            ordered = self.rank_pop_cov()
        elif method == 'popcov-blowup':
            ordered = self.rank_pop_cov(True)
        elif method == 'hybrid':
            self.count_kmers(cache_from=cache_from, cache_to=cache_to)
            ordered, ordered_blowup = self.rank_hybrid(query_batch_size=query_batch_size)
        else:
            raise RuntimeError('Bad ordering method: "%s"' % method)

        if ordered:
            with open(out_file, 'w') as f:
                f.write('\t'.join([tup[0] + ',' + str(tup[1] + 1) for tup in ordered]))

        if ordered_blowup:
            with open(out_file + '.blowup', 'w') as f:
                f.write('\t'.join([tup[0] + ',' + str(tup[1] + 1) for tup in ordered_blowup]))

    def count_kmers(self, cache_from=None, cache_to=None):
        logging.info('  Counting k-mers')
        self.count_kmers_ref(cache_from=cache_from, cache_to=cache_to)
        self.count_kmers_added(cache_from=cache_from, cache_to=cache_to)
        self.avg_read_prob(cache_from=cache_from, cache_to=cache_to)

    def rank_hybrid(self, threshold=0.5, query_batch_size=1000000):
        """
        `threshold` for blowup avoidance
        """
        logging.info('Hybrid-ranking variants')
        logging.info('  Computing hybrid scores')
        if self.hap_parser is not None:
            self.hap_parser.reset_chunk()

        var_wgts = []
        num_v = num_v_cum = 0
        bar = 100
        incr = 1.1

        def process_batch(batch, vecs, idss, variants, var_wgts_chrom):
            assert len(batch) == len(vecs)
            query_ref = self.h_ref.query_batch(batch)
            query_aug = self.h_added.query_batch(batch)
            for c_ref, c_aug, pc, vec, ids in zip(query_ref, query_aug, batch, vecs, idss):
                if c_ref == -1:
                    assert c_aug == -1
                    continue
                assert c_ref >= 0
                if c_aug == 0:
                    c_aug = 1
                p = self.prob_read(variants, ids, vec)
                if hist_ref is not None:
                    hist_ref.add(c_ref)
                if hist_aug is not None:
                    hist_aug.add(c_aug)
                c_total = c_ref + c_aug

                # Average relative probability of this read's other mappings
                avg_wgt = c_ref * self.wgt_ref + (c_aug - 1) * self.wgt_added
                hybrid_wgt = (p - avg_wgt) / (c_total)
                for j in range(len(ids)):
                    if vec[j]:
                        var_wgts_chrom[ids[j]] -= hybrid_wgt

        #hist_ref, hist_aug = Quantiler(), Quantiler()
        hist_ref, hist_aug = None, None
        for chrom, variants in self.variants.items():
            num_v = len(variants)
            var_wgts_chrom = [0] * num_v
            batch, vecs, idss = [], [], []
            for v in range(num_v):
                if num_v_cum + v >= bar:
                    bar = max(bar + 1, bar*incr)
                    msofar = (num_v_cum + v)
                    pct = msofar * 100.0 / self.num_v
                    logging.info('    %d out of %d variants; %0.3f%% done' % (msofar, self.num_v, pct))
                pos = variants.poss[v]
                r = self.r
                k = 1
                while v + k < num_v and variants.poss[v + k] < pos + r:
                    k += 1
                if k > self.max_v_in_window:
                    alt_freqs = [(variants.sum_probs(v + j), v + j) for j in range(1, k)]
                    ids = [v] + [f[1] for f in sorted(alt_freqs, reverse=True)[:self.max_v_in_window - 1]]
                    ids = list(sorted(ids))
                else:
                    ids = list(range(v, v + k))
                it = PseudocontigIterator(self.genome[chrom], variants, ids, r)
                for pc in it:
                    batch.append(pc)
                    vecs.append(it.vec[:])
                    idss.append(ids)
                    if len(batch) >= query_batch_size:
                        process_batch(batch, vecs, idss, variants, var_wgts_chrom)
                        batch, vecs, idss = [], [], []
            if len(batch) > 0:
                process_batch(batch, vecs, idss, variants, var_wgts_chrom)
            for v in range(num_v):
                var_wgts.append([var_wgts_chrom[v], chrom, variants.poss[v]])
            num_v_cum += num_v
            assert len(var_wgts) == num_v_cum

        #logging.info('Ref quantiles: ' + str(hist_ref))
        #logging.info('Aug quantiles: ' + str(hist_aug))

        assert num_v == self.num_v == len(var_wgts)

        # Uncomment the following 2 lines to write SNP hybrid scores to an
        # intermediate file

        #with open('hybrid_wgts.txt', 'w') as f_hybrid:
        #    f_hybrid.write(','.join([str(w) for w in var_wgts]))

        # Uncomment the following 2 lines and comment out all of the function
        # above this line to resume computation from an intermediate hybrid scores file

        #with open('hybrid_wgts.txt', 'r') as f_hybrid:
        #    var_wgts = [float(w) for w in f_hybrid.readline().rstrip().split(',')]

        logging.info('    Sorting')
        ordered = [(v[1], v[2]) for v in sorted(var_wgts)]
        # Tuples of (chrom, offset) ordered according to hybrid ranking,
        # before any blowup avoidance

        # Compute scores with blowup avoidance
        upper_tier, lower_tier = [], []

        # Normalize weights to [0.01,1]
        logging.info('    Normalizing weights')
        tmp_wgts = [-w for w in map(itemgetter(0), var_wgts)]
        assert max(tmp_wgts) > min(tmp_wgts)
        min_wgt = min(tmp_wgts)
        range_wgts = max(tmp_wgts) - min_wgt
        for i in range(self.num_v):
            new_wgt = (-var_wgts[i][0] - min_wgt)*0.99 / range_wgts + 0.01
            assert -0.01 <= new_wgt <= 1.01
            var_wgts[i][0] = new_wgt

        logging.info('    Blowup re-ranking')
        num_v = 0
        for chrom, variants in self.variants.items():
            for i in range(len(variants)):
                # var_wgts has to stay in order relative to self.variants
                # until at least here
                wgt = var_wgts[i+num_v][0]
                first = last = i
                while first > 0 and (variants.poss[i] - variants.poss[first-1]) < self.r:
                    first -= 1
                while last < (self.num_v-1) and (variants.poss[last+1] - variants.poss[i]) < self.r:
                    last += 1
                neighbors = last - first
                assert -0.01 <= wgt <= 1.01
                if wgt > threshold:
                    upper_tier.append((wgt, neighbors, chrom, i))
                else:
                    lower_tier.append((wgt, neighbors, chrom, i))
            num_v += len(variants)

        if len(upper_tier) == 0:
            logging.warning('  Upper tier empty')
        if len(lower_tier) == 0:
            logging.warning('  Lower tier empty')

        return ordered, self.rank_dynamic_blowup(upper_tier, lower_tier)

    def compute_hybrid(self, chrom_seq, variants, first_var, var_wgts, hist_ref, hist_aug):
        pos = variants.poss[first_var]
        num_v = len(variants)
        r = self.r

        #if self.variants[first_var].pos < pos or self.variants[first_var].pos >= pos+r:
        #    return

        # Number of variants in window starting at this one
        k = 1
        while first_var+k < num_v and variants.poss[first_var+k] < pos + r:
            k += 1

        if k > self.max_v_in_window:
            alt_freqs = [(variants.sum_probs(first_var+j), first_var+j) for j in range(1, k)]
            ids = [first_var] + [f[1] for f in sorted(alt_freqs, reverse=True)[:self.max_v_in_window-1]]
            ids = list(sorted(ids))
        else:
            ids = list(range(first_var, first_var+k))

        queries_per_batch = 10000
        batch, vecs = [], []
        it = PseudocontigIterator(chrom_seq, variants, ids, r)
        while True:
            for pc in it:
                batch.append(pc)
                vecs.append(it.vec)
                if len(batch) >= queries_per_batch:
                    break
            if len(batch) == 0:
                break
            assert len(batch) == len(vecs)
            query_ref = self.h_ref.query_batch(batch)
            query_aug = self.h_added.query_batch(batch)
            for c_ref, c_aug, pc, vec in zip(query_ref, query_aug, batch, vecs):
                if c_ref == -1:
                    assert c_aug == -1
                    continue
                assert c_ref >= 0
                assert c_aug >= 0
                p = self.prob_read(variants, ids, vec)
                if hist_ref is not None:
                    hist_ref.add(c_ref)
                if hist_aug is not None:
                    hist_aug.add(c_aug)
                if c_aug == 0:
                    print('Error! Read %s from added pseudocontigs could not be found (SNPs %d - %d)' %
                          (pc, first_var, first_var+k))
                    for j in range(first_var, first_var+k):
                        print('%s: %d, %s --> %s' % (self.variants[j].chrom, self.variants[j].pos,
                                                     self.variants[j].orig, ','.join(self.variants[j].alts)))
                    exit()
                c_total = c_ref + c_aug

                # Average relative probability of this read's other mappings
                avg_wgt = c_ref * self.wgt_ref + (c_aug-1) * self.wgt_added
                hybrid_wgt = (p - avg_wgt) / (c_total)
                for j in range(len(ids)):
                    if vec[j]:
                        var_wgts[ids[j]] -= hybrid_wgt
            batch, vecs = [], []

    def rank_pop_cov(self, with_blowup=False, threshold=0.5):
        if with_blowup:
            logging.info('Popcov+-ranking variants')
            upper_tier, lower_tier = [], []
            logging.info('  Creating initial tiers')
            for chrom, variants in self.variants.items():
                num_v = len(variants)
                for i in range(num_v):
                    wgt = variants.sum_probs(i)
                    first = last = i
                    while first > 0 and variants.poss[i] - variants.poss[first-1] < self.r:
                        first -= 1
                    while last < (num_v-1) and variants.poss[last+1] - variants.poss[i] < self.r:
                        last += 1
                    neighbors = last - first
                    if wgt > threshold:
                        upper_tier.append((wgt, neighbors, chrom, i))
                    else:
                        lower_tier.append((wgt, neighbors, chrom, i))
            if len(upper_tier) == 0:
                logging.warning('  Upper tier empty')
            if len(lower_tier) == 0:
                logging.warning('  Lower tier empty')
            logging.info('  Applying dynamic re-ranking')
            ordered = self.rank_dynamic_blowup(upper_tier, lower_tier)
        else:
            logging.info('Popcov-ranking variants')
            # Variant weight is the sum of frequencies of alternate alleles
            var_wgts = []
            for chrom, variants in self.variants.items():
                num_v = len(variants)
                var_wgts += [(-variants.sum_probs(i), chrom, variants.poss[i]) for i in range(num_v)]
            var_wgts.sort()
            ordered = [(v[1], v[2]) for v in var_wgts]

        return ordered

    def rank_dynamic_blowup(self, upper_tier, lower_tier, penalty=0.5):
        """
        Variants in tiers should be tuples, each of the form (weight, # neighbors, index in self.variants) 
        penalty: Weight multiplier for each variant every time a nearby variant is added to the graph
        """

        if not upper_tier and not lower_tier:
            return []

        threshold = penalty
        ordered = []

        if not upper_tier:
            max_val = max(lower_tier)[0]
            while max_val <= threshold:
                threshold *= penalty
            upper_tier, new_lower = [], []
            for i in range(len(lower_tier)):
                if lower_tier[i][0] > threshold:
                    upper_tier.append(lower_tier[i])
                else:
                    new_lower.append(lower_tier[i])
            lower_tier = new_lower[:]

        while upper_tier:
            # elements are [wgt, neighbors, chrom, i]
            upper_tier.sort(key=lambda x:(-x[0], x[1]))

            # Maps (chrom, i) coordinate in self.variants to id in upper/lower tier list
            vmap = {}
            for i, tup in enumerate(upper_tier):
                vmap[(tup[2], tup[3])] = (0, i)
            for i, tup in enumerate(lower_tier):
                vmap[(tup[2], tup[3])] = (1, i)

            for var_id, v in enumerate(upper_tier):
                if v[0] < 0:
                    continue
                vchrom, vid = v[2], v[3]
                variants = self.variants[vchrom]
                num_v = len(variants)
                pos = variants.poss[vid]
                ordered.append((v[2], pos))

                # Update other SNP weights
                first = v[3]
                last = first
                while first > 0 and (pos - variants.poss[first-1]) < self.r:
                    first -= 1
                while last < (num_v-1) and (variants.poss[last+1] - pos) < self.r:
                    last += 1

                if last > first:
                    for j in range(first, last+1):
                        if j == vid or (vchrom, j) not in vmap:
                            continue

                        tier, id = vmap[(vchrom, j)]
                        if tier == 0:  # upper
                            if id <= var_id:
                                continue
                            demote = upper_tier[id]
                            lower_tier.append((demote[0] * penalty, demote[1], demote[2], demote[3]))
                            vmap[(vchrom, j)] = (1, len(lower_tier)-1)
                            upper_tier[id] = (-1, demote[1], demote[2], demote[3])
                        else:
                            demote = lower_tier[id]
                            lower_tier[id] = (demote[0] * penalty, demote[1], demote[2], demote[3])

            if not lower_tier:
                break
            max_val = max(lower_tier)[0]
            if max_val > threshold:
                raise RuntimeError('Missed a point above threshold! max_val was '
                                   '%f and threshold was %f' % (max_val, threshold))
            while max_val <= threshold:
                threshold *= penalty
            upper_tier, new_lower = [], []
            for i in range(len(lower_tier)):
                if lower_tier[i][0] > threshold:
                    upper_tier.append(lower_tier[i])
                else:
                    new_lower.append(lower_tier[i])
            lower_tier = new_lower[:]

        return ordered


def go(args):
    r = args.window_size or 35
    o = args.output or 'ordered.txt'
    max_v = args.prune or r
    counter_type = args.counter or 'KMC3'

    logging.info('Reading genome')
    genome = read_genome(args.reference, target_chrom=args.chrom)

    logging.info('Parsing 1ksnp')
    vars = parse_1ksnp(args.vars, G=genome, target_chrom=args.chrom)

    ranker = VarRanker(genome, vars, r, args.phasing, max_v, counter_type, args.temp)
    if args.pseudocontigs:
        raise RuntimeError('--pseudocontigs not supported')
        #ranker.seen_pcs(o)
    else:
        if args.cache_to is not None:
            if not os.path.exists(args.cache_to):
                os.makedirs(args.cache_to)
        ranker.rank(args.method, o, query_batch_size=args.query_batch,
                    cache_from=args.cache_from, cache_to=args.cache_to)


if __name__ == '__main__':
    format_str = '%(asctime)s:%(levelname)s:%(message)s'
    level = logging.INFO
    logging.basicConfig(format=format_str, datefmt='%m/%d/%y-%H:%M:%S', level=level)

    if '--version' in sys.argv:
        print('ERG v' + VERSION)
        sys.exit(0)

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--method', type=str, required=True,
        help='Variant ranking method. Currently supported ranking methods: '
             '[popcov | popcov-blowup | hybrid]\n\'hybrid\' will produce '
             'hybrid ranking files both with and without blowup avoidance,')
    parser.add_argument('--reference', type=str, required=True, 
        help='Path to fasta file containing reference genome')
    parser.add_argument('--temp', type=str, required=False, 
        help='Set base temporary directory')
    parser.add_argument("--vars", type=str, required=True, metavar='<path>',
        help="1ksnp file containing variant information")
    parser.add_argument('--chrom', type=str, metavar='<name>',
        help='Name of chromosome from reference genome to process. If not '
             'present, process all chromosomes.')
    parser.add_argument('--window-size', type=int, metavar='<int>',
        help='Radius of window (i.e. max read length) to use. Larger values '
             'will take longer. Default: 35')
    parser.add_argument('--query-batch', type=int, default=100000000, metavar='<int>',
        help='Number of queries to batch when querying k-mer counts; '
             'default: 100000000')
    parser.add_argument('--pseudocontigs', action="store_true",
        help=argparse.SUPPRESS) # help='Rank pseudocontigs rather than SNPs')
    parser.add_argument('--phasing', type=str, required=False, metavar='<path>',
        help="File containing phasing information for each individual")
    parser.add_argument('--output', type=str, required=False, metavar='<path>',
        help="File to write output ranking to. Default: 'ordered.txt'")
    parser.add_argument('--cache-to', type=str, required=False,
        metavar='<path>', help="Cache results in given directory")
    parser.add_argument('--cache-from', type=str, required=False,
        metavar='<path>', help="Use cached results from given directory")
    parser.add_argument('--counter', type=str, required=False, default='KMC3',
        help='Type of counter to use; options are: "Simple"; '
             '"KMC3,<threads>,<gb>,<batch_size>"')
    parser.add_argument('--prune', type=int, required=False,
        help='In each window, prune haplotypes by only processing up to '
             'this many variants. We recommend including this argument when '
             'ranking with the hybrid strategy for window sizes over 35.')

    args = parser.parse_args(sys.argv[1:])
    go(args)
