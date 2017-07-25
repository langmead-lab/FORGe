#! /usr/bin/env python2.7

'''
Rank a set of variants for inclusion in a graph genome, from highest to lowest priority
'''

import sys
import argparse
import jellyfish
import io
import variant
from util import *

VERSION = '0.0.1'

class VarRanker:
    def __init__(self, genome, variants, r, phasing, debug=False):
        self.genome = genome
        self.chrom_lens = dict()
        for chrom, seq in genome.items():
            self.chrom_lens[chrom] = len(seq)

        self.variants = variants
        self.num_v = len(variants)
        self.r = r

        if phasing:
            self.hap_parser = io.HaplotypeParser(phasing)

        self.debug = debug

        self.h_ref = None
        self.h_added = None

        self.wgt_ref = None
        self.wgt_added = None

        self.curr_vars = None

    def avg_read_prob(self):
        if self.wgt_ref and self.wgt_added:
            return

        # Average probability (weighted by genome length) of a specific read from the linear genome being chosen 
        total_prob_ref = sum(self.chrom_lens.values())
        count_ref = 0

        # Average probability (weighted by genome length) of a specific read from the added pseudocontigs being chosen 
        total_prob_added = 0
        count_added = 0

        if self.hap_parser:
            self.hap_parser.reset_chunk()

        var_i = 0
        amb = 0.0
        for chrom, seq in self.genome.items():
            count_ref += len(seq) - self.r + 1
            for i in range(self.chrom_lens[chrom] - self.r + 1):
                read = seq[i:i+self.r]
                if 'N' in read or 'M' in read or 'R' in read:
                    continue

                # Set [var_i, var_j) to the range of variants contained in the current read
                while var_i < self.num_v and self.variants[var_i].chrom == chrom and self.variants[var_i].pos < i:
                    var_i += 1
                var_j = var_i
                while var_j < self.num_v and self.variants[var_i].chrom == chrom and self.variants[var_j].pos < i+self.r:
                    var_j += 1
                num_vars = var_j - var_i

                if num_vars == 0:
                    continue

                counts = [self.variants[id].num_alts for id in range(var_i,var_j)]
                vec = get_next_vector(num_vars, counts, [0]*num_vars)
                while vec:
                    #new_read = list(read)
                    #for v in range(num_vars):
                    #    if vec[v] > 0:
                    #        new_read = list() + self.variants[var_i+v].alts[vec[v]-1] + list()

                    p = self.prob_read(self.variants, range(var_i, var_j), vec)
                    total_prob_ref -= p
                    total_prob_added += p
                    count_added += 1

                    vec = get_next_vector(num_vars, counts, None)

        self.wgt_ref = float(total_prob_ref) / count_ref
        self.wgt_added = float(total_prob_added) / count_added
        print('Avg probability of reads in ref:  %f' % self.wgt_ref)
        print('Avg probability of added reads:   %f' % self.wgt_added)

    def count_kmers_ref(self):
        if self.h_ref:
            return

        # Create new Jellyfish counter and count all kmers in reference genome
        jellyfish.MerDNA.k(self.r)
        self.h_ref = jellyfish.HashCounter(1024, 5)

        for chrom in self.genome.values():
            mers = jellyfish.string_canonicals(chrom)
            for m in mers:
                self.h_ref.add(m, 1)

    def count_kmers_added(self):
        if self.h_added:
            return

        jellyfish.MerDNA.k(self.r)
        self.h_added = jellyfish.HashCounter(1024, 5)

        for i in range(self.num_v):
            chrom = self.variants[i].chrom
            pos = self.variants[i].pos

            # Number of variants in window starting at this one
            k = 1
            while i+k < self.num_v and self.variants[i+k].chrom == chrom and self.variants[i+k].pos < pos+self.r:
                k += 1

            it = PseudocontigIterator(self.genome[chrom], self.variants[i:i+k], self.r)

            pseudocontig = it.next()
            while pseudocontig:
                # Add to jellyfish
                mers = jellyfish.string_canonicals(pseudocontig)
                for m in mers:
                    self.h_added.add(m, 1)

                pseudocontig = it.next()

    def prob_read(self, variants, var_ids, vec):
        '''
            Probability that a read contains the allele vector vec
        '''

        if self.hap_parser:
            if not self.curr_vars or not (self.curr_vars == var_ids):
                self.curr_vars = var_ids
                self.counts = [variants[v].num_alts for v in var_ids]
                self.freqs = self.hap_parser.get_freqs(var_ids, self.counts)

        f = self.freqs[vec_to_id(vec, self.counts)]
        return f

    def rank(self, method, out_file):
        ordered = None
        print(method)
        if method == 'popcov':
            ordered = self.rank_pop_cov()
        elif method == 'popcov-blowup':
            ordered = self.rank_pop_cov(True)
        elif method == 'amb':
            ordered = self.rank_ambiguity()

        if ordered:
            with open(out_file, 'w') as f:
                f.write('\t'.join([self.variants[i].chrom + ',' + str(self.variants[i].pos+1) for i in ordered]))

    def rank_ambiguity(self):
        self.count_kmers_ref()
        self.count_kmers_added()

        self.avg_read_prob()

        var_wgts = [(self.ambiguity(i), i) for i in range(self.num_v)]
        var_wgts.sort()
        ordered = [v[1] for v in var_wgts]

        return ordered

    def rank_pop_cov(self, with_blowup=False, threshold=0.5):
        if with_blowup:
            upper_tier = []
            lower_tier = []

            for i in range(self.num_v):
                wgt = sum(self.variants[i].probs)

                first = i
                last = i
                while first > 0 and self.variants[first-1].chrom == self.variants[i].chrom and (self.variants[i].pos - self.variants[first-1].pos) < self.r:
                    first -= 1
                while last < (self.num_v-1) and self.variants[last+1].chrom == self.variants[i].chrom and (self.variants[last+1].pos - self.variants[i].pos) < self.r:
                    last += 1
                neighbors = last - first

                if wgt > threshold:
                    upper_tier.append((wgt, neighbors, i))
                else:
                    lower_tier.append((wgt, neighbors, i))
            ordered = self.rank_dynamic_blowup(upper_tier, lower_tier)
        else:
            # Variant weight is the sum of frequencies of alternate alleles
            var_wgts = [(-sum(self.variants[i].probs), i) for i in range(self.num_v)]
            var_wgts.sort()
            ordered = [v[1] for v in var_wgts]

        return ordered

    def rank_dynamic_blowup(self, upper_tier, lower_tier, penalty=0.5):
        '''
        Variants in tiers should be tuples, each of the form (weight, # neighbors, index in self.variants) 
        penalty: Weight multiplier for each variant every time a nearby variant is added to the graph
        '''

        threshold = penalty

        ordered = []
        tier_num = 0
        while upper_tier:
            if self.debug:
                tier_num += 1
                print('Processing tier %d (> %f), %d SNPs (%d remaining)' % (tier_num, threshold, len(upper_tier), len(lower_tier)))

            upper_tier.sort(key=lambda x:(-x[0], x[1]))

            # Maps id in self.variants to id in upper/lower tier list
            vmap = [0] * self.num_v
            for i in range(len(upper_tier)):
                vmap[upper_tier[i][2]] = (0,i)
            for i in range(len(lower_tier)):
                vmap[lower_tier[i][2]] = (1,i)

            for var_id in range(len(upper_tier)):
                v = upper_tier[var_id]
                var = self.variants[v[2]]
                if v[0] < 0:
                    continue

                chrom = var.chrom
                pos = var.pos
                #ordered.append(str(pos+1))
                ordered.append(v[2])

                # Update other SNP weights
                first = v[2]
                last = first
                while first > 0 and self.variants[first-1].chrom == chrom and (pos - self.variants[first-1].pos) < self.r:
                    first -= 1
                while last < (self.num_v-1) and self.variants[last+1].chrom == chrom and (self.variants[last+1].pos - pos) < self.r:
                    last += 1

                if (last > first):
                    for j in range(first, last+1):
                        if j == v[2] or not vmap[j]:
                            continue

                        if vmap[j][0] == 0:
                            id = vmap[j][1]
                            if id <= var_id:
                                continue
                            lower_tier.append((upper_tier[id][0] * penalty, upper_tier[id][1], upper_tier[id][2]))
                            vmap[j] = (1, len(lower_tier)-1)
                            upper_tier[id] = (-1, upper_tier[id][1], upper_tier[id][2])
                        else:
                            id = vmap[j][1]
                            lower_tier[id] = (lower_tier[id][0] * penalty, lower_tier[id][1], lower_tier[id][2])

            if not lower_tier:
                break
            max_val = max(lower_tier)[0]
            if max_val > threshold:
                print('Error! Missed a point above threshold!')
                exit()
            while max_val <= threshold:
                threshold *= penalty
            upper_tier = []
            new_lower = []
            for i in range(len(lower_tier)):
                if lower_tier[i][0] > threshold:
                    upper_tier.append(lower_tier[i])
                else:
                    new_lower.append(lower_tier[i])
            lower_tier = new_lower[:]
        
        return ordered

def go(args):
    if args.window_size:
        r = args.window_size
    else:
        r = 35

    genome = io.read_genome(args.reference, args.chrom)
    vars = io.parse_1ksnp(args.vars)

    ranker = VarRanker(genome, vars, r, args.phasing)
    ranker.rank(args.method, 'ordered.txt')


if __name__ == '__main__':

    if '--version' in sys.argv:
        print('ERG v' + VERSION)
        sys.exit(0)

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--method', type=str, required=True,
        help='Variant ranking method. Currently supported ranking methods: [popcov | amb | hybrid | blowup | popcov-blowup]')
    parser.add_argument('--reference', type=str, required=True, 
        help='Path to fasta file containing reference genome')
    parser.add_argument("--vars", type=str, required=True,
        help="Path to 1ksnp file containing variant information")
    parser.add_argument('--chrom', type=str,
        help="Name of chromosome from reference genome to process. If not present, process all chromosomes.")
    parser.add_argument('--window-size', type=int,
        help="Radius of window (i.e. max read length) to use. Larger values will take longer. Default: 35")
    parser.add_argument('--phasing', type=str, required=False,
        help="Path to file containing phasing information for each individual")

    args = parser.parse_args(sys.argv[1:])
    go(args)
