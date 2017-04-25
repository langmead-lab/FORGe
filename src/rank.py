#! /usr/bin/env python2.7

'''
Rank a set of variants for inclusion in a graph genome, from highest to lowest priority
'''

import sys
import argparse
import jellyfish
import read_iterator
import io
import variant
from util import *

VERSION = '0.0.1'

class VarRanker:
    def __init__(self, genome, vars, r, phasing, debug=False):
        self.genome = genome
        self.chrom_lens = dict()
        for chrom, seq in genome.items():
            self.chrom_lens[chrom] = len(seq)

        self.vars = vars
        self.num_v = len(vars)
        self.r = r

        if phasing:
            self.hap_parser = io.HaplotypeParser(phasing)

        self.debug = debug

        self.h_ref = None
        self.h_added = None

    def ambiguity(self, var_id):


    def avg_read_prob(self):
        # Average probability (weighted by genome length) of a specific read from the linear genome being chosen 
        total_prob_linear = sum([len(s) for s in self.chrom_lens.values()])

        # Average probability (weighted by genome length) of a specific read from the added pseudocontigs being chosen 
        total_prob_added = 0
        count_added = 0

        if self.hap_parser:
            self.hap_parser.reset_chunk()

        var_i = 0
        amb = 0.0
        for chrom, seq in self.genomes.items():
            for i in range(self.chrom_lens[chrom] - r + 1):
                read = seq[i:i+self.r]
                if 'N' in read or 'M' in read or 'R' in read:
                    continue

                # Set [var_i, var_j) to the range of variants contained in the current read
                while var_i < self.num_v and self.vars[var_i][0] < i:
                    var_i += 1
                var_j = var_i
                while var_j < self.num_v and self.vars[var_j][0] < i+self.r:
                    var_j += 1
                num_vars = var_j - var_i

                if num_vars == 0:
                    continue
                else:
                    counts = [self.vars[id].num_alts for id in range(var_i,var_j)]
                    vec = io.get_next_vector(num_vars, counts, None)
                    while vec:
                        new_read = list(read)
                        for v in range(num_vars):
                            if vec[v] > 0:
                                new_read = list() + self.vars[var_i+v].alts[vec[v]-1] + list()

    def count_kmers_ref(self):
        if self.h_ref:
            return

        # Create new Jellyfish counter and count all kmers in reference genome
        jellyfish.MerDNA(self.r)
        self.h_ref = jellyfish.HashCounter(1024, 5)

        for chrom in self.genome.values:
            mers = jellyfish.string_canonicals(chrom)
            for m in mers:
                self.h_ref.add(m, 1)

    def count_kmers_added(self):
        if self.h_added:
            return

        jellyfish.MerDNA(self.r)
        self.h_added = jellyfish.HashCounter(1024, 5)

        for i in range(self.num_v):
            chrom = self.vars[i].chrom
            pos = self.vars[i].pos

            # Number of variants in window starting at this one
            k = 1
            counts = [self.vars[i].num_alts]
            while i+k < self.num_v and self.vars[i+k].chrom == chrom and self.vars[i+k].pos < pos+r:
                counts.append(self.vars[i+k].num_alts)
                k += 1

            segment_start = max(self.vars[i].pos - self.r + 1, 0)
            segment_end = min(self.vars[i].pos + self.r, self.chrom_lens[chrom])
            segment = self.genome[chrom][segment_start:segment_end]
            chunks = [segment[:self.vars[i].pos - segment_start]]
            for j in range(1, k):
                chunks.append(segment[self.vars[i+j-1].pos + len(self.vars[i+j-1].orig) - segment_start : self.vars[i+j].pos - segment_start])
            chunks.append(segments[self.vars[i+j].pos + len(self.vars[i+j].orig) - segment_start :])

            # Allele vector
            v = [1] + [0] * (k-1)
            while v:
                # j is the greatest index in v s.t. v[j] > 0
                j = k-1
                while v[j] == 0:
                    j -= 1

                # Create new pseudocontig
                start = max(self.vars[i+j].pos - self.r + 1, 0)
                end = min(pos + self.r - 1, self.chrom_lens[chrom])
                E = list(self.genome[chrom][start:end+1])
                for j in range(k):
                    if v[j] > 0:
                        var = self.vars[i+j]
                        E = E[:var.pos] + var.alts[j-1] + E[var.pos + len(var.orig):]
                pseudocontig = ''.join(E)

                # Add to jellyfish
                mers = jellyfish.string_canonicals(pseudocontig)
                for m in mers:
                    h_added.add(m, 1)

                v = self.get_next_vector(k, counts, v)

    def prob_read(self, vars, var_ids, vec):
        '''
            Probability that a read contains the allele vector vec
        '''

        if self.hap_parser:
            if not self.curr_vars or not (self.curr_vars == var_ids):
                self.curr_vars = var_ids
                self.counts = [v.num_alts for v in vars]
                self.freqs = hap_parser.get_freqs(var_ids, self.counts)

        f = self.freqs[vec_to_id(vec, self.counts)]
        return f

    def rank(self, method, out_file):
        ordered = None
        if method == 'pop_cov':
            ordered = self.rank_pop_cov()
        elif method == 'pop_cov_blowup':
            ordered = self.rank_pop_cov(True)
        elif method == 'ambiguity':
            ordered = self.rank_ambiguity()

        if ordered:
            with open(out_file, 'w') as f:
                f.write(','.join([str(i) for i in ordered]))

    def rank_ambiguity(self):
        self.count_kmers_ref()
        self.count_kmers_added()

        var_wgts = [(self.ambiguity(i), i) for i in range(self.num_v)]
        var_wgts.sort()
        ordered = [v[1] for v in var_wgts]

        return ordered

    def rank_pop_cov(self, with_blowup=False):
        if with_blowup:
            upper_tier = []
            lower_tier = []
            for i in range(self.num_v):
                wgt = sum(self.vars[i].probs)
                first = i
                last = i
                while first > 0 and self.vars[first-1][4] == self.vars[i][4] and (self.vars[i][0] - self.vars[first-1][0]) < self.r:
                    first -= 1
                while last < (self.num_v-1) and self.vars[last+1][4] == self.vars[i][4] and (self.vars[last+1][0] - self.vars[i][0]) < self.r:
                    last += 1
                neighbors = last - first

                if wgt > threshold:
                    upper_tier.append((wgt, neighbors, i))
                else:
                    lower_tier.append((wgt, neighbors, i))
            ordered = self.rank_dynamic_blowup(upper_tier, lower_tier)
        else:
            # Variant weight is the sum of frequencies of alternate alleles
            var_wgts = [(sum(self.vars[i].probs), i) for i in range(self.num_v)]
            var_wgts.sort()
            ordered = [v[1] for v in var_wgts]

        return ordered

    def rank_dynamic_blowup(self, upper_tier, lower_tier, penalty=0.5):
        '''
        Variants in tiers should be tuples, each of the form (weight, # neighbors, index in self.vars) 
        penalty: Weight multiplier for each variant every time a nearby variant is added to the graph
        '''

        threshold = penalty

        ordered = []
        tier_num = 0
        while upper_tier:
            tier_num += 1
            if self.debug:
                print('Processing tier %d (> %f), %d SNPs (%d remaining)' % (tier_num, threshold, len(upper_tier), len(lower_tier)))

            upper_tier.sort(key=lambda x:(-x[0], x[1]))

            # Maps id in self.vars to id in upper/lower tier list
            vmap = [0] * self.num_v
            for i in range(len(upper_tier)):
                vmap[upper_tier[i][2]] = (0,i)
            for i in range(len(lower_tier)):
                vmap[lower_tier[i][2]] = (1,i)

            for var_id in range(len(upper_tier)):
                v = upper_tier[var_id]
                var = self.vars[v[2]]
                if v[0] < 0:
                    continue

                chrom = var.chrom
                pos = var.pos
                ordered.append(str(pos+1))

                # Update other SNP weights
                first = v[2]
                last = first
                while first > 0 and self.vars[first-1][4] == chrom and (pos - self.vars[first-1][0]) < r:
                    first -= 1
                while last < (self.num_v-1) and self.vars[last+1][4] == chrom and (self.vars[last+1][0] - pos) < r:
                    last += 1

                if (last > first):
                    for j in range(first, last+1):
                        if j == v[2] or not vmap[j]:
                            continue

                        if vmap[j][0] == 0:
                            id = vmap[j][1]
                            if id <= var_id:
                                continue
                            lower_tier.append((upper_tier[id][0] / 2, upper_tier[id][1], upper_tier[id][2]))
                            vmap[j] = (1, len(lower_tier)-1)
                            upper_tier[id] = (-1, upper_tier[id][1], upper_tier[id][2])
                        else:
                            id = vmap[j][1]
                            lower_tier[id] = (lower_tier[id][0]/2, lower_tier[id][1], lower_tier[id][2])

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
    vars = io.read_variants(args.vars)

    ranker = VarRanker(genome, vars, r, args.phasing)


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
        help="Name of chromosome from reference genome to process. If not present, process all chromosomes."
    parser.add_argument('--window-size', type=int,
        help="Radius of window (i.e. max read length) to use. Larger values will take longer. Default: 35"
    parser.add_argument('--phasing', type=str, required=False,
        help="Path to file containing phasing information for each individual")

    args = parser.parse_args(sys.argv[1:])
    go(args)
