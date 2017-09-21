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

        variants = self.variants

        # Average probability (weighted by genome length) of a specific read from the linear genome being chosen 
        total_prob_ref = sum(self.chrom_lens.values())
        count_ref = 0

        # Average probability (weighted by genome length) of a specific read from the added pseudocontigs being chosen 
        total_prob_added = 0
        count_added = 0

        if self.hap_parser:
            self.hap_parser.reset_chunk()

        num_v = len(variants)
        r = self.r

        var_i = 0
        amb = 0.0
        for chrom, seq in self.genome.items():
            print('Processing chrom %s' % chrom)
            num_reads = self.chrom_lens[chrom] - r + 1
            count_ref += num_reads
            for i in range(num_reads):
                read = seq[i:i+r]
                if 'N' in read or 'M' in read or 'R' in read:
                    continue

                # Set [var_i, var_j) to the range of variants contained in the current read
                while var_i < num_v and variants[var_i].chrom == chrom and variants[var_i].pos < i:
                    var_i += 1
                var_j = var_i
                while var_j < num_v and variants[var_i].chrom == chrom and variants[var_j].pos < i+r:
                    var_j += 1
                num_vars = var_j - var_i

                if num_vars == 0:
                    continue

                counts = [variants[id].num_alts for id in range(var_i,var_j)]
                vec = get_next_vector(num_vars, counts, [0]*num_vars)
                while vec:
                    #new_read = list(read)
                    #for v in range(num_vars):
                    #    if vec[v] > 0:
                    #        new_read = list() + variants[var_i+v].alts[vec[v]-1] + list()

                    p = self.prob_read(variants, range(var_i, var_j), vec)
                    total_prob_ref -= p
                    total_prob_added += p
                    count_added += 1

                    vec = get_next_vector(num_vars, counts, vec)

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

    def prob_read(self, variants, var_ids, vec, debug=False):
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
        ordered_blowup = None
        print(method)
        if method == 'popcov':
            ordered = self.rank_pop_cov()
        elif method == 'popcov-blowup':
            ordered = self.rank_pop_cov(True)
        elif method == 'amb':

            with open('debug.txt', 'w') as f_debug:
                f_debug.write('Ranking by amgiguity...\n')

            ordered, ordered_blowup = self.rank_ambiguity()

        with open('debug.txt', 'w') as f_debug:
            f_debug.write('Finished ranking, writing...\n')

        if ordered:
            with open(out_file, 'w') as f:
                f.write('\t'.join([self.variants[i].chrom + ',' + str(self.variants[i].pos+1) for i in ordered]))
        if ordered_blowup:
            with open(out_file+'.blowup', 'w') as f:
                f.write('\t'.join([self.variants[i].chrom + ',' + str(self.variants[i].pos+1) for i in ordered_blowup]))

    def rank_ambiguity(self, threshold=0.5):
        print('Counting kmers in ref')

        with open('debug.txt', 'w') as f_debug:
            f_debug.write('  Counting kmers in ref...\n')

        self.count_kmers_ref()
        print('Counting added kmers')
        self.count_kmers_added()

        with open('debug.txt', 'w') as f_debug:
            f_debug.write('  Counting added kmers...\n')

        print('Finished counting kmers')
        print('')

        with open('debug.txt', 'w') as f_debug:
            f_debug.write('  Computing avg read prob...\n')


        self.avg_read_prob()

        with open('debug.txt', 'w') as f_debug:
            f_debug.write('  Computing ambiguities...\n')

        print('Computing ambiguities')
        if self.hap_parser:
            self.hap_parser.reset_chunk()

        var_wgts = [0] * self.num_v
        #for chrom, seq in self.genome.items():
        #    for var_id in range(self.num_v):
        #        if self.variants[var_id].chrom == chrom:
        #            break

        #    for i in range(len(seq) - self.r + 1):
        #        while var_id < self.num_v and self.variants[var_id].chrom == chrom and self.variants[var_id].pos < i:
        #            var_id += 1
        #        if var_id == self.num_v or not self.variants[var_id].chrom == chrom:
        #            break
        #        elif self.variants[var_id].pos < i+self.r:
        #            self.compute_ambiguity(chrom, i, var_id, var_wgts)

        for v in range(self.num_v):
            self.compute_ambiguity(v, var_wgts)

        with open('amb_wgts.txt', 'w') as f_amb:
            f_amb.write(','.join([str(w) for w in var_wgts]))

        var_ambs = [(var_wgts[i], i) for i in range(self.num_v)]
        var_ambs.sort()
        ordered = [v[1] for v in var_ambs]

        # Compute blowup ranking as well
        upper_tier = []
        lower_tier = []

        # Normalize weights to [0.01,1]
        min_wgt = min(var_wgts)
        range_wgts = max(var_wgts) - min_wgt
        for i in range(self.num_v):
            var_wgts[i] = (var_wgts[i] - min_wgt)*0.99 / range_wgts + 0.01


        for i in range(self.num_v):
            wgt = var_wgts[i]

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
        ordered_blowup = self.rank_dynamic_blowup(upper_tier, lower_tier)

        return ordered, ordered_blowup

    def compute_ambiguity(self, first_var, var_ambs):
        r = self.r
        chrom = self.variants[first_var].chrom
        pos = self.variants[first_var].pos

        #if self.variants[first_var].pos < pos or self.variants[first_var].pos >= pos+r:
        #    return

        # Number of variants in window starting at this one
        k = 1
        while first_var+k < self.num_v and self.variants[first_var+k].chrom == chrom and self.variants[first_var+k].pos < pos+r:
            k += 1

        it = PseudocontigIterator(self.genome[chrom], self.variants[first_var:first_var+k], r)

        pseudocontig = it.next()
        while pseudocontig:
            vec = it.curr_vec

            p = self.prob_read(self.variants, range(first_var, first_var+k), vec, debug=True)
            for i in range(len(pseudocontig) - self.r + 1):
                mer = jellyfish.MerDNA(pseudocontig[i:i+r])
                mer.canonicalize()
                c_linear = self.h_ref[mer]
                if not c_linear:
                    c_linear = 0
                c_added = self.h_added[mer]
                if not c_added:
                    c_added = 0
                    if c_added == 0:
                        print('Error! Read %s from added pseudocontigs could not be found (SNPs %d - %d)' % (pseudocontig[i:i+r], first_var, first_var+k))
                        for j in range(first_var, first_var+k):
                            print('%s: %d, %s --> %s' % (self.variants[j].chrom, self.variants[j].pos, self.variants[j].orig, ','.join(self.variants[j].alts)))
                        exit()
                c_total = c_linear + c_added

                if c_total == 0:
                    print('Variants %d -%d / %d' % (first_var, first_var+k-1, self.num_v))
                    print('Vector:       ' + str(vec))
                    print('Pseudocontig: ' + str(pseudocontig))
                    print('Read:         ' + str(pseudocontig[i:i+r]))
                    exit()

                # Average relative probability of this read's other mappings
                avg_wgt = c_linear * self.wgt_ref + (c_added-1) * self.wgt_added
                amb_wgt = (p - avg_wgt) / (c_total)
                for j in range(k):
                    if vec[j]:
                        #amb_added += p - (p / float(c_total))
                        var_ambs[first_var+j] = amb_wgt

            pseudocontig = it.next()

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

    with open('debug.txt', 'w') as f_debug:
        f_debug.write('Reading genome...\n')

    genome = io.read_genome(args.reference, args.chrom)

    with open('debug.txt', 'w') as f_debug:
        f_debug.write('Parsing variants...\n')

    vars = io.parse_1ksnp(args.vars)

    with open('debug.txt', 'w') as f_debug:
        f_debug.write('Ranking variants...\n')

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
