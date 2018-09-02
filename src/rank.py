#!/usr/bin/env python

'''
Rank a set of variants for inclusion in a graph genome, from highest to lowest priority
'''

import sys
import argparse
import iohelp
from util import *

VERSION = '0.0.1'

class VarRanker:
    def __init__(self, genome, variants, r, phasing, max_v):
        self.genome = genome
        self.chrom_lens = dict()
        for chrom, seq in genome.items():
            self.chrom_lens[chrom] = len(seq)

        self.variants = variants
        self.num_v = len(variants)
        self.r = r

        self.phasing = phasing
        self.hap_parser = iohelp.HaplotypeParser(phasing) if phasing else None

        self.max_v_in_window = max_v

        self.h_ref = None
        self.h_added = None

        self.wgt_ref = None
        self.wgt_added = None

        self.curr_vars = None
        self.freqs = {}

    def avg_read_prob(self):
        #self.wgt_ref = 0.778096
        #self.wgt_added = 0.002113

        if self.wgt_ref and self.wgt_added:
            return

        variants = self.variants

        # Average probability (weighted by genome length) of a specific read from the linear genome being chosen 
        total_prob_ref = 0
        count_ref = 0

        # Average probability (weighted by genome length) of a specific read from the added pseudocontigs being chosen 
        total_prob_added = 0
        count_added = 0

        if self.hap_parser:
            self.hap_parser.reset_chunk()

        num_v = len(variants)
        r = self.r

        var_i = 0
        last_i, last_j, last_pref, last_added = -1, -1, -1, -1
        for chrom, seq in self.genome.items():
            print('Processing chrom %s' % chrom)
            num_reads = self.chrom_lens[chrom] - r + 1
            count_ref += num_reads
            for i in range(num_reads):
                read = seq[i:i+r]
                if 'N' in read or 'M' in read or 'R' in read:
                    count_ref -= 1
                    continue

                #total_prob_ref += 1

                # Set [var_i, var_j) to the range of variants contained in the current read
                while var_i < num_v and variants[var_i].chrom == chrom and variants[var_i].pos < i:
                    var_i += 1
                var_j = var_i
                while var_j < num_v and variants[var_i].chrom == chrom and variants[var_j].pos < i+r:
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
                        last_added *= (variants[c].num_alts + 1)
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
        print('Avg probability of reads in ref:  %f' % self.wgt_ref)
        print('Avg probability of added reads:   %f' % self.wgt_added)

    def count_kmers_ref(self):
        import dna_jellyfish
        if self.h_ref:
            return

        # Create new Jellyfish counter and count all kmers in reference genome
        dna_jellyfish.MerDNA.k(self.r)
        self.h_ref = dna_jellyfish.HashCounter(1024, 5)

        for chrom in self.genome.values():
            assert chrom is not None
            mers = dna_jellyfish.string_canonicals(chrom)
            assert mers is not None
            for m in mers:
                self.h_ref.add(m, 1)

    def count_kmers_added(self):
        import dna_jellyfish

        if self.h_added:
            return

        total = 0
        total_r = 0

        dna_jellyfish.MerDNA.k(self.r)
        self.h_added = dna_jellyfish.HashCounter(1024, 5)

        for i in range(self.num_v):
            chrom = self.variants[i].chrom
            pos = self.variants[i].pos

            # Number of variants in window starting at this one
            k = 1
            while i+k < self.num_v and self.variants[i+k].chrom == chrom and self.variants[i+k].pos < pos+self.r:
                k += 1

            if k > self.max_v_in_window:
                alt_freqs = [(sum(self.variants[i+j].probs), i+j) for j in range(1, k)]
                ids = [f[1] for f in sorted(alt_freqs, reverse=True)[:self.max_v_in_window-1]]
                it = PseudocontigIterator(self.genome[chrom], [self.variants[i]]+[self.variants[v] for v in ids], self.r)
            else:
                it = PseudocontigIterator(self.genome[chrom], self.variants[i:i+k], self.r)

            pseudocontig = it.next()
            while pseudocontig:
                # Add to jellyfish
                mers = dna_jellyfish.string_canonicals(pseudocontig)
                for m in mers:
                    self.h_added.add(m, 1)

                #total += 1
                #total_r += len(pseudocontig) - self.r + 1

                pseudocontig = it.next()

        #print('%d total pseudocontigs' % total)
        #print('%d total reads' % total_r)

    def prob_read(self, variants, var_ids, vec):
        '''
            Probability that a read contains the allele vector vec
            Simultaneously computes probabilities for all haplotypes of the given variants to save time
        '''

        if not self.curr_vars or not (self.curr_vars == var_ids):
            assert len(var_ids) > 0
            self.curr_vars = var_ids
            counts = [variants[v].num_alts for v in var_ids]
            num_v = len(var_ids)
            self.counts = counts
            if self.hap_parser:
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
                            p *= variants[var_ids[i]].probs[vtmp[i]-1]
                        else:
                            p *= (1 - sum(variants[var_ids[i]].probs))
                    self.freqs[vec_to_id(vtmp, counts)] = p
                    vtmp = get_next_vector(num_v, counts, vtmp)
                    done = vtmp is None

        f = self.freqs[vec_to_id(vec, self.counts)]
        return f

    def prob_read_ref(self, variants, var_ids):
        '''
            Probability that a read is from the reference genome, plus-one smoothed
            Faster than prob_read() when other haplotype probs are unneeded
        '''

        if self.hap_parser:
            self.curr_vars = var_ids
            self.counts = [variants[v].num_alts for v in var_ids]
            return self.hap_parser.get_ref_freq(var_ids, self.counts)
        else:
            p = 1
            for i in range(len(var_ids)):
                p *= (1 - sum(variants[var_ids[i]].probs))
            return p

    def seen_pcs(self, out_prefix):
        pcs = []
        with open('temp.txt', 'r') as f:
            for line in f:
                row = [int(a) for a in line.rstrip().split('\t')]
                pcs.append((row[0], row[1:]))
        iohelp.write_pcs(self.variants, pcs, out_prefix)
        exit()

        if not self.phasing:
            print('Phasing file is required to rank pseudocontigs')

        pcs = []
        with open(self.phasing, 'r') as f:
            numH = len(f.readline().rstrip().split(','))
            f.seek(0)

            block_firsts = [-1] * numH
            # Last alternate allele in block
            block_lasts = [-1] * numH
            block_vecs = [None] * numH

            v = 0
            for line in f:
                if (v+1) % 100000 == 0:
                    print('Processing line %d of %d' % (v+1, self.num_v))

                curr_pcs = []
                row = [int(a) for a in line.rstrip().split(',')]
                for i in range(numH):
                    if block_vecs[i]:
                        if (self.variants[v].pos - self.variants[block_lasts[i]].pos) < self.r:
                            if row[i]:
                                block_lasts[i] = v
                            block_vecs[i].append(row[i])
                        else:
                            pc = (block_firsts[i], block_vecs[i])
                            if not pc in curr_pcs:
                                curr_pcs.append(pc)
                            if row[i]:
                                block_firsts[i] = v
                                block_lasts[i] = v
                                block_vecs[i] = [row[i]]    
                            else:
                                block_firsts[i] = -1
                                block_lasts[i] = -1
                                block_vecs[i] = None
                    elif row[i]:
                        block_firsts[i] = v
                        block_lasts[i] = v
                        block_vecs[i] = [row[i]]
                pcs += curr_pcs
                v += 1
            curr_pcs = []
            for i in range(numH):
                if block_vecs[i]:
                    pc = (block_firsts[i], block_vecs[i])
                    if not pc in curr_pcs:
                        curr_pcs.append(pc)
            pcs += curr_pcs

        print('Sorting and writing')
        pcs = sorted(list(pcs))
        with open('temp.txt', 'w') as f:
            for pc in pcs:
                f.write(str(pc[0]) + '\t' + '\t'.join([str(a) for a in pc[1]]) + '\n')

        iohelp.write_pcs(self.variants, pcs, out_prefix)

    def rank_pcs(self, out_prefix, pcts):
        if not self.hap_parser:
            print('Phasing file is required to rank pseudocontigs')

        pc_wgts = []
        
        print('Computing weights...')
        #for i in range(self.num_v):
        for i in range(10):
            print(i)
            chrom = self.variants[i].chrom
            pos = self.variants[i].pos

            # Number of variants in window starting at this one
            k = 1
            while i+k < self.num_v and self.variants[i+k].chrom == chrom and self.variants[i+k].pos < pos+self.r:
                k += 1

            var_ids = range(i, i+k)
            counts = [self.variants[v].num_alts for v in var_ids]

            # Get allele vectors present in population, where first allele is always alt
            pcs = self.hap_parser.get_seen_pcs(var_ids, counts)

            start = str(pos+1)
            for pc in pcs:
                vec = pc[0]
                if not vec[0]:
                    continue

                last_id = k-1
                while not vec[last_id]:
                    last_id -= 1
                nreads = self.r - (self.variants[i+last_id].pos - self.variants[i].pos)

                freq = pc[1]
                wgt = freq * nreads
                pc_wgts.append((wgt, i, vec))

        print('Finished calculating weights, sorting...')
        pc_wgts.sort(reverse=True)

        print('Writing haplotypes...')
        used_vars = []
        for v in self.variants:
            used_vars.append([0] * v.num_alts)

        last_n = 0
        total_pcs = len(pc_wgts)
        for pct in pcts:
            print('Writing %d%% of pseudocontigs' % pct)

            curr_n = int((total_pcs * pct) / 100.0)
            for pc_id in range(last_n, curr_n):
                i = pc_wgts[pc_id][1]
                vec = pc_wgts[pc_id][2]
                
                for k in range(len(vec)):
                    if vec[k] > 0:
                        used_vars[i+k][vec[k]-1] = 1

            iohelp.write_pcs(self.variants, used_vars, sorted([(p[1],p[2]) for p in pc_wgts[:curr_n]]), out_prefix)

            last_n = curr_n

    def rank(self, method, out_file):
        ordered = None
        ordered_blowup = None
        print(method)
        if method == 'popcov':
            ordered = self.rank_pop_cov()
        elif method == 'popcov-blowup':
            ordered = self.rank_pop_cov(True)
        elif method == 'hybrid':
            ordered, ordered_blowup = self.rank_hybrid()

        if ordered:
            with open(out_file, 'w') as f:
                f.write('\t'.join([self.variants[i].chrom + ',' + str(self.variants[i].pos+1) for i in ordered]))
        if ordered_blowup:
            with open(out_file+'.blowup', 'w') as f:
                f.write('\t'.join([self.variants[i].chrom + ',' + str(self.variants[i].pos+1) for i in ordered_blowup]))

    def rank_hybrid(self, threshold=0.5):

        print('Counting kmers in ref')
        #time1 = time.time()
        self.count_kmers_ref()
        #time2 = time.time()
        #print('Ref counting time: %f m' % ((time2-time1)/60))
        print('Counting added kmers')
        self.count_kmers_added()
        #time3 = time.time()
        #print('Added counting time: %f m' % ((time3-time2)/60))

        print('Finished counting kmers')
        print('')

        self.avg_read_prob()
        #time4 = time.time()
        #print('Avg prob time: %f m' % ((time4-time3)/60))

        print('Computing hybrid scores')
        if self.hap_parser:
            self.hap_parser.reset_chunk()

        var_wgts = [0] * self.num_v

        for v in range(self.num_v):
            if v % 100000 == 0:
                print('Processing %d / %d variants' % (v, self.num_v))
            self.compute_hybrid(v, var_wgts)

        # Uncomment the following 2 lines to write SNP hybrid scores to an intermediate file
        #with open('hybrid_wgts.txt', 'w') as f_hybrid:
        #    f_hybrid.write(','.join([str(w) for w in var_wgts]))
        
        # Uncomment the following 2 lines and comment out all of the function above this line to resume computation from an intermediate hybrid scores file
        #with open('hybrid_wgts.txt', 'r') as f_hybrid:
        #    var_wgts = [float(w) for w in f_hybrid.readline().rstrip().split(',')]

        var_scores = [(var_wgts[i], i) for i in range(self.num_v)]
        var_scores.sort()
        ordered = [v[1] for v in var_scores]

        # Compute blowup ranking as well
        upper_tier = []
        lower_tier = []

        # Normalize weights to [0.01,1]
        var_wgts = [-w for w in var_wgts]
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

    def compute_hybrid(self, first_var, var_wgts):
        import dna_jellyfish

        r = self.r
        chrom = self.variants[first_var].chrom
        pos = self.variants[first_var].pos

        #if self.variants[first_var].pos < pos or self.variants[first_var].pos >= pos+r:
        #    return

        # Number of variants in window starting at this one
        k = 1
        while first_var+k < self.num_v and self.variants[first_var+k].chrom == chrom and self.variants[first_var+k].pos < pos+r:
            k += 1

        #if k > 14:
        #    sys.stdout.write('Processing variant %d with %d neighbors' % (first_var, k))

        if k > self.max_v_in_window:
            alt_freqs = [(sum(self.variants[first_var+j].probs), first_var+j) for j in range(1, k)]
            ids = [first_var] + [f[1] for f in sorted(alt_freqs, reverse=True)[:self.max_v_in_window-1]]
            it = PseudocontigIterator(self.genome[chrom], [self.variants[v] for v in ids], self.r)
        else:
            ids = range(first_var, first_var+k)
            it = PseudocontigIterator(self.genome[chrom], self.variants[first_var:first_var+k], r)

        pseudocontig = it.next()
        while pseudocontig:
            vec = it.curr_vec

            p = self.prob_read(self.variants, ids, vec)
            for i in range(len(pseudocontig) - self.r + 1):
                mer = dna_jellyfish.MerDNA(pseudocontig[i:i+r])
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
                hybrid_wgt = (p - avg_wgt) / (c_total)
                for j in range(len(ids)):
                    if vec[j]:
                        var_wgts[ids[j]] -= hybrid_wgt

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

        if not upper_tier and not lower_tier:
            return []

        threshold = penalty

        ordered = []
        tier_num = 0

        if not upper_tier:
            max_val = max(lower_tier)[0]
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

        while upper_tier:
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
    if args.output:
        o = args.output
    else:
        o = 'ordered.txt'
    if args.prune:
        max_v = args.prune
    else:
        max_v = r

    genome = iohelp.read_genome(args.reference, args.chrom)

    vars = iohelp.parse_1ksnp(args.vars)

    ranker = VarRanker(genome, vars, r, args.phasing, max_v)
    if args.pseudocontigs:
        ranker.seen_pcs(o)
    else:
        ranker.rank(args.method, o)


if __name__ == '__main__':

    if '--version' in sys.argv:
        print('ERG v' + VERSION)
        sys.exit(0)

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--method', type=str, required=True,
        help='Variant ranking method. Currently supported ranking methods: [popcov | popcov-blowup | hybrid]\n\'hybrid\' will produce hybrid ranking files both with and without blowup avoidance,')
    parser.add_argument('--reference', type=str, required=True, 
        help='Path to fasta file containing reference genome')
    parser.add_argument("--vars", type=str, required=True,
        help="Path to 1ksnp file containing variant information")
    parser.add_argument('--chrom', type=str,
        help="Name of chromosome from reference genome to process. If not present, process all chromosomes.")
    parser.add_argument('--window-size', type=int,
        help="Radius of window (i.e. max read length) to use. Larger values will take longer. Default: 35")
    parser.add_argument('--pseudocontigs', action="store_true", help='Rank pseudocontigs rather than SNPs')
    parser.add_argument('--phasing', type=str, required=False,
        help="Path to file containing phasing information for each individual")
    parser.add_argument('--output', type=str, required=False,
        help="Path to file to write output ranking to. Default: 'ordered.txt'")
    parser.add_argument('--prune', type=int, required=False,
        help='In each window, prune haplotypes by only processing up to this many variants. We recommend including this argument when ranking with the hybrid strategy for window sizes over 35.')

    args = parser.parse_args(sys.argv[1:])
    go(args)
