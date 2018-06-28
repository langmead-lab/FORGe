#!/usr/bin/env python

import sys
import argparse
import bisect

# Update a genome with a set of SNPs and write to a new file


def write_snps(snps, locs, outfile):
    chrom_start = dict()
    curr_chr = None
    for i in range(len(locs)):
        l = locs[i]
        if not l[0] == curr_chr:
            if l[0] in chrom_start:
                print('Error! Reads are not sorted!')
            chrom_start[l[0]] = i
            curr_chr = l[0]

    f_out = open(outfile, 'w')

    last_loc = None
    unique_count = 0

    skipped_dels = 0
    skipped_ins = 0

    curr_chr = None
    curr_id = 0
    num_alts = 0
    count_added = 0
    num_target = len(locs)
    found_snps = [0] * num_target

    freqs = dict()
    min_freq = 1
    with open(snps, 'r') as f:
        num_lines = 0
        for line in f:
            num_lines += 1
            row = line.rstrip().split('\t')
            chrom = row[0]
            loc = int(row[1])

            if not chrom == curr_chr:
                curr_id = chrom_start[chrom]
                curr_chr = chrom

            while curr_id < num_target and loc > locs[curr_id][1]:
                curr_id += 1
                num_alts = 0

            if curr_id < num_target and loc == locs[curr_id][1]:
                found_snps[curr_id] = 1

                orig = row[2]
                alt = row[3]

                freq = float(row[4])
                if freq in freqs:
                    freqs[freq] += 1
                else:
                    freqs[freq] = 1
                if freq < min_freq:
                    min_freq = freq

                # Convert locations from 1-indexed to 0-indexed
                if len(orig) == 1 and len(alt) == 1:
                    # SNP
                    f_out.write(row[7] + '.' + str(num_alts) + '\tsingle\t' + row[0] + '\t' + str(loc-1) + '\t' + row[3] + '\n')
                elif len(orig) > len(alt):
                    # Deletion
                    if orig[:len(alt)] == alt:
                        f_out.write(row[7] + '.' + str(num_alts) + '\tdeletion\t' + row[0] + '\t' + str(loc-1+len(alt)) + '\t' + str(len(orig)-len(alt)) + '\n')
                    else:
                        skipped_dels += 1
                elif len(orig) < len(alt):
                    # Insertion
                    if orig == alt[:len(orig)]:
                        f_out.write(row[7] + '.' + str(num_alts) + '\tinsertion\t' + row[0] + '\t' + str(loc-1+len(orig)) + '\t' + alt[len(orig):] + '\n')
                    else:
                        skipped_ins += 1

                if not loc == last_loc:
                    unique_count += 1
                last_loc = loc
                num_alts += 1

            #elif chrom > locs[curr_id][0] or (chrom == locs[curr_id][0] and loc > locs[curr_id][1]):
            #    print('Couldn\'t find %s, %d!' % (locs[curr_id][0], locs[curr_id][1]))
            #    print('%s, %d' % (chrom, loc))
            #    exit()

    f_out.close()
    print('Skipped %d deletions, %d insertions' % (skipped_dels, skipped_ins))
    print('Found %d / %d target SNPs' % (unique_count, num_target))

    #print(min_freq)
    #print(', '.join([str(f)+': '+str(freqs[f]) for f in sorted(freqs.keys())]))

def read_sorted(sorted_snps, pct):
    with open(sorted_snps, 'r') as f:
        #locs = f.readline().rstrip().split('\t')
        #locs = [(int(l.split(',')[1]), int(l.split(',')[0])) for l in locs]
        row = f.readline().rstrip().split('\t')
        locs = [(r.split(',')[0], int(r.split(',')[1])) for r in row]
        print('%d total variants' % len(locs))
        num = int(len(locs) * pct / 100)
        return sorted(locs[:num])

if __name__ == '__main__':

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--snps", type=str, required=True, help="Path to 1ksnp file containing mutation information")
    parser.add_argument("--out", type=str, required=True, help="Path to output file")
    parser.add_argument("--sorted-snps", type=str, required=False, help="Path to file containing sorted SNP locations, output by SNP ranking program")
    parser.add_argument("--pct", type=float, required=False, help="Percentage of SNPs to include")

    args = parser.parse_args(sys.argv[1:])

    if args.pct <= 0:
        print('Not writing a file for 0\% of SNPs')
        exit()
    elif args.pct > 100:
        args.pct = 100


    if args.sorted_snps:
        locs = read_sorted(args.sorted_snps, args.pct)
    else:
        locs = None
    write_snps(args.snps, locs, args.out)

