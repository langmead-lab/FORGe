#! /usr/bin/env python2.7

__author__ = 'jacobpritt'

import sys
import argparse

'''
Parse a vcf file and optional filters and write the SNPs in 1ksnp format
'''

def read_genome(filename):
    G = dict()
    c = None

    with open(filename, 'r') as f:
        for line in f:
            # Skip header line
            if line[0] == '>':
                if c:
                    G[chrom] = c
                chrom = line.split(' ')[0][1:]
                chrom = chrom.rstrip()
                c = ''
                continue

            c += line.rstrip()
    if c:
        G[chrom] = c

    return G

def read_filters(filename):
    filters = []
    with open(filename, 'r') as f:
        for line in f:
            l = line.rstrip()
            if len(l) > 0:
                filters.append(l)
    return filters

def get_mutation_type(info):
    '''
        Return the value of the VT attribute from the info field
    '''

    attrs = info.split(';')
    for a in attrs:
        if a[:3] == 'VT=':
            return a[3:]

def process_vcf(filename):
    S = []
    with open(filename, 'r') as f:
        labels = None

        last = -1
        last_row = None

        for line in f:
            # Skip header lines
            if line[0] == '#' and line[1] == '#':
                continue

            if not labels:
                labels = line.rstrip().split('\t')

            else:
                row = line.rstrip().split('\t')
                type = get_mutation_type(row[7])
                if type == 'SNP':
                    chrom = row[0]
                    name = row[2]
                    loc = int(row[1])
                    alt = row[4]

                    if loc == last:
                        print('Duplicate rows:')
                        print(last_row)
                        print(line.rstrip())
                        print('')

                    if len(alt) > 1:
                        print('Multiallelic:')
                        print(line.rstrip())
                        print('')


                    last_row = line.rstrip()




def parse_vcf(filename, individuals=None, ingroup=None, outgroup=None, indiv_hap=None, indels=False):
    filters = ingroup or outgroup

    S = []

    print_filtered = False
    if print_filtered:
        f_out = open(filename+'.filtered', 'w')

    if individuals:
        f_ind = open(individuals, 'w')

    if indiv_hap:
        f_hapA = open(indiv_hap+'_hapA.1ksnp', 'w')
        f_hapB = open(indiv_hap+'_hapB.1ksnp', 'w')

    with open(filename, 'r') as f:
        labels = None

        line_id = 0

        for line in f:
            # Skip header lines
            if line[0] == '#' and line[1] == '#':
                if print_filtered:
                    f_out.write(line)
                continue

            if not labels:
                if print_filtered:
                    f_out.write(line)

                labels = line.rstrip().split('\t')

                indiv_col = None
                if indiv_hap:
                    for i in range(9, len(labels)):
                        if labels[i] == indiv_hap:
                            indiv_col = i
                    if not indiv_col:
                        print('Error! Couldn\'t find individual %s in VCF' % indiv)
                        exit()

                if ingroup:
                    cols = []
                    for i in range(9, len(labels)):
                        if labels[i] in ingroup:
                            cols.append(i)
                elif outgroup:
                    cols = []
                    for i in range(9, len(labels)):
                        if not labels[i] in outgroup:
                            cols.append(i)
            else:
                row = line.rstrip().split('\t')
                type = get_mutation_type(row[7])
                if type == 'SNP' or (indels and type == 'INDEL'):
                    chrom = row[0]
                    name = row[2]
                    loc = int(row[1])

                    if len(S) > 0 and S[-1][0] == loc:
                        # Rows with duplicate SNP locations are probably the result of some error in combining multiallelic SNPs
                        # We ignore the second row and treat it as a SNP with a single alt allele
                        continue

                    orig = row[3]
                    #if not len(orig) == 1:
                    #    print('Error! Row:')
                    #    print(row[:5])
                    alts = row[4].split(',')

                    # number of times the alternate allele appears
                    counts = [0] * len(alts)

                    if individuals:
                        alleles1 = []
                        #alleles2 = []

                    if not filters:
                        total = (len(row)-9) * 2
                        for i in range(len(row)-9):
                            c = row[9+i]
                            v1 = int(c[0])
                            if v1 > 0:
                                counts[v1-1] += 1
                            if len(c) > 1:
                              v2 = int(c[2])
                              if v2 > 0:
                                  counts[v2-1] += 1

                            if individuals:
                                alleles1.append(str(v1))
                                if len(c) > 1:
                                  alleles1.append(str(v2))
                                #alleles2.append(str(v2))
                    else:
                        total = len(cols) * 2
                        for i in range(len(cols)):
                            c = cols[i]
                            v1 = int(row[c][0])
                            if v1 > 0:
                                counts[v1-1] += 1
                            if len(row[c]) > 1:
                              v2 = int(row[c][2])
                              if v2 > 0:
                                counts[v2-1] += 1

                            if individuals:
                                alleles1.append(str(v1))
                                if len(row[c]) > 1:
                                  alleles1.append(str(v2))

                    if indiv_col:
                        alleleA = int(row[indiv_col][0])
                        if alleleA > 0:
                            new_row = [chrom, str(loc), orig, alts[alleleA-1], str(float(counts[alleleA-1])/total), '99', '1', name]
                            f_hapA.write('\t'.join(new_row) + '\n')

                        alleleB = int(row[indiv_col][2])
                        if alleleB > 0:
                            new_row = [chrom, str(loc), orig, alts[alleleB-1], str(float(counts[alleleB-1])/total), '99', '1', name]
                            f_hapB.write('\t'.join(new_row) + '\n')

                    if sum(counts)> 0:
                        if print_filtered:
                            f_out.write(line)

                        if len(S) == 0:
                            S = [(loc, alts, counts, total, chrom, name, orig)]
                        else:
                            i = len(S)
                            while S[i-1][0] > loc:
                                i -= 1
                            if i == len(S):
                                S.append((loc, alts, counts, total, chrom, name, orig))
                            else:
                                S = S[:i] + [(loc, alts, counts, total, chrom, name, orig)] + S[i:]

                        if individuals:
                            f_ind.write(','.join(alleles1) + '\n')
                            #f_ind.write(','.join(alleles2) + '\n')
                    line_id += 1
    if print_filtered:
        f_out.close()

    if individuals:
        f_ind.close()

    if indiv_hap:
        f_hapA.close()
        f_hapB.close()

    #print(S)
    for i in range(len(S)):
        counts = S[i][2]
        total = float(S[i][3])
        probs = [0] * len(counts)
        for j in range(len(counts)):
            probs[j] = float(counts[j]) / total
        S[i] = (S[i][0], S[i][1], probs, S[i][4], S[i][5], S[i][6])
    #print('')
    #print(S)

    return S

def write_1ksnp(S, G, filename):
    count_dup0 = 0
    count_dup1 = 0
    count_dup2 = 0
    with open(filename, 'w') as f:
        for snp in S:
            for i in range(len(snp[1])):
                if G[snp[3]][snp[0]] == snp[5]:
                    count_dup0 += 1
                if G[snp[3]][snp[0]-1] == snp[5]:
                    count_dup1 += 1
                if G[snp[3]][snp[0]+1] == snp[5]:
                    count_dup2 += 1

                # TODO: Find out what '99' column is
                #row = [snp[3], str(snp[0]), G[snp[3]][snp[0]-1], snp[1][i], str(snp[2][i]), '99', str(len(snp[1])), snp[4]]
                row = [snp[3], str(snp[0]), snp[5], snp[1][i], str(snp[2][i]), '99', str(len(snp[1])), snp[4]]
                f.write('\t'.join(row) + '\n')

    print('0-indexed: %d / %d' % (count_dup0, len(S)))
    print('1-indexed: %d / %d' % (count_dup1, len(S)))


def write_individuals(individuals, filename):
    with open(filename, 'w') as f:
        for ind in individuals:
            f.write(','.join(ind) + '\n')

def vcf_to_1ksnp(args):
    G = read_genome(args.reference)

    if args.ingroup and args.outgroup:
        print('Warning: ingroup and outgroup arguments are mutually exclusive: using only ingroup file')

    ingroup = None
    if args.ingroup:
        ingroup = read_filters(args.ingroup)

    outgroup = None
    if args.outgroup:
        outgroup = read_filters(args.outgroup)

    S = parse_vcf(args.vcf, args.individuals, ingroup, outgroup, args.indiv_hap, args.include_indels)
    print('Found %d SNPs' % len(S))

    write_1ksnp(S, G, args.out)

if __name__ == '__main__':

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--reference', type=str, required=True,
        help='Path to fasta file containing reference genome')
    parser.add_argument("--vcf", type=str, required=True, help="Path to VCF file containing mutation information")
    parser.add_argument("--ingroup", type=str, required=False, help="Optional path to file containing individuals to include, one on each line")
    parser.add_argument("--outgroup", type=str, required=False, help="Optional path to file containing individuals to exclude, one on each line")
    parser.add_argument("--out", type=str, required=True, help="Path to output file")
    parser.add_argument("--individuals", type=str, help="If present, write alleles for each individual to this file")
    parser.add_argument("--indiv-hap", type=str, help="Name of individual for which to write haplotype SNPs to file")
    parser.add_argument("--include-indels", action="store_true", help="Process both SNPs and INDELs")

    args = parser.parse_args(sys.argv[1:])

    vcf_to_1ksnp(args)
