#! /usr/bin/env python3

import sys
import argparse

# Update a genome with a set of SNPs and write to a new file


def get_mutation_type(info):
    '''
        Return the value of the VT attribute from the info field
    '''

    attrs = info.split(';')
    for a in attrs:
        if a[:3] == 'VT=':
            return a[3:]


def update_genome(indiv, seq, label, vcf, out_prefix, indels=None):
    hapA = list(seq[:])
    hapB = list(seq[:])
    fA = open(out_prefix + '_hapA.fa', 'w')
    fA.write(label)
    fB = open(out_prefix + '_hapB.fa', 'w')
    fB.write(label)

    f_vars = open('indiv_varsA.txt', 'w')

    if indels:
        f_indelsA = open(indels+'A.txt', 'w')
        f_indelsB = open(indels+'B.txt', 'w')

    with open(vcf, 'r') as f:
        labels = None

        line_id = 0

        offset = 0
        for line in f:
            # Skip header lines
            if line[0] == '#' and line[1] == '#':
                continue

            if not labels:
                labels = line.rstrip().split('\t')

                col = None
                for i in range(9, len(labels)):
                    if labels[i] == indiv:
                        col = i
                if not col:
                    print('Error! Couldn\'t find individual %s in VCF' % indiv)
                    exit()
            else:
                row = line.rstrip().split('\t')
                type = get_mutation_type(row[7])
                if type == 'SNP' or (indels and type == 'INDEL'):
                    chrom = row[0]
                    name = row[2]
                    loc = int(row[1])

                    orig = row[3]
                    alts = row[4].split(',')

                    alleleA = int(row[col][0])
                    alleleB = int(row[col][2])

                    if alleleA > 0:
                        if indels:
                            if type == 'INDEL':
                                origLen = len(orig)
                                altLen = len(alts[alleleA-1])
                                if origLen > altLen:
                                    f_indelsA.write(chrom + '\t' + str(loc+altLen) + '\t' + str(altLen-origLen) + '\n')
                                    f_vars.write(name + '\tdeletion\t' + chrom + '\t' + str(loc-1+origLen) + '\t' + str(origLen - altLen) + '\n')
                                elif len(alts[alleleA-1]) > 1:
                                    f_indelsA.write(chrom + '\t' + str(loc+origLen) + '\t' + str(altLen-origLen) + '\n')
                                    f_vars.write(name + '\tinsertion\t' + chrom + '\t' + str(loc-1+origLen) + '\t' + alts[alleleA-1][origLen:] + '\n')
                                else:
                                    print(orig)
                                    print(alts[allelleA-1])
                                    exit()
                            else:
                                f_vars.write(name + '\tsingle\t' + chrom + '\t' + str(loc-1) + '\t' + alts[alleleA-1] + '\n')
                            offset = add_alt(hapA, loc-1, orig, alts[alleleA-1], offset)
                        else:
                            hapA[loc+offset-1] = alts[alleleA-1]
                    if alleleB > 0:
                        if indels:
                            if type == 'INDEL':
                                origLen = len(orig)
                                altLen = len(alts[alleleB-1])
                                if origLen > altLen:
                                    f_indelsB.write(chrom + '\t' + str(loc+altLen) + '\t' + str(altLen-origLen) + '\n')
                                elif len(alts[alleleB-1]) > 1:
                                    f_indelsB.write(chrom + '\t' + str(loc+origLen) + '\t' + str(altLen-origLen) + '\n')
                                else:
                                    print(orig)
                                    print(alts[allelleB-1])
                                    exit()
                            offset = add_alt(hapB, loc-1, orig, alts[alleleB-1], offset)
                        else:
                            hapB[loc+offset-1] = alts[alleleB-1]

                    line_id += 1

    if indels:
        f_indelsA.close()
        f_indelsB.close()
    f_vars.close()

    for i in range(0, len(seq), 60):
        fA.write(''.join(hapA[i:i+60]) + '\n')
        fB.write(''.join(hapB[i:i+60]) + '\n') 

    fA.close()
    fB.close()

def add_alt(genome, loc, orig, alt, offset):
    loc += offset

    if len(orig) == 1 and len(alt) == 1:
        # SNP
        genome[loc] = alt
    elif len(orig) > len(alt):
        # Deletion
        for i in range(len(alt)):
            genome[loc+i] = alt[i]
        del genome[loc+len(alt):loc+len(orig)]
        offset -= (len(orig) - len(alt))
    elif len(alt) > len(orig):
        # Insertion
        for i in range(len(orig)):
            genome[loc+i] = alt[i]
        genome[loc+len(orig):loc+len(orig)] = list(alt[len(orig):])
        offset += len(alt) - len(orig)
    else:
        # len(orig)=len(alt)>1 : Weird combination of SNP/In/Del
        for i in range(len(alt)):
            genome[loc+i] = alt[i]

    return offset

def read_chrom(ref, chrom):
    with open(ref, 'r') as f_in:
        label = None
        seq = ''
        for line in f_in:
            if line[0] == '>':
                if label:
                    return label, seq

                curr_chrom = line[1:].split(' ')[0]
                if curr_chrom == chrom:
                    label = line
            elif label:
                seq += line.rstrip()

        return label, seq

if __name__ == '__main__':

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ref', type=str, required=True, help='Path to fasta file containing reference genome')
    parser.add_argument("--vcf", type=str, required=True, help="Path to VCF file containing mutation information")
    parser.add_argument("--chrom", type=str, required=True, help="Chromosome to process")
    parser.add_argument("--out-prefix", type=str, required=True, help="Path to output prefix")
    parser.add_argument("--name", type=str, help="Name of individual in VCF to process")
    parser.add_argument("--include-indels", type=str, help="If present, extract both SNPs and INDELs. Argument is file to write list of INDELs included")

    args = parser.parse_args(sys.argv[1:])

    label,genome = read_chrom(args.ref, args.chrom)
    update_genome(args.name, genome, label, args.vcf, args.out_prefix, args.include_indels)

