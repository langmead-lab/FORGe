#!/usr/bin/env python

import sys


def go():
    in_sam = sys.argv[1]
    in_indels = sys.argv[2]
    out_sam = sys.argv[3]

    # Map individual position to reference position
    '''
    max_len = 60000000
    indiv_to_ref = [0] * max_len
    indiv_i = 0
    ref_i = 0
    with open(in_indels, 'r') as f:
        for line in f:
            row = line.rstrip().split('\t')
            pos = int(row[1])
            length = -int(row[2])
            while indiv_i < pos:
                indiv_to_ref[indiv_i] = ref_i
                ref_i += 1
                indiv_i += 1
            if length > 0:
                for i in range(length):
                    indiv_to_ref[indiv_i] = ref_i-1
                    indiv_i += 1
            else:
                ref_i -= length

        while indiv_i < max_len:
            indiv_to_ref[indiv_i] = ref_i
            ref_i += 1
            indiv_i += 1
    '''

    max_len = 60000000
    indiv_to_ref = [0] * max_len
    indiv_i = 0
    ref_i = 0

    indels_ref = []
    indels_indiv = []
    with open(in_indels, 'r') as f:
        for line in f:
            row = line.rstrip().split('\t')
            pos = int(row[1])+1
            length = -int(row[2])
            while ref_i < pos:
                indiv_to_ref[indiv_i] = ref_i
                ref_i += 1
                indiv_i += 1
            if length > 0:
                indels_indiv.append(ref_i-1)
                for i in range(length):
                    indels_ref.append(indiv_i)
                    indiv_to_ref[indiv_i] = ref_i-1
                    indiv_i += 1
            else:
                indels_ref.append(indiv_i)
                for i in range(-length):
                    indels_indiv.append(ref_i+i)
                ref_i -= length

        while indiv_i < max_len:
            indiv_to_ref[indiv_i] = ref_i
            ref_i += 1
            indiv_i += 1

    '''
    a = 49188557
    print(a)
    for i in range(len(indels_indiv)):
        if indels_indiv[i] > a:
            print(indels_indiv[i-3:i+3])
            break
    print('')

    b = 4959275
    print(b)
    for i in range(len(indels_ref)):
        if indels_ref[i] > b:
            print(indels_ref[i-3:i+3])
            break
    print('')
    exit()
    '''

    '''
    indels = []
    with open(in_indels, 'r') as f:
        for line in f:
            row = line.rstrip().split('\t')
            indels.append((int(row[1]), int(row[2])))

    def pos_offset(a):
        off = 0
        for i in indels:
            if i[0] > a:
                return off
            else:
                off += i[1]
        return off

    def calc_offsets():
        max_len = 60000000
        offsets = [0] * max_len
        curr_offset = 0
        curr_i = 0
        max_i = len(indels)
        for i in range(max_len):
            while curr_i < max_i and indels[curr_i][0] <= i:
                curr_offset += indels[curr_i][1]
                curr_i += 1
            offsets[i] = curr_offset

        return offsets

    offsets = calc_offsets()
    '''

    with open(in_sam, 'r') as f_in:
        with open(out_sam, 'w') as f_out:
            for line in f_in:
                if line[0] == '@':
                    f_out.write(line)
                    continue
                row = line.rstrip().split('\t')
                chrom = row[2]
                row[3] = str(indiv_to_ref[int(row[3])]-1)
                #pos = int(row[3])
                #row[3] = str(pos + offsets[pos] - 1)
                f_out.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    go()
