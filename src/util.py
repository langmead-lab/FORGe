'''
Utility functions
'''


from variant import Variant


def get_next_vector(k, counts, v=None):
    '''
    Loop through all allele vectors in order
    counts: Array of length k storing the number of alternate alleles for each variant
    If v is None, return an array of all zeroes. If v is the last array, return None
    '''

    if not v:
        return [0] * k
    else:
        j = k-1
        while v[j] == counts[j]:
            if j == 0:
                return None

            v[j] = 0
            j -= 1
        v[j] += 1
        return v


def vec_to_id(v, counts):
    id = 0
    for i in range(len(v)):
        if v[i] > counts[i]:
            print('Error in allele vector! Vector is ' + str(v) + ' but counts are ' + str(counts))
            traceback.print_stack()
            exit()

        id = id * (counts[i]+1) + v[i]

    return id


class PseudocontigIterator:
    ''' Loop over all pseudocontigs that contain a certain set of variants '''

    def __init__(self, chrom_seq, vars, r, vec=None):
        self.seq = chrom_seq
        self.vars = vars
        self.r = r
        self.k = len(vars)
        self.counts = [v.num_alts for v in self.vars]

        if vec:
            self.vec = vec
        else:
            self.vec = [1] + [0] * (self.k-1)

        self.read_chunks = []
        for i in range(1, len(vars)):
            self.read_chunks.append(self.seq[vars[i-1].pos + len(vars[i-1].orig) : vars[i].pos])

        self.valid = True
        for chunk in self.read_chunks:
            if 'N' in chunk or 'M' in chunk or 'R' in chunk:
                self.valid = False
                break

    def next(self, debug=False):
        '''
            Return the current read and iterate to the next one
        '''

        if not self.vec:
            return None

        if debug:
            print(self.vec)

        valid = True

        first_alt_base = -1
        last_alt_base = -1
        if self.vec[0] == 0:
            read = self.vars[0].orig
        else:
            first_alt_base = 0
            read = self.vars[0].alts[self.vec[0]-1]
            last_alt_base = len(read)
        for i in range(1, self.k):
            if self.vec[i] == 0:
                read += self.read_chunks[i-1] + self.vars[i].orig
            else:
                if first_alt_base < 0:
                    first_alt_base = len(read)
                read += self.read_chunks[i-1] + self.vars[i].alts[self.vec[i]-1]
                last_alt_base = len(read)

        if first_alt_base < 0:
            print('Error! Tried to process empty vector')
            exit()

        # trim/extend ends of read to the right length
        end_offset = first_alt_base + self.r
        if end_offset <= last_alt_base:
            end_offset = last_alt_base+1
        rlen = len(read)
        if end_offset > rlen:
            pos = self.vars[-1].pos + len(self.vars[-1].orig)
            suffix = self.seq[pos : pos+end_offset-rlen]
            if 'N' in suffix or 'M' in suffix or 'R' in suffix:
                valid = False
            read += suffix
        else:
            read = read[:end_offset]

        start_offset = last_alt_base - self.r
        if start_offset < 0:
            pos = self.vars[0].pos
            assert pos >= 0
            prefix = self.seq[max(pos + start_offset, 0):pos]
            if 'N' in prefix or 'M' in prefix or 'R' in prefix:
                valid = False

            read = prefix + read
        else:
            read = read[start_offset:]
        self.start = self.vars[0].pos + start_offset

        self.curr_vec = self.vec[:]
        
        self.vec = get_next_vector(self.k, self.counts, self.vec)

        if valid:
            return read
        else:
            return self.next()

class ReadIterator:
    ''' Loop over all reads in a genome '''

    def __init__(self, chrom, seq, vars, r):
        self.chrom = chrom
        self.seq = seq
        self.chrom_len = len(chrom_seq)
        self.vars = vars
        self.num_v = len(vars)
        self.r = r

        self.i = 0
        for j in range(self.num_v):
            if self.vars[j].chrom == chrom:
                break
        self.first_var = j
        pos = self.vars[self.first_var].pos
        while j < self.num_v and self.vars[j].chrom == chrom and self.vars[j].pos - self.i < r:
            j += 1
        self.last_var = j
        self.curr_vars = vars[self.first_var:self.last_var]
        self.vec = [0] * (self.last_var - self.first_var)

    def next(self):
        if self.i >= self.chrom_len:
            return None

        read = self.seq[self.i]

        # Update vector
        self.vec = get_next_vector(self.last_var-self.first_var, self.curr_vars, self.vec)
        if not self.vec:
            if self.vars[self.first_var].pos == self.i:
                self.first_var += 1
            self.i += 1
            if self.vars[self.last_var].chrom == self.chrom and self.vars[self.last_var].pos - self.i < self.r:
                self.last_var += 1
            self.curr_vars = vars[self.first_var:self.last_var]
            self.vec = [0] * (self.last_var - self.first_var)


def test_pc_iter_1():
    seq = 'AAAAAAAAA'
    #      012345678
    #          T

    variants = [Variant('t', 't', 4, 'A', 'T', 0.25)]

    it = PseudocontigIterator(seq, variants, 4)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 1 == len(pcs)
    assert 'AAATAAA' == pcs[0]

    it = PseudocontigIterator(seq, variants, 5)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 1 == len(pcs)
    assert 'AAAATAAAA' == pcs[0]

    it = PseudocontigIterator(seq, variants, 6)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 1 == len(pcs)
    assert 'AAAATAAAA' == pcs[0]

    it = PseudocontigIterator(seq, variants, 9)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 1 == len(pcs)
    assert 'AAAATAAAA' == pcs[0]


def test_pc_iter_2():
    seq = 'AAAAAAAAAA'
    #      0123456789
    #          CG

    variants = [Variant('t', 't', 4, 'A', 'C', 0.25),
                Variant('t', 't', 5, 'A', 'G', 0.25)]
    it = PseudocontigIterator(seq, variants, 4)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 2 == len(pcs)
    assert 'AAACAAA' in pcs
    assert 'AACGAA' in pcs

    it = PseudocontigIterator(seq, variants, 5)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 2 == len(pcs)
    assert 'AAAACAAAA' in pcs
    assert 'AAACGAAA' in pcs


def test_pc_iter_3():
    seq = 'AAAAAAAAAAA'
    #      01234567890
    #          CGT

    variants = [Variant('t', 't', 4, 'A', 'C', 0.25),
                Variant('t', 't', 5, 'A', 'G', 0.25),
                Variant('t', 't', 6, 'A', 'T', 0.25)]
    it = PseudocontigIterator(seq, variants, 4)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 'AAACAAA' in pcs
    assert 'AACGAA' in pcs
    assert 'ACATA' in pcs
    assert 'ACGTA' in pcs


def test_pc_iter_4():
    seq = 'AAANAAAAA'
    #      012345678
    #          T

    variants = [Variant('t', 't', 4, 'A', 'T', 0.25)]
    it = PseudocontigIterator(seq, variants, 4)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 0 == len(pcs)


def test_pc_iter_deletion_1():
    seq = 'AAAAAAAAA'
    #      012345678
    #          a

    variants = [Variant('t', 't', 4, 'A', [''], 0.25)]
    it = PseudocontigIterator(seq, variants, 4)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 1 == len(pcs)
    assert 'AAAAAAAA' == pcs[0]


def test_pc_iter_deletion_2():
    seq = 'AAAAAAAAA'
    #      012345678
    #         xxx

    variants = [Variant('t', 't', 3, 'AAA', [''], 0.25)]
    it = PseudocontigIterator(seq, variants, 4)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 1 == len(pcs)
    assert 'AAAAAA' == pcs[0]


def test_pc_iter_insertion_1():
    seq = 'AAAAAAAAA'
    #      012345678
    #          ^
    #          TT

    variants = [Variant('t', 't', 4, 'A', ['TT'], 0.25)]
    it = PseudocontigIterator(seq, variants, 4)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 1 == len(pcs)
    assert 'AATTAA' == pcs[0]


def test_pc_iter_insertion_2():
    seq = 'AAAAAAAAA'
    #      012345678
    #          ^
    #          TT

    variants = [Variant('t', 't', 4, '', ['TT'], 0.25)]
    it = PseudocontigIterator(seq, variants, 4)
    pc = it.next()
    pcs = []
    while pc:
        pcs.append(pc)
        pc = it.next()
    assert 1 == len(pcs)
    assert 'AATTAA' == pcs[0]
