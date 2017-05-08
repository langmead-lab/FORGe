'''
Utility functions
'''

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

    def __init__(self, chrom_seq, vars, r):
        self.seq = chrom_seq
        self.vars = vars
        self.r = r
        self.k = len(vars)
        self.counts = [v.num_alts for v in self.vars]

        self.vec = [1] + [0] * (self.k-1)

        self.read_chunks = []
        for i in range(1, len(vars)):
            self.read_chunks.append(self.seq[vars[i-1].pos + len(vars[i-1].orig) : vars[i].pos])

    def next(self, var_id, debug=False):
        if not self.vec:
            return None

        if debug:
            print(self.vec)

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
            read += self.seq[pos : pos + end_offset - rlen]
        else:
            read = read[:end_offset]

        start_offset = last_alt_base - self.r
        if start_offset < 0:
            pos = self.vars[0].pos
            read = self.seq[pos+start_offset:pos] + read
            if debug:
                print(len(self.seq))
                print(self.seq[92110:92170])
                print(read)
        else:
            read = read[start_offset:]
        self.start = self.vars[0].pos + start_offset

        self.vec = get_next_vector(self.k, self.counts, self.vec)

        return read

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

            
