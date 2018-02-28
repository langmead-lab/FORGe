#! /usr/bin/env python2.7

class Variant:
    def __init__(self, first_var, counts, vec):
        if not len(counts) == len(vec):
            print('Error in pseudocontig inititaion! Counts and vector length do not match!')

        self.first = first_var
        self.counts = counts
        self.vec = vec
        self.vec_id = self.vec_to_id(vec, counts)
        self.num_v = len(counts)

    def vec_to_id(self, v, counts):
        id = 0
        for i in range(len(v)):
            if v[i] > counts[i]:
                print('Error in allele vector! Vector is ' + str(v) + ' but counts are ' + str(counts))
                traceback.print_stack()
                exit()

            id = id * (counts[i]+1) + v[i]

        return id

    def id_to_vec(self, id, counts):
        v = [0] * len(counts)
        for i in range(len(counts)-1, 0, -1):
            v[i] = id % (counts[i]+1)
            id = (id - v[i]) / (counts[i]+1)
        v[0] = id
        return v
