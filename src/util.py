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

