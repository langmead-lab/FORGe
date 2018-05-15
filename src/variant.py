#! /usr/bin/env python2.7

class VariantSet:
    """
    Columnar storage seems to keep memory footprint lower.  Also, since common
    case is to have 1 ALT, we keep that ALT in the alt1s and prob1s lists.
    The altrests and probrests arrays have non-None entries only for variants
    that have more than 1 ALT.

    Another idea is to serialize chromosomes to disk until they're needed
    again.
    """

    def __init__(self):
        self.poss = []
        self.origs = []
        self.alt1s = []
        self.prob1s = []
        self.altrests = []
        self.probrests = []

    def __len__(self):
        return len(self.poss)

    def add_var(self, pos, orig, alt, prob):
        intern(alt)
        intern(orig)
        new_id = len(self.poss)
        assert orig != alt
        assert prob <= 1.01
        self.poss.append(pos)
        self.origs.append(orig)
        self.alt1s.append(alt)
        self.prob1s.append(prob)
        self.altrests.append(None)
        self.probrests.append(None)
        return new_id

    def to_string(self, var_id):
        alt_strs = ['%s:%0.4f' % (self.alt1s[var_id], self.prob1s[var_id])]
        if self.altrests[var_id] is not None:
            alts, probs = self.altrests[var_id], self.probrests[var_id]
            for i, alt in enumerate(alts):
                alt_strs.append('%s:%0.4f' % (alt, probs[i]))
        return ('%d %s ' % (self.poss[var_id], self.origs[var_id])) + ' '.join(alt_strs)

    def add_alt_to_last(self, alt, prob):
        intern(alt)
        if len(self.poss) == 0:
            raise RuntimeError('No variants in set')
        if self.altrests[-1] is None:
            self.altrests[-1] = [alt]
            self.probrests[-1] = [prob]
        else:
            self.altrests[-1].append(alt)
            self.probrests[-1].append(prob)
        assert self.sum_probs(-1) <= 1.01, self.to_string(-1)

    def get_alt(self, var_id, alt_id):
        if var_id >= len(self.poss):
            raise RuntimeError('var_id %d exceeds # variants %d' % (var_id, len(self.poss)))
        if alt_id == 0:
            return self.alt1s[var_id]
        assert self.altrests[var_id] is not None, (var_id, alt_id)
        return self.altrests[var_id][alt_id-1]

    def get_prob(self, var_id, alt_id):
        if var_id >= len(self.poss):
            raise RuntimeError('var_id %d exceeds # variants %d' % (var_id, len(self.poss)))
        if alt_id == 0:
            return self.prob1s[var_id]
        assert self.probrests[var_id] is not None
        return self.probrests[var_id][alt_id-1]

    def num_alts(self, var_id):
        return 1 + (0 if self.altrests[var_id] is None else len(self.altrests[var_id]))

    def sum_probs(self, var_id):
        ret = self.prob1s[var_id]
        if self.probrests[var_id] is not None:
            tot = ret + sum(self.probrests[var_id])
            assert tot <= 1.01, self.to_string(-1)
            return tot
        return ret
