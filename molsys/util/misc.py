from functools import cmp_to_key
import string
digs = string.ascii_uppercase # max base == 26

def argsorted(seq, cmp=None, key=None, reverse=False, sort_flag=False):
    """Return the index that would sort a sequence. (python2.7 fashion)

    Instead of numpy.argsort, it allows to define optional arguments of sorted
    seq: the sequence to be argsorted
    cmp (function): custom comparison of two arguments of the sequence which should
        return negative/zero/positive number whether the first argument is
        considered lower/equal to/greater than the second
    key (function): function used to extract a comparison key from each list element
    reverse (bool): as if sorting had reverse comparison
    sort_flag (bool): if True: sort sequence in place as the standard `sort` method
    """
    if key is None:
        try:
            argsorted = sorted(
                range(len(seq)), cmp=cmp, reverse=reverse,
                key=seq.__getitem__)
        except TypeError as e: ### python3
            if cmp is not None:
                raise TypeError(e)
            argsorted = sorted(
                range(len(seq)), reverse=reverse,
                key=seq.__getitem__)
    else:
        try:
            argsorted = sorted(
                range(len(seq)), cmp=cmp, reverse=reverse,
                key=lambda x: key(seq.__getitem__(x)) )
        except TypeError as e: ### python3
            if cmp is not None:
                raise TypeError(e)
            argsorted = sorted(
                range(len(seq)), reverse=reverse,
                key=lambda x: key(seq.__getitem__(x)) )
    if sort_flag:
        try:
            seq.sort(cmp=cmp, key=key, reverse=reverse)
        except TypeError as e: ### python3
            if cmp is not None:
                raise TypeError(e)
            seq.sort(key=key, reverse=reverse)
    return argsorted

def normalize_ratio(cratio, total):
    """ TBI: update documentation! this is not only for colors! [RA]
    return normalized color ratio so that:
        the total number of elements (edges or vertices) is colored and
        the actual (float) ratio among colors is close to the given (int) ratio.
    It implements D'Hondt's quotients internally. Credits: https://github.com/rg3
    N.B.: in case of "ties", the latter color gets the higher ratio of
    elements. (this is consistent AND reproducible!)
    
    :Parameters:
     - cratio (list of ints): color ratio
     - total (int): total number of elements to color with given ratio
    :Returns:
     - norm_cratio (list of ints): color ratio normalized wrt. elements
    """
    def dhondt_quotient(cratio, subtotal):
        return float(cratio) / (subtotal + 1)

    # calculate the quotients matrix (list in this case)
    dict_cratio = dict(enumerate(cratio))
    quot = []
    ret ={}
    for p in dict_cratio:
        ret[p] = 0
        for s in range(0, total):
            q = dhondt_quotient(dict_cratio[p], s)
            quot.append((q, p))

    # sort the quotients by value
    quot.sort(reverse=True)

    # take the highest quotients with the assigned parties
    for s in range(0, total):
        ret[quot[s][1]] += 1
    norm_cratio = ret.values()
    return norm_cratio # already ordered

def triplenats_on_sphere(trisum, trimin=1):
    """returns triplets of natural numbers on a sphere
    trisum(int):the summation of the triples must be equal to trisum
    trimin(int):minimum allowed natural per triplet element (default: 1)"""
    trinat = []
    for itri in itertools.product(range(trimin, trisum), repeat=3):
        if sum(itri) == trisum:
            trinat.append(itri)
    return trinat

def int2base(value, base=None, maximum=None):
    """credits: A. Martelli"""
    if base is None:
        base = len(digs)
    if maximum is not None:
        return int2base_spaced(value, base, maximum)
    assert 0 <= base <= len(digs), "Not enough digits"
    if value < 0:
        sign = -1
    elif value == 0:
        return digs[0]
    else:
        sign = 1

    value *= sign
    digits = []

    while value:
        digits.append(digs[int(value % base)])
        value = int(value / base)

    if sign < 0:
        digits.append('-')

    digits.reverse()

    return ''.join(digits)

def int2base_spaced(value, base=None, maximum=None):
    if base is None:
        base = len(digs)
    if maximum is None:
        return int2base(value, base)
    assert base <= len(digs), "Not enough digits"

    vb = int2base(value, base)
    mb = int2base(maximum-1, base)

    addl = max([len(mb) - len(vb), 0])
    add = digs[0]*addl

    return add+vb

