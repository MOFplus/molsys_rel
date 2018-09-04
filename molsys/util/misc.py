import profilehooks

def argsorted(seq, cmp=None, key=None, reverse=False, sort_flag=False):
    if key is None:
        argsorted = sorted(
            range(len(seq)), cmp=cmp, reverse=reverse,
            key=seq.__getitem__)
    else:
        argsorted = sorted(
            range(len(seq)), cmp=cmp, reverse=reverse,
            key=lambda x: key(seq.__getitem__(x)) )
    if sort_flag:
        seq.sort(cmp=cmp, key=key, reverse=reverse)
    return argsorted
