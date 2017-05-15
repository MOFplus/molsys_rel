# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 18:09:01 2017

@author: rochus

          aftype

          a class for an aftype (atomtype and fragmenttype)

"""

# generates the missing rich comparison methods
from functools import total_ordering

@total_ordering
class aftype(object):

    def __init__(self, atype, fragtype):
        self.atype = atype
        self.fragtype = fragtype
        return

    @property
    def atype(self):
        return self._atype

    @atype.setter
    def atype(self, at):
        self._atype = at
        if not "_" in at:
            self._truncated = True
            self._atype_trunc = at
        else:
            self._truncated = False
            self._atype_trunc = at.split("_")[0]
        return

    @property
    def fragtype(self):
        return self._fragtype

    @fragtype.setter
    def fragtype(self, ft):
        self._fragtype = ft

    def __repr__(self):
        return "%s@%s" % (self._atype, self._fragtype)

    def __eq__(self, other):
        assert type(other) is aftype
        if self._truncated or other._truncated:
            return (self._atype_trunc == other._atype_trunc) and (self._fragtype == other._fragtype)
        else:
            return (self._atype == other._atype) and (self._fragtype == other._fragtype)

    def __lt__(self, other):
        assert type(other) is aftype
        return ("%s@%s" % (self._atype, self._fragtype)) < ("%s@%s" % (other._atype, other._fragtype))

    def __gt__(self, other):
        assert type(other) is aftype
        return ("%s@%s" % (self._atype, self._fragtype)) > ("%s@%s" % (other._atype, other._fragtype))



def aftype_sort(afl, ic):
    """
    helper function to sort a list of aftype objects according to the type (ic)
    """
    if ic == "bnd":
        afl.sort()
    elif ic == "ang":
        if afl[0] > afl[2]: afl.reverse()
    elif ic == "dih":
        if afl[1] > afl[2]:
            afl.reverse()
        elif afl[1] == afl[2]:
            if afl[0] > afl[3]: afl.reverse()
    elif ic == "oop":
        plane = afl[1:]
        plane.sort()
        afl[1:] = plane
    return afl

class afdict(object):
    """
    this is a "pseudo" dicitionary using two lists for key and value
    the reason is that the keys can not be hashed since strange comparsion is used.
    this is definetly less efficient than a real dictionary but carries only a few entries
    and is just used to store parameters for easy lookup

    one important limitation: you can only set a key that does not exist!
    so to change an entry you need to delete it first and then set it new.
    this is to prevent having multiple keys giving the same value.

    with apenditem it is however possible to append to an existing value (must be appendable)

    further more: if you use truncated aftypes during lookup it is possible that more than
    one full key in the afdict matches. only the first stored will be returned
    """


    def __init__(self):
        self._keys = []
        self._values = []
        return

    def __setitem__(self, key, value):
        if key in self._keys:
            raise KeyError, "key %s exists in afdict" % str(key)
        self._keys.append(key)
        self._values.append(value)
        return

    def appenditem(self, key, item):
        assert key in self
        assert type(self[key]) == type(list())
        self[key].append(item)
        return

    def __getitem__(self, key):
        try:
            idx = self._keys.index(key)
        except ValueError:
            raise KeyError, "key %s not in afdict" % str(key)
        return self._values[idx]

    def __contains__(self, item):
        return item in self._keys

    def __repr__(self):
        maxlen = 0
        keystring = []
        for k in self._keys:
            ks = str(k)
            if len(ks)> maxlen: maxlen = len(ks)
            keystring.append(ks)
        form = "%%-%ds = %%s\n" % (maxlen+3)
        out = "\n"
        for k,v in zip(keystring, self._values):
            out += form%(k,v)
        return out



if __name__ == "__main__":
    a = aftype("c3_c3", "ph")
    b = aftype("c3_c2h1", "ph")
    c = aftype("c3", "ph")
    d = aftype("c3", "co2")

    print a == b
    print a == c
    print a == d

    l = [a,b,c]
    l.sort
    print l

    print aftype_sort(l, "ang")
    exit()

    print "tuple comparison"
    t1 = (a,b)
    t2 = (a,c)
    t3 = (c,c)
    print t1 == t2
    print t1 == t3

    print "test afdict"
    # afd is non-political and means afdict!!!
    afd = afdict()

    afd[t1] = [str(t1)]
    afd[(a,d)] = [str((a,d,))]

    print afd[t1]
    print afd[t3]
    print (c,c) in afd
    afd.appenditem(t1, "test")

    print afd




