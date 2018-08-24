# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 18:09:01 2017

@author: rochus

          aftype

          a class for an aftype (atomtype and fragmenttype)

"""

# generates the missing rich comparison methods
from functools import total_ordering
# import regular expressions module
import re

@total_ordering
class aftype(object):

    def __init__(self, atype, fragtype):
        # boolean for * fragtype
        self._wild_ft = False
        # boolean for * atype
        self._wild_at = False
        # boolean for pure elem atype
        self._pure = False
        # boolean for truncated atype
        self._truncated = False
        self.atype = atype
        self.fragtype = fragtype
        return

    @property
    def atype(self):
        return self._atype

    @atype.setter
    def atype(self, at):
        self._atype = at
        match = re.search("[0123456789]", at)
        if at == "*":
            self._wild_at = True
        # check for trucated type like c3
        elif not "_" in at:
            # check for pure elem type like c
            if not match:
                # we have a pure elem type
                self._pure = True
                self._atype_pure = at
                self._atype_trunc = at
            else:
                self._truncated = True
                self._atype_trunc = at
                self._atype_pure = at[0:match.start()]
        else:
            self._atype_trunc = at.split("_")[0]
            self._atype_pure = at[0:match.start()]
        return

    @property
    def fragtype(self):
        return self._fragtype

    @fragtype.setter
    def fragtype(self, ft):
        self._fragtype = ft
        if ft == "*":
            self._wild_ft == True

    def __repr__(self):
        return "%s@%s" % (self._atype, self._fragtype)

    # comparison methods for all possibilities
    def full_compare(self,other):
        return (self._atype == other._atype) and (self._fragtype == other._fragtype)

    def trunc_compare(self, other):
        if self._truncated or other._truncated:
            return (self._atype_trunc == other._atype_trunc) and (self._fragtype == other._fragtype)
        else:
            return self.full_compare(other)

    def pure_compare(self,other):
        if self._pure or other._pure:
            return (self._atype_pure == other._atype_pure) and (self._fragtype == other._fragtype)
        else:
            return self.trunc_compare(other)

    def wildat_compare(self,other):
        if self._wild_at or other._wild_at:
            return self._fragtype == other._fragtype
        else:
            return self.pure_compare(other)
    
    def wildft_compare(self,other):
        if self.fragtype == "*" or other._wild_ft == "*":
            if self._wild_at or other._wild_at:
                return True
            elif self._pure or other._pure:
                return self._atype_pure == other._atype_pure
            elif self._truncated or other._truncated:
                return self._atype_trunc == other._atype_trunc
            else:
                return self._atype == other._atype          
        else:
            return self.wildat_compare(other)
  
    def __eq__(self,other, level = "full"):
        if level == "full":
            return self.full_compare(other)
        elif level == "trunc":
            return self.trunc_compare(other)
        elif level == "pure":
            return self.pure_compare(other)
        elif level == "wild_at":
            return self.wildat_compare(other)
        elif level == "wild_ft":
            return self.wildft_compare(other)

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

    def index(self, key, comp_level = "full"):
        # we need to write a custom search/index method for the _keys list
        # since we have to use different comparison leves because of
        # wildcards
        # first loop over keys, list of tuples of aftypes
        for i,k in enumerate(self._keys):
            # now check if k and key have the same length
            if len(k) == len(key):
                # now we have to compare in a second loop each entry of
                # key and k in respect to self._comp_level
                for j in range(len(k)):
                    if not k[j].__eq__(key[j], comp_level):
                        break
                    elif j+1 == len(k):
                        # found it
                        return i
        # have not found it
        return -1


    def __setitem__(self, key, value):
        # in principle we have to check already here if the given key
        # contains wildcards, then we have to in pri
        idx = self.index(key)
        if idx >= 0:
            raise KeyError("key %s exists in afdict" % str(key))
        self._keys.append(key)
        self._values.append(value)
        return

    def appenditem(self, key, item):
        idx = self.index(key)
        if idx >= 0:
            assert type(self._values[idx]) == type(list())
            self._values[idx].append(item)
        else:
            raise KeyError("key %s not in afdict" % str(key))
        return

    def __getitem__(self, key, wildcards = False):
        if wildcards:
            idx = self.index(key, "full")
            if idx >= 0:
                return self._values[idx]
            else:
                idx = self.index(key, "trunc")
                if idx >= 0:
                    return self._values[idx]
                else:
                    idx = self.index(key, "pure")
                    if idx >= 0:
                        return self._values[idx]
                    else:
                        idx = self.index(key, "wild_at")
                        if idx >= 0:
                            return self._values[idx]
                        else:
                            idx = self.index(key, "wild_ft")
                            if idx >= 0:
                                return self._values[idx]
                            else:
                                raise KeyError("key %s not in afdict" % str(key))
        else:
            idx = self.index(key)
            if idx >= 0:
                return self._values[idx]
            else:
                raise KeyError("key %s not in afdict" % str(key))

    def __contains__(self, item, wildcards = False):
        raise NotImplementedError("Do not use this method -> use index() instead!")

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
    e = aftype("*", "*")

    print(a == b)
    print(a == c)
    print(a == d)

    l = [a,b,c]
    l.sort
    print(l)

    print(aftype_sort(l, "ang"))
    #exit()

    print("tuple comparison")
    t1 = (a,b)
    t2 = (a,c)
    t3 = (c,c)
    print(t1 == t2)
    print(t1 == t3)

    print("test afdict")
    # afd is non-political and means afdict!!!
    afd = afdict()

    afd[t1] = [str(t1)]
    afd[(a,d)] = [str((a,d,))]

    print(afd[t1])
    #print(afd[t3])
    #print(c,c) in afd
    afd.appenditem(t1, "test")

    print(afd)




