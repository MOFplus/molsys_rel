### COLOR DICTIONARIES ###
# based on default molden element colors
ecolor2elem = [
    "b" ,"f" ,"n" ,"o" ,"c" ,"he","ne","ge","li","s" ,"cl","p" ,"al","si",
]
elem2ecolor = dict(ke[::-1] for ke in enumerate(ecolor2elem))
maxecolor = len(ecolor2elem)
vcolor2elem = [
    "n" ,"o" ,"b" ,"f" ,"c" ,"he","ne","ge","li","s" ,"cl","p" ,"si","al",
]
elem2vcolor = dict(kv[::-1] for kv in enumerate(vcolor2elem))
maxvcolor = len(vcolor2elem)

# string conversion tools: elem+atype+color <-> string #
def elematypecolor2string(elem, atype, color):
    """
    return formatted string from element, atomtype, and color
    """
    return "%s_%s/%s" % (elem, atype, color)

eac2str = elematypecolor2string #nickname

def string2elematypecolor(st):
    """
    return element, atomtype, and color from a given formatted string.

    Color is after the last slash
    Element is before the first underscore
    Atype is everything in the middle
    """
    self.assert_eacstr(st)
    colorsplit = st.split("/")
    rest, color = colorsplit[:-1], colorsplit[-1]
    rest = "".join(rest)
    elemsplit = rest.split("_")
    elem, atype = elemsplit[0], elemsplit[1:]
    return elem, atype, color

str2eac = string2elematypecolor #nickname

def assert_eacstr(st):
    """
    check format of element_atomtype/color string
    """
    assert st.index("_"), "string is not well-formatted: no \"_\" found"
    assert st.rindex("/"), "string is not well-formatted: no \"/\" found"
    assert st.index("_") < st.rindex("/"), "string is not well-formatted: first \"_\" after last \"/\""
    return

