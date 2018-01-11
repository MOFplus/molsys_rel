# -*- coding: utf-8 -*-
from __future__ import absolute_import

try:
    import graph_tool
except ImportError:
    graph = None
else:
    from .graph import graph

try:
    import os
    import sys
    f = open(os.devnull, 'w') #cross-platform and version
    sys.stderr = f
    import pandas ### MODULE WITH ANNOYING WARNING MESSAGES
    import chemcoord
except:
    zmat = None
else:
    from .zmat import zmat
finally:
    sys.stderr = sys.__stderr__

try:
    import spglib
except:
    spg = None
else:
    from .spg import spg

try:
    import ff_gen.ric_new
except:
    ric = None
else:
    from .ric import ric

from .fragments import fragments
from .bb import bb
from .base import base
from .ff import ff
from .molecules import molecules

__all__=["graph", "fragments", "bb", "zmat", "spg", "ff", "molecules", "ric" "base"]
