# -*- coding: utf-8 -*-

try:
    import graph_tool
except ImportError:
    graph = "None"
else:
    from graph import graph
__all__=["graph"]
