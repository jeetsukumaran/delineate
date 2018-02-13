#! /usr/bin/env python

import dendropy

class LineageEdge(dendropy.Edge):

    def __init__(self, **kwargs):
        dendropy.Edge.__init__(self, **kwargs)

class LineageNode(dendropy.Node):

    def edge_factory(cls, **kwargs):
        return LineageEdge(**kwargs)

    def __init__(self, **kwargs):
        dendropy.Node.__init__(self, **kwargs)

class LineageTree(dendropy.Tree):

    def node_factory(cls, **kwargs):
        return LineageNode(**kwargs)
