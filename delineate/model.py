#! /usr/bin/env python

import random
import dendropy

class ContinuousParameter(object):
    def __init__(self, *args, **kwargs):
        pass

class GammaDistributedParameter(ContinuousParameter):

    def __init__(self, rng, shape, scale, initial_value=None):
        self.rng = rng
        self.shape = shape
        self.scale = scale
        if initial_value is not None:
            self.value = initial_value
        else:
            self.value = self.draw()

    def draw(self):
        return self.rng.gammavariate(self.shape, self.scale)

class LineageEdge(dendropy.Edge):

    def __init__(self, tree, **kwargs):
        self.tree = tree
        dendropy.Edge.__init__(self, **kwargs)

class LineageNode(dendropy.Node):

    def __init__(self, tree, **kwargs):
        self.tree = tree
        dendropy.Node.__init__(self, **kwargs)

    def edge_factory(self, **kwargs):
        return LineageEdge(tree=self.tree, **kwargs)

class LineageTree(dendropy.Tree):

    def __init__(self, *args, **kwargs):
        dendropy.Tree.__init__(self, *args, **kwargs)

    def node_factory(self, *args, **kwargs):
        return LineageNode(tree=self, **kwargs)

    def parameterize(self, *args, **kwargs):
        self.rng = kwargs.pop("rng", None)
        if self.rng is None:
            self.rng = random.Random(kwargs.pop("random_seed", None))
        self.species_conversion_rate = kwargs.pop("species_conversion_rate", None)
        if self.species_conversion_rate is None:
            self.species_conversion_rate = GammaDistributedParameter(
                    rng=self.rng,
                    shape=1.0,
                    scale=0.1)


