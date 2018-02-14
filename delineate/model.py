#! /usr/bin/env python

import random
import math
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
        self._probability_of_no_speciation = None
        self._log_probability_of_no_speciation = None
        self._probability_of_any_speciation = None
        self._log_probability_of_any_speciation = None
        dendropy.Edge.__init__(self, **kwargs)

    def _get_log_probability_of_no_speciation(self):
        if self._log_probability_of_no_speciation is None:
            self._log_probability_of_no_speciation = -1.0 * self.tree.species_conversion_rate.value * self.length
        return self._log_probability_of_no_speciation
    log_probability_of_no_speciation = property(_get_log_probability_of_no_speciation)

    def _get_probability_of_no_speciation(self):
        if self._probability_of_no_speciation is None:
            self._probability_of_no_speciation = math.exp(self.log_probability_of_no_speciation)
        return self._probability_of_no_speciation
    probability_of_no_speciation = property(_get_probability_of_no_speciation)

    def _get_probability_of_any_speciation(self):
        if self._probability_of_any_speciation is None:
            self._probability_of_any_speciation = 1.0 - self.probability_of_no_speciation
        return self._probability_of_any_speciation
    probability_of_any_speciation = property(_get_probability_of_any_speciation)

    def _get_log_probability_of_any_speciation(self):
        if self._log_probability_of_any_speciation is None:
            self._log_probability_of_any_speciation = math.log(self.probability_of_any_speciation)
        return self._log_probability_of_any_speciation
    log_probability_of_any_speciation = property(_get_log_probability_of_any_speciation)

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

    def build(self, *args, **kwargs):
        self.rng = kwargs.pop("rng", None)
        if self.rng is None:
            self.rng = random.Random(kwargs.pop("random_seed", None))
        self.species_conversion_rate = kwargs.pop("species_conversion_rate", None)
        if self.species_conversion_rate is None:
            self.species_conversion_rate = GammaDistributedParameter(
                    rng=self.rng,
                    shape=1.0,
                    scale=0.1)


