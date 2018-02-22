#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
        # The probability that all of the tips that are descendants of this node have
        # had at least one "good" speciation event between themselves and this
        # node.
        self.algvar_x = None
        # The probability that none of the tips that are descendants of
        # this node have had a good speciation event separating them from this
        # node.
        self.algvar_y = None
        # Probability that all of the descendants of this node that are member
        # of the "good species" clade are separated from this node by 0 good
        # speciation events and all of the descendants of this node that are
        # not members of the good species clade are separated from this node by
        # at least 1 good speciation event.
        self.algvar_z = None
        # The status of node i with respect to whether its descendants are in
        # the focal clade or not:
        #   (a) 1=“ancestor of only members of c”,
        #   (b) 2=“ancestor of only non-members of c”, and
        #   (c) 3=“ancestor of both members of c and non-members of c”
        # where "c" = the clade of candidate good species.
        self.algvar_s = None
        dendropy.Node.__init__(self, **kwargs)

    def edge_factory(self, **kwargs):
        return LineageEdge(tree=self.tree, **kwargs)

    def calc_probability_of_good_species_clade(self, taxa):
        Z = 0.0
        if not self._child_nodes:
            self.algvar_x = 0.0
            self.algvar_y = 1.0
            if self.taxon in taxa:
                self.algvar_s = 1
            else:
                self.algvar_s = 2
        else:
            assert len(self._child_nodes) == 2
            ch1, ch2 = self._child_nodes
            if ch1.algvar_s == ch2.algvar_s:
                self.algvar_s = ch1.algvar_s
            else:
                self.algvar_s = 3
            if self.algvar_s == 1:
                # Set: yi = yr(i)yl(i)  (1 - pr(i))   (1 - pl(i))
                self.algvar_y = ch1.algvar_y * ch2.algvar_y * ch1.edge.probability_of_no_speciation * ch2.edge.probability_of_no_speciation
                assert self.algvar_y == ch1.algvar_y * ch2.algvar_y * (1 - ch1.edge.probability_of_any_speciation) * (1 - ch2.edge.probability_of_any_speciation) ## debug
            elif self.algvar_s == 2:
                ch1.algvar_t = ch1.edge.probability_of_any_speciation + (ch1.edge.probability_of_no_speciation) * ch1.algvar_x
                assert ch1.algvar_t == ch1.edge.probability_of_any_speciation + (1 - ch1.edge.probability_of_any_speciation) * ch1.algvar_x # debug
                ch2.algvar_t = ch2.edge.probability_of_any_speciation + (ch2.edge.probability_of_no_speciation) * ch2.algvar_x
                assert ch2.algvar_t == ch2.edge.probability_of_any_speciation + (1 - ch2.edge.probability_of_any_speciation) * ch2.algvar_x # debug
                self.algvar_x = ch1.algvar_t * ch2.algvar_t
            elif self.algvar_s == 3:
                if ch1.algvar_s != 3 and ch2.algvar_s != 3:
                    if ch1.algvar_s == 1:
                        assert ch2.algvar_s != 1
                        ch_f = ch1
                        ch_g = ch2
                    elif ch2.algvar_s == 1:
                        assert ch1.algvar_s != 1
                        ch_f = ch2
                        ch_g = ch1
                    else:
                        raise ValueError
                    ch_g.algvar_t = ch_g.edge.probability_of_any_speciation + ch_g.edge.probability_of_no_speciation * ch_g.algvar_x
                    ch_f.algvar_t = ch_f.edge.probability_of_no_speciation * ch_f.algvar_y
                    Z += (ch_f.algvar_y * ch_f.edge.probability_of_any_speciation)
                    self.algvar_z = ch_f.algvar_t * ch_g.algvar_t
                elif ch1.algvar_s == 3 and ch2.algvar_s == 3:
                    raise ValueError
                else:
                    assert (ch1.algvar_s == 3 or ch2.algvar_s == 3) and not (ch1.algvar_s == 3 and ch2.algvar_s == 3) # debug
                    if ch1.algvar_s == 3:
                        assert ch2.algvar_s != 3
                        ch_f = ch1
                        ch_g = ch2
                    elif ch2.algvar_s == 3:
                        assert ch1.algvar_s != 3
                        ch_f = ch2
                        ch_g = ch1
                    else:
                        raise ValueError
                    Z += (ch_f.algvar_z * ch_f.edge.probability_of_any_speciation)
                    ch_g.algvar_t = ch_g.edge.probability_of_any_speciation + (ch_g.edge.probability_of_no_speciation * ch_g.algvar_x)
                    ch_f.algvar_t = ch_f.edge.probability_of_no_speciation * ch_f.algvar_z
                    self.algvar_z = ch_f.algvar_t * ch_g.algvar_t
            else:
                raise ValueError
        return Z

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
        self.encode_bipartitions()

    def probability_of_good_species_clade(self, taxa):
        ### Prep
        if not isinstance(taxa, set):
            taxa = set(taxa)
        ### Initialization
        Z = 0.0
        # for leaf in self.leaf_node_iter():
        #     leaf.algvar_x = 0.0
        #     leaf.algvar_y = 0.0
        #     if leaf.taxon in taxa:
        #         leaf.algvar_s = 1
        #     else
        #         leaf.algvar_s = 2
        ### Sweep
        for ndi in self.postorder_node_iter():
            Z += ndi.calc_probability_of_good_species_clade(taxa)
        ### Finalization
        if self.seed_node.algvar_s == 1:
            Z = self.seed_node.algvar_y
        elif self.seed_node.algvar_s == 3:
            Z += self.seed_node.algvar_z
        else:
            raise ValueError
        return Z

